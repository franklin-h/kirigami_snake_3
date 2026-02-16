% plot_envelope_fixedcols.m
% Assumes CSV columns: [Time, Displacement(mm), Force(N)] with no name lookups.
% New:
%   - Mean force–displacement curve across trials
%   - Shaded ±1 SE envelope around the mean curve
% Existing:
%   - Choose representative file by index
%   - Stress–strain (strain = disp/L0; stress = Force/Area -> MPa)
%   - LaTeX axis/tick labels; no titles

clear; clc; close all;

% --- Settings --------------------------------------------------------------
L0 = 65;               % initial length in mm for strain
filePattern = '*.csv'; % file pattern to gather trials
showAllTrialsLight = false; % set true to draw all trials faintly
nGrid = 600;           % number of points for common displacement grid
% --------------------------------------------------------------------------

% Folder selection (cancel => current folder)
dataFolder = uigetdir(pwd, 'Select folder containing CSV trial files');
if isequal(dataFolder,0), dataFolder = pwd; end

files = dir(fullfile(dataFolder, filePattern));
if isempty(files)
    error('No CSV files matching "%s" found in %s', filePattern, dataFolder);
end

% Stable order
[~, order] = sort({files.name});
files = files(order);

fprintf('\nFound %d CSV files in "%s":\n', numel(files), dataFolder);
for k = 1:numel(files), fprintf('%2d) %s\n', k, files(k).name); end

repIdx = input('\nEnter the number of the representative file to plot: ');
if isempty(repIdx) || ~isscalar(repIdx) || repIdx < 1 || repIdx > numel(files)
    error('Representative index must be an integer between 1 and %d.', numel(files));
end

A_mm2 = input('Enter cross-sectional area (in mm^2) for stress calc (required): ');
if isempty(A_mm2) || ~isscalar(A_mm2) || A_mm2 <= 0
    error('Cross-sectional area must be a positive scalar (mm^2).');
end

% Read all trials (first 3 columns only)
n = numel(files);
allT = cell(n,1);  % time
allU = cell(n,1);  % displacement (mm)
allF = cell(n,1);  % force (N)

for k = 1:n
    fpath = fullfile(dataFolder, files(k).name);
    [t,u,f] = readTrialCSV_fixedCols(fpath); % assumes columns [t,u,f]
    allT{k} = t(:);
    allU{k} = u(:);
    allF{k} = f(:);
end

% Common time grid: representative trial's time (kept; not used for mean curve)
tRef = allT{repIdx};
m = numel(tRef);
Umat = nan(m, n);
Fmat = nan(m, n);
for k = 1:n
    Umat(:,k) = interp1(allT{k}, allU{k}, tRef, 'linear', NaN);
    Fmat(:,k) = interp1(allT{k}, allF{k}, tRef, 'linear', NaN);
end

% Stats vs time (not used for mean F–u curve, but kept)
U_mu   = mean(Umat, 2, 'omitnan'); U_sd = std(Umat, 0, 2, 'omitnan');
U_min  = min (Umat, [], 2, 'omitnan'); U_max = max(Umat, [], 2, 'omitnan');
F_mu_t = mean(Fmat, 2, 'omitnan'); F_sd_t = std(Fmat, 0, 2, 'omitnan'); %#ok<NASGU>
F_min  = min (Fmat, [], 2, 'omitnan'); F_max = max(Fmat, [], 2, 'omitnan');

U_rep = Umat(:,repIdx);
F_rep = Fmat(:,repIdx);

% Stress–strain per trial (kept)
strain = cell(n,1);
stress = cell(n,1);
for k = 1:n
    epsK = allU{k} ./ L0;        % unitless
    sigK = allF{k} ./ A_mm2;     % MPa (N/mm^2)
    msk  = isfinite(epsK) & isfinite(sigK);
    strain{k} = epsK(msk);       % keep acquisition order (time order)
    stress{k} = sigK(msk);
end

% ===================== NEW: Mean F–u curve with ±1 SE =====================
% Build a common displacement grid over the overlapping range across trials.
uMin = -inf; uMax = inf;
U_sorted = cell(n,1); F_sorted = cell(n,1);
for k = 1:n
    Uk = allU{k}; Fk = allF{k};
    msk = isfinite(Uk) & isfinite(Fk);
    Uk = Uk(msk); Fk = Fk(msk);
    % Sort by displacement and drop duplicate U to satisfy interp1
    [Uk, idx] = sort(Uk);
    Fk = Fk(idx);
    [Uk, iu] = unique(Uk, 'stable');
    Fk = Fk(iu);

    if numel(Uk) >= 2
        U_sorted{k} = Uk;
        F_sorted{k} = Fk;
        uMin = max(uMin, Uk(1));
        uMax = min(uMax, Uk(end));
    end
end

if ~isfinite(uMin) || ~isfinite(uMax) || uMax <= uMin
    error('Displacement ranges do not overlap across trials; cannot compute mean curve.');
end

uGrid = linspace(uMin, uMax, nGrid).';           % common displacement vector
Fgrid = nan(numel(uGrid), n);
for k = 1:n
    if ~isempty(U_sorted{k})
        Fgrid(:,k) = interp1(U_sorted{k}, F_sorted{k}, uGrid, 'linear', NaN);
    end
end

F_mu = mean(Fgrid, 2, 'omitnan');                % pointwise mean force
F_sd = std (Fgrid, 0, 2, 'omitnan');             % sample SD (N-1)
N_eff = sum(isfinite(Fgrid), 2);                  % available trials per point
% F_se = F_sd ./ sqrt(N_eff);                       % standard error
F_se = F_sd; 
mask = (N_eff >= 2) & isfinite(F_mu) & isfinite(F_se);

% ===================== Plotting ==========================================
stress_strain = figure('Color','w'); hold on; box on;
ax = gca;
ax.FontName = 'Helvetica';
ax.FontSize = 25;
ax.XLabel.FontSize=30;
ax.YLabel.FontSize=30; 
ax.TickLabelInterpreter = 'none';

% --- 6-color order: 3 from 'dye', then 3 from default 'lines' ------------
try
    cDye = dye(3);
catch
    axTmp = axes('Visible','off');
    colororder(axTmp,'dye');
    cDye = colororder(axTmp); 
    cDye = cDye(1:3,:);
    delete(axTmp);
end
cLines = lines(3);
ax.ColorOrder      = [cDye; cLines];
ax.ColorOrderIndex = 1;

% Labels for legend
labels = {'Plastic Sample 1','Plastic Sample 2','Plastic Sample 3', ...
          'Wood Sample 1','Wood Sample 2','Wood Sample 3'};

% Plot individual trials
for k = 1:n
    if k <= numel(labels), dn = labels{k}; else, dn = sprintf('Sample %d', k); end
    lw = showAllTrialsLight * 1 + (~showAllTrialsLight) * 4;
    % plot(allU{k}, allF{k}, '-', 'LineWidth', lw, ...
    %      'DisplayName', dn, 'HandleVisibility','on');
end

% --- NEW: SE band (plot first), then mean curve ---------------------------
if any(mask)
    uBand  = uGrid(mask);
    muBand = F_mu(mask);
    seBand = F_se(mask);

    % Shaded ±1 SE region
    fill([uBand; flipud(uBand)], ...
         [muBand + seBand; flipud(muBand - seBand)], ...
         [0 0 0], 'FaceAlpha', 0.15, 'EdgeColor','none', ...
         'DisplayName', '±1 SE');

    % Mean curve on top
    plot(uGrid, F_mu, 'LineWidth', 2, 'DisplayName','Mean');
else
    warning('Insufficient overlap to compute SE band (need at least 2 trials per grid point).');
    plot(uGrid, F_mu, 'k-', 'LineWidth', 2, 'DisplayName','Mean F');
end

% Legend
leg = legend('Location','north', 'Interpreter','none');
leg.FontName = 'Helvetica';

ylim([0, 10]);
xlim([0, 44]);
xlabel("Displacement (mm)");
ylabel("Force (N)");

% Save figure
outPng = fullfile(dataFolder, 'force_displacement_mean_SE.png');
exportgraphics(stress_strain, outPng, 'Resolution', 400);
fprintf('Saved figure: %s\n', outPng);

% ========================= Helper Function ================================
function [t,u,f] = readTrialCSV_fixedCols(fpath)
% Reads CSV and returns [time, displacement(mm), force(N)] from columns 1:3.
% No header/units parsing; robust to header rows since readmatrix usually
% skips them for numeric output. Cleans NaN rows and sorts by time.

    M = readmatrix(fpath);             % returns numeric matrix
    if size(M,2) < 3
        error('File "%s" must have at least three columns [Time, Disp, Force].', fpath);
    end
    M = M(:,1:3);
    % Drop rows with non-finite values in the first 3 columns
    M = M(all(isfinite(M),2), :);
    if isempty(M)
        error('No numeric data found in the first three columns of "%s".', fpath);
    end
    t = M(:,1); u = M(:,2); f = M(:,3);

    % Sort by time (ascending) and remove duplicate timestamps
    [t, isort] = sort(t);
    u = u(isort); f = f(isort);
    [t, iu] = unique(t, 'stable');
    u = u(iu); f = f(iu);
end
