% plot_envelope_fixedcols.m
% Assumes CSV columns: [Time, Displacement(mm), Force(N)] with no name lookups.
% Features:
%   - Choose representative file by index
%   - ±1 SD shaded band, min/max lines, representative curve
%   - Stress–strain (strain = disp/L0; stress = Force/Area -> MPa)
%   - LaTeX axis/tick labels; no titles

clear; clc; close all;

% --- Settings --------------------------------------------------------------
L0 = 65;               % initial length in mm for strain
filePattern = '*.csv'; % file pattern to gather trials
showAllTrialsLight = false; % set true to draw all trials faintly
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

% Common time grid: representative trial's time
tRef = allT{repIdx};
m = numel(tRef);
Umat = nan(m, n);
Fmat = nan(m, n);

for k = 1:n
    Umat(:,k) = interp1(allT{k}, allU{k}, tRef, 'linear', NaN);
    Fmat(:,k) = interp1(allT{k}, allF{k}, tRef, 'linear', NaN);
end

% Stats vs time
U_mu   = mean(Umat, 2, 'omitnan'); U_sd = std(Umat, 0, 2, 'omitnan');
U_min  = min (Umat, [], 2, 'omitnan'); U_max = max(Umat, [], 2, 'omitnan');

F_mu   = mean(Fmat, 2, 'omitnan'); F_sd = std(Fmat, 0, 2, 'omitnan');
F_min  = min (Fmat, [], 2, 'omitnan'); F_max = max(Fmat, [], 2, 'omitnan');

U_rep = Umat(:,repIdx);
F_rep = Fmat(:,repIdx);

strain = cell(n,1);
stress = cell(n,1);
for k = 1:n
    epsK = allU{k} ./ L0;        % unitless
    sigK = allF{k} ./ A_mm2;     % MPa (N/mm^2)
    msk  = isfinite(epsK) & isfinite(sigK);
    strain{k} = epsK(msk);       % keep acquisition order (time order)
    stress{k} = sigK(msk);
end

stress_strain = figure('Color','w'); hold on; box on;
colororder("dye")
ax = gca;
ax.FontName = 'Helvetica';              % <- use Helvetica
ax.FontSize = 25;
ax.XLabel.FontSize=30;
ax.YLabel.FontSize=30; 
ax.TickLabelInterpreter = 'none';   % <- NOT 'latex'


% --- after creating the axes (right after: ax = gca; ...) ---------------
% Build a 6-color order: 3 from 'dye', then 3 from default 'lines'
try
    cDye = dye(3);                  % works in newer MATLAB versions
catch
    % Fallback way to fetch 'dye' colors if dye() isn't available
    axTmp = axes('Visible','off');
    colororder(axTmp,'dye');
    cDye = colororder(axTmp); 
    cDye = cDye(1:3,:);
    delete(axTmp);
end
cLines = lines(3);
ax.ColorOrder      = [cDye; cLines];
ax.ColorOrderIndex = 1;             % start at the first color

% ------------------------------------------------------------------------

% Use explicit legend labels (will also control DisplayName)
labels = {'Plastic Sample 1','Plastic Sample 2','Plastic Sample 3', ...
          'Wood Sample 1','Wood Sample 2','Wood Sample 3'};

% Plot every trial (color comes from ColorOrder automatically)
for k = 1:n
    if k <= numel(labels)
        dn = labels{k};
    else
        dn = sprintf('Sample %d', k); % fallback if more than 6 files
    end
    plot(allU{k}, allF{k}, '-', 'LineWidth', 4, ...
         'DisplayName', dn, 'HandleVisibility','on');
end

% Legend (unchanged location, just picks up the new names)
leg = legend('Location','north', 'Interpreter','none');
leg.FontName = 'Helvetica';


% grid on; 
% legend('Interpreter','latex','Location','north');
ylim([0, 10]);
xlim([0, 44])
xlabel("Displacement (mm)")
ylabel("Force (N)")
% Save figure
% outPng = fullfile(dataFolder, 'stress_strain_all_curves.png');
% legend show 
% legend(Location='north')
title = strcat(dataFolder,".png"); 
exportgraphics(stress_strain, title, 'Resolution', 400);

% fprintf('Saved stress–strain figure to: %s\n', outPng);

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
