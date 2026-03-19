% plot_envelope_fixedcols.m
% Assumes CSV columns: [Time, Displacement(mm), Force(N)] with no name lookups.
% Features:
%   - Choose representative file by index
%   - ±1 SD shaded band, min/max lines, representative curve
%   - Stress–strain (strain = disp/L0; stress = Force/Area -> MPa)
%   - LaTeX axis/tick labels; no titles
%

% 
% MOD: Legend shows line samples (via proxy lines) instead of squares for surfaces.

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
files = files(order(end:-1:1));

fprintf('\nFound %d CSV files in "%s":\n', numel(files), dataFolder);
for k = 1:numel(files), fprintf('%2d) %s\n', k, files(k).name); end

repIdx = input('\nEnter the number of the representative file to plot (Deprecated, enter anything): ');
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

% --- Displacement vs Time --------------------------------------------------
% figure('Color','w'); hold on; box on;
% ax = gca; ax.TickLabelInterpreter = 'latex';
%
% xBand = [tRef; flipud(tRef)];
% yBand = [U_mu - U_sd; flipud(U_mu + U_sd)];
% fill(xBand, yBand, [0.6 0.8 1.0], 'EdgeColor','none', 'FaceAlpha',0.35, ...
%     'DisplayName','Mean $\pm$ 1 SD');
%
% if showAllTrialsLight
%     for k = 1:n
%         plot(tRef, Umat(:,k), 'Color',[0.7 0.7 0.7], 'HandleVisibility','off');
%     end
% end
%
% xlabel('Time', 'Interpreter','latex');
% ylabel('Displacement', 'Interpreter','latex');
% grid on; legend('Interpreter','latex','Location','best');

% --- Force vs Time ---------------------------------------------------------
% figure('Color','w'); hold on; box on;
% ax = gca; ax.TickLabelInterpreter = 'latex';
%
% xBand = [tRef; flipud(tRef)];
% yBand = [F_mu - F_sd; flipud(F_mu + F_sd)];
% fill(xBand, yBand, [0.8 0.9 0.6], 'EdgeColor','none', 'FaceAlpha',0.35, ...
%     'DisplayName','Mean $\pm$ 1 SD');
%
% if showAllTrialsLight
%     for k = 1:n
%         plot(tRef, Fmat(:,k), 'Color',[0.75 0.75 0.75], 'HandleVisibility','off');
%     end
% end
%
% xlabel('Time', 'Interpreter','latex');
% ylabel('Force', 'Interpreter','latex');
% grid on; legend('Interpreter','latex','Location','best');

% --- Stress–Strain (here using Displacement vs Force) with per-curve gradient
% Settings you can tweak:
maxPtsPerCurve = 2000;      % downsample target per curve (set [] to disable)
lightFrac       = 0.1;      % 0..1, how close to white the start is (0=no lighten)
lineWid         = 1;        % line width

% Convert to strain/stress (keep time order) — optional, not used below
strain = cell(n,1); stress = cell(n,1);
for k = 1:n
    epsK = allU{k} ./ L0;
    sigK = allF{k} ./ A_mm2;   % MPa
    msk  = isfinite(epsK) & isfinite(sigK);
    epsK = epsK(msk); sigK = sigK(msk);

    % Optional downsampling for speed
    if ~isempty(maxPtsPerCurve) && numel(epsK) > maxPtsPerCurve
        idx  = round(linspace(1,numel(epsK), maxPtsPerCurve));
        epsK = epsK(idx); sigK = sigK(idx);
    end
    strain{k} = epsK; stress{k} = sigK;
end

stress_strain = figure('Color','w');

set(stress_strain, 'Renderer','opengl', 'GraphicsSmoothing','off'); % speed

% ---- Fonts/interpreters: use Helvetica everywhere, no LaTeX ----
set(stress_strain, 'DefaultTextFontName','Helvetica', ...
                   'DefaultAxesFontName','Helvetica', ...
                   'DefaultTextInterpreter','none', ...
                   'DefaultAxesTickLabelInterpreter','none');

ax = gca; hold(ax,'on'); box(ax,'on');
ax.TickLabelInterpreter = 'none';
ax.FontName = 'Helvetica';
ax.FontSize = 25;
ax.XLabel.FontSize = 30;
ax.YLabel.FontSize = 30;
legend('AutoUpdate','off');  % legend building can be expensive

% colororder("dye")
% --- Force the 3 series colors: pink, red, burgundy ---
colororder(ax, [ ...
    1.0000 0.7529 0.7961;  % pink    (#FFC0CB)
    1.0000 0.0000 0.0000;  % red     (#FF0000)
    0.5020 0.0000 0.1255]); % burgundy (#800020)

% colororder(ax, [ ...
%     0.6784 0.8471 0.9020;  % light blue  (#ADD8E6)
%     0.2549 0.4118 0.8824;  % royal blue  (#4169E1)
%     0.0000 0.0000 0.5020]);% navy        (#000080)
CO  = ax.ColorOrder; nCO = size(CO,1);
% mixWithWhite = @(base,a) 1 - (1-base).*a;   % lighten towards white
mixWithBlack = @(base,a) base .* a;          % darken towards black (a in [0..1])

% ---- Legend labels (customize as needed) ---------------------------------
legend_labels = ["12 mm", "8 mm", "4 mm"];
if numel(legend_labels) ~= n
    % Fallback to file names if counts don't match
    legend_labels = string({files.name});
end

% Handles
hSurf = gobjects(n,1);     % visible gradient "lines" (as surfaces)
hLeg  = gobjects(n,1);     % legend proxy lines (NaN lines)

for k = 1:n
    % Choose data: using displacement & force for axes
    x = allU{k}; y = allF{k};
    % If you want strain/stress instead, swap the previous line:
    % x = strain{k}; y = stress{k};

    if numel(x) < 2, continue; end

    % base = CO(mod(k-1,nCO)+1, :);           % default MATLAB color for this trial
    % c0   = mixWithWhite(base, lightFrac);    % starting (lighter) color
    % t    = linspace(1,0,numel(x))';          % 1 = base (dark) start, 0 = light end
    % cols = c0 + (base - c0).*t;              % dark -> light


    base = CO(mod(k-1,nCO)+1, :);
    
    % Option A (recommended): end near black, controlled by lightFrac
    % lightFrac = 0.1 means endpoint is 10% of base (almost black)
    % c1   = mixWithBlack(base, lightFrac);        % ending (near-black) color
    % t    = linspace(0,1,numel(x))';              % 0=start, 1=end
    % cols = base + (c1 - base).*t;                % base -> near-black

    tLight = 0.9;     % 0..1, how light the START is (0=white, 1=base)
    tBlack = 0.50;     % 0..1, how dark the END is (0=black, 1=base)
    
    cStart = 1 - (1 - base)*tLight;  % lighter tint toward white
    cEnd   = base * tBlack;          % toward black
    
    t    = linspace(0,1,numel(x))';
    cols = cStart + (cEnd - cStart).*t;  % start(light) -> end(black)

    % Build a single SURFACE as a colored polyline (fast, one object/curve)
    X = [x x];
    Y = [y y];
    Z = zeros(numel(x),2);
    C = zeros(numel(x),2,3);    % truecolor (M×N×3)
    C(:,1,:) = cols;
    C(:,2,:) = cols;

    hSurf(k) = surface('XData',X,'YData',Y,'ZData',Z,'CData',C, ...
                       'FaceColor','none','EdgeColor','interp', ...
                       'LineWidth',lineWid, ...
                       'DisplayName',legend_labels(k), ...
                       'HitTest','off','PickableParts','none');

    % Hide the surface from the legend (otherwise it shows a square)
    hSurf(k).Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Add an invisible line as the legend proxy (shows as a line sample)
    hLeg(k) = plot(ax, NaN, NaN, '-', ...
                   'LineWidth', lineWid, 'Color', base, ...
                   'DisplayName', legend_labels(k));
end

xlabel('Displacement (mm)', 'Interpreter','none', 'FontName','Helvetica'); % was 'latex'
ylabel('Force (N)',        'Interpreter','none', 'FontName','Helvetica'); % was 'latex'
grid off;
ylim([-0.3, 0.6]); xlim([0, 12.2]);

% Legend (use proxy lines so legend shows line samples)
valid = isgraphics(hLeg);
leg = legend(hLeg(valid), legend_labels(valid), 'Interpreter','none', 'Location','best');
leg.FontName = 'Helvetica';
leg.NumColumns = min(nnz(valid), 3);
leg.ItemTokenSize = [30, 18];  % optional: lengthen the line token

% Save figure
exportgraphics(stress_strain, strcat(dataFolder,".png"), 'Resolution', 400);

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
