% plot_envelope_fixedcols_noColorChange.m
% Assumes CSV columns: [Time, Displacement(mm), Force(N)] with no name lookups.
% Plots Displacement vs Force for each file with NO per-curve gradient/color change
% (each curve is a single solid color). Legend shows line samples.

clear; clc; close all;

% --- Settings --------------------------------------------------------------
L0 = 65;               % initial length in mm for strain (kept; optional)
filePattern = '*.csv'; % file pattern to gather trials
showAllTrialsLight = false; % kept; optional (not used below)
lineWid = 1;           % line width
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

repIdx = input('\nEnter the number of the representative file to plot: ');
if isempty(repIdx) || ~isscalar(repIdx) || repIdx < 1 || repIdx > numel(files)
    error('Representative index must be an integer between 1 and %d.', numel(files));
end %#ok<NASGU>  % repIdx kept for compatibility; not used below

A_mm2 = input('Enter cross-sectional area (in mm^2) for stress calc (required): ');
if isempty(A_mm2) || ~isscalar(A_mm2) || A_mm2 <= 0
    error('Cross-sectional area must be a positive scalar (mm^2).');
end %#ok<NASGU>  % A_mm2 kept for compatibility; not used below

% Read all trials (first 3 columns only)
n = numel(files);
allU = cell(n,1);  % displacement (mm)
allF = cell(n,1);  % force (N)

for k = 1:n
    fpath = fullfile(dataFolder, files(k).name);
    [~,u,f] = readTrialCSV_fixedCols(fpath); % assumes columns [t,u,f]
    allU{k} = u(:);
    allF{k} = f(:);
end

% ===================== Plotting: NO gradient ==============================
fig = figure('Color','w');
set(fig, 'Renderer','opengl', 'GraphicsSmoothing','off'); % speed

% Fonts/interpreters: use Helvetica everywhere, no LaTeX
set(fig, 'DefaultTextFontName','Helvetica', ...
         'DefaultAxesFontName','Helvetica', ...
         'DefaultTextInterpreter','none', ...
         'DefaultAxesTickLabelInterpreter','none');

ax = gca; hold(ax,'on'); box(ax,'on');
ax.TickLabelInterpreter = 'none';
ax.FontName = 'Helvetica';
ax.FontSize = 25;
ax.XLabel.FontSize = 30;
ax.YLabel.FontSize = 30;

% Use a standard color order (solid colors only)
colororder(ax, 'dye');
CO  = ax.ColorOrder; 
nCO = size(CO,1);

% Legend labels (customize as needed)
legend_labels = ["12 mm", "8 mm", "4 mm"];
if numel(legend_labels) ~= n
    legend_labels = string({files.name});
end

hLine = gobjects(n,1);

for k = 1:n
    x = allU{k};
    y = allF{k};
    if numel(x) < 2, continue; end

    base = CO(mod(k-1,nCO)+1, :); % one constant color for this curve

    hLine(k) = plot(ax, x, y, '-', ...
                    'Color', base, ...
                    'LineWidth', lineWid, ...
                    'DisplayName', legend_labels(k));
end

xlabel('Displacement (mm)', 'Interpreter','none', 'FontName','Helvetica');
ylabel('Force (N)',        'Interpreter','none', 'FontName','Helvetica');
grid off;

ylim([-0.3, 0.6]);
xlim([0, 12.2]);

% Legend (line samples automatically)
valid = isgraphics(hLine);
leg = legend(hLine(valid), legend_labels(valid), ...
             'Interpreter','none', 'Location','best');
leg.FontName = 'Helvetica';
leg.NumColumns = min(nnz(valid), 3);
leg.ItemTokenSize = [30, 18];

% Save figure
outPng = fullfile(dataFolder, 'force_displacement_noColorChange.png');
exportgraphics(fig, strcat(dataFolder,".png"), 'Resolution', 400);
fprintf('Saved figure: %s\n', outPng);

% ========================= Helper Function ================================
function [t,u,f] = readTrialCSV_fixedCols(fpath)
% Reads CSV and returns [time, displacement(mm), force(N)] from columns 1:3.
% No header/units parsing; robust to header rows since readmatrix usually
% skips them for numeric output. Cleans NaN rows and sorts by time.

    M = readmatrix(fpath);
    if size(M,2) < 3
        error('File "%s" must have at least three columns [Time, Disp, Force].', fpath);
    end
    M = M(:,1:3);
    M = M(all(isfinite(M),2), :);
    if isempty(M)
        error('No numeric data found in the first three columns of "%s".', fpath);
    end
    t = M(:,1); u = M(:,2); f = M(:,3);

    [t, isort] = sort(t);
    u = u(isort); f = f(isort);
    [t, iu] = unique(t, 'stable');
    u = u(iu); f = f(iu);
end
