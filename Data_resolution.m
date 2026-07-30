clear; close all; clc;

%% ============================================================
% SEPARATE proxy-resolution plots
%
% - one figure per proxy record
% - main y-range: 0–1.5 kyr
% - broken y-axis only where needed
% - 250 yr reference line = red dashed
% - 500 yr reference line = yellow/orange dashed
% - x-axis shown as index / age
% - all figures are generated sequentially, remain open, and are not saved
%% ============================================================

%% =========================
% Global font settings
%% =========================
set(groot, ...
    'defaultAxesFontName', 'Cambria Math', ...
    'defaultTextFontName', 'Cambria Math', ...
    'defaultLegendFontName', 'Cambria Math');

%% =========================
% Folder setup
%% =========================
rootDir = pwd;
inputDir = fullfile(rootDir, 'input');

if exist(inputDir, 'dir')
    rawFiles = dir(fullfile(inputDir, '*.txt'));
else
    inputDir = rootDir;
    rawFiles = dir(fullfile(inputDir, '*.txt'));
end

if isempty(rawFiles)
    error(['No .txt files found. Put the proxy files into an input ', ...
           'folder or next to this script.']);
end

fprintf('Found %d .txt files.\n', numel(rawFiles));

%% =========================
% Style settings
%% =========================
plotBackground = [0.975 0.935 0.945];
proxyLineColor = [0 0 0];

% Reference lines
line250Color = [0.80 0.10 0.10];
line500Color = [1.00 0.55 0.00];
referenceLineStyle = '--';
referenceLineW = 1.1;

% Axis and separator colours
separatorColor = [0.35 0.35 0.35];
axisColorY = [0.35 0.35 0.35];
axisColorX = [0.35 0.35 0.35];
tickLabelColor = [0.35 0.35 0.35];

% Line widths
proxyLineWidth = 0.5;
axisLineW = 0.45;
separatorLineW = 0.25;
breakMarkLineW = 0.45;

% Font sizes
fontSizeAxesMain = 11;
fontSizeAxesSmall = 9;
fontSizeTitle = 13;
fontSizeLabel = 12;

% Grid
gridAlpha = 0.22;
gridColor = [1 1 1];

% Y-axis ranges
baseYLim = [0 1.5];
midYLim = [1.5 5];

% Reference values in kyr
line250 = 0.25;
line500 = 0.50;

nXTicks = 5;

%% =========================
% Process each file separately
%% =========================
for f = 1:numel(rawFiles)

    filePath = fullfile(rawFiles(f).folder, rawFiles(f).name);
    [~, fileBaseName, ~] = fileparts(filePath);

    %% -------------------------
    % Clean proxy display names
    %% -------------------------
    switch fileBaseName

        case {'EPICA_deuterium', 'EPICA'}
            displayName = 'EPICA';

        case {'Fuji_deuterium', 'Dome_Fuji', 'Dome Fuji', 'Fuji'}
            displayName = 'Dome Fuji';

        case {'Speleothem_D18O', 'Speleothem'}
            displayName = 'Speleothem';

        case {'1094_N. pachydermaD18O', 'ODP1094', ...
              'ODP 1094', '1094'}
            displayName = 'ODP 1094';

        case {'1089_G.bulloidesD18O', 'ODP1089', ...
              'ODP 1089', '1089'}
            displayName = 'ODP 1089';

        case {'983_Nps%', 'ODP983', 'ODP 983', '983'}
            displayName = 'ODP 983';

        case {'PS75_076_lnZrRb', 'PS75_076', ...
              'PS75-076', 'PS75/076'}
            displayName = 'PS75/076';

        case {'U1385_logCaTi', 'U1385'}
            displayName = 'U1385';

        case {'U1308_logSiSr', 'U1308'}
            displayName = 'U1308';

        otherwise
            displayName = strrep(fileBaseName, '_', ' ');
    end

    fprintf('Reading %s\n', rawFiles(f).name);

    %% -------------------------
    % Read data
    %% -------------------------
    try
        data = readmatrix(filePath);
    catch
        warning('Skipping %s because it could not be read.', ...
            rawFiles(f).name);
        continue;
    end

    if isempty(data) || size(data, 2) < 2
        warning('Skipping %s because it has fewer than two columns.', ...
            rawFiles(f).name);
        continue;
    end

    data = data(:, 1:2);
    data = data(all(isfinite(data), 2), :);

    if size(data, 1) < 3
        warning('Skipping %s because it has too few valid rows.', ...
            rawFiles(f).name);
        continue;
    end

    %% -------------------------
    % Sort and clean age values
    %% -------------------------
    age = data(:, 1);
    age = sort(age);

    age = unique(age, 'stable');

    if numel(age) < 3
        warning(['Skipping %s because too few unique age values ', ...
                 'remain.'], rawFiles(f).name);
        continue;
    end

    %% -------------------------
    % Calculate temporal spacing
    %% -------------------------
    dAge = diff(age);

    intervalAge = ...
        (age(1:end-1) + age(2:end)) ./ 2;

    x = (1:numel(dAge))';

    good = isfinite(dAge) & dAge >= 0;

    dAge = dAge(good);
    intervalAge = intervalAge(good);
    x = x(good);

    if isempty(dAge)
        warning('Skipping %s because no valid age differences remain.', ...
            rawFiles(f).name);
        continue;
    end

    maxGap = max(dAge);

    %% -------------------------
    % Determine axis configuration
    %% -------------------------
    if maxGap <= baseYLim(2)
        plotType = 1;

    elseif maxGap <= midYLim(2)
        plotType = 2;

    else
        plotType = 3;
    end

    fprintf('Plotting %s, maximum gap %.2f kyr\n', ...
        displayName, maxGap);

    %% =========================
    % Create separate figure
    %% =========================
    fig = figure( ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [100 100 900 520]);

    tilePos = [0.12 0.14 0.83 0.78];

    axLow = [];
    axMid = [];
    axTop = [];

    %% -------------------------
    % Create axes
    %% -------------------------
    switch plotType

        case 1
            axLow = axes(fig, ...
                'Position', tilePos);

        case 2
            lowFrac = 0.79;
            overlap = 0.0035 * tilePos(4);

            axLow = axes(fig, ...
                'Position', ...
                [tilePos(1), ...
                 tilePos(2), ...
                 tilePos(3), ...
                 tilePos(4) * lowFrac + overlap]);

            axTop = axes(fig, ...
                'Position', ...
                [tilePos(1), ...
                 tilePos(2) + tilePos(4) * lowFrac - overlap, ...
                 tilePos(3), ...
                 tilePos(4) * (1 - lowFrac) + overlap]);

        case 3
            lowFrac = 0.65;
            midFrac = 0.19;
            topFrac = 0.16;
            overlap = 0.0035 * tilePos(4);

            axLow = axes(fig, ...
                'Position', ...
                [tilePos(1), ...
                 tilePos(2), ...
                 tilePos(3), ...
                 tilePos(4) * lowFrac + overlap]);

            axMid = axes(fig, ...
                'Position', ...
                [tilePos(1), ...
                 tilePos(2) + tilePos(4) * lowFrac - overlap, ...
                 tilePos(3), ...
                 tilePos(4) * midFrac + 2 * overlap]);

            axTop = axes(fig, ...
                'Position', ...
                [tilePos(1), ...
                 tilePos(2) + ...
                 tilePos(4) * (lowFrac + midFrac) - overlap, ...
                 tilePos(3), ...
                 tilePos(4) * topFrac + overlap]);
    end

    axList = axLow;

    if ~isempty(axMid)
        axList = [axList axMid];
    end

    if ~isempty(axTop)
        axList = [axList axTop];
    end

    %% -------------------------
    % Plot data on all axes
    %% -------------------------
    for a = 1:numel(axList)

        ax = axList(a);
        hold(ax, 'on');

        plot(ax, x, dAge, ...
            'Color', proxyLineColor, ...
            'LineWidth', proxyLineWidth);

        ax.Color = plotBackground;
        ax.Box = 'off';
        ax.Layer = 'top';
        ax.LineWidth = axisLineW;
        ax.FontName = 'Cambria Math';

        ax.XColor = axisColorX;
        ax.YColor = axisColorY;

        grid(ax, 'on');

        ax.GridColor = gridColor;
        ax.GridAlpha = gridAlpha;

        xlim(ax, [min(x) max(x)]);
    end

    set(axLow, ...
        'FontSize', fontSizeAxesMain);

    if ~isempty(axMid)
        set(axMid, ...
            'FontSize', fontSizeAxesSmall);
    end

    if ~isempty(axTop)
        set(axTop, ...
            'FontSize', fontSizeAxesSmall);
    end

    %% -------------------------
    % Y-axis limits and ticks
    %% -------------------------
    ylim(axLow, baseYLim);
    yticks(axLow, [0 0.25 0.5 1.0 1.5]);

    if plotType == 2

        upperMax = ceil(maxGap * 1.05 * 10) / 10;

        if upperMax <= 1.5
            upperMax = 2;
        end

        ylim(axTop, [1.5 upperMax]);

        topTicks = makeUpperTicks(1.5, upperMax, 3);
        yticks(axTop, topTicks);

    elseif plotType == 3

        ylim(axMid, midYLim);
        yticks(axMid, 3);

        upperMax = ceil(maxGap * 1.05);

        if upperMax <= 5
            upperMax = 6;
        end

        ylim(axTop, [5 upperMax]);

        topTicks = makeUpperTicks(5, upperMax, 3);
        yticks(axTop, topTicks);
    end

    %% -------------------------
    % Reference lines
    %% -------------------------
    yline(axLow, line250, referenceLineStyle, ...
        'Color', line250Color, ...
        'LineWidth', referenceLineW, ...
        'HandleVisibility', 'off');

    yline(axLow, line500, referenceLineStyle, ...
        'Color', line500Color, ...
        'LineWidth', referenceLineW, ...
        'HandleVisibility', 'off');

    %% -------------------------
    % Separator lines
    %% -------------------------
    if plotType == 2

        drawSeparatorLine( ...
            axLow, ...
            baseYLim(2), ...
            separatorColor, ...
            separatorLineW);

    elseif plotType == 3

        drawSeparatorLine( ...
            axLow, ...
            baseYLim(2), ...
            separatorColor, ...
            separatorLineW);

        drawSeparatorLine( ...
            axMid, ...
            midYLim(2), ...
            separatorColor, ...
            separatorLineW);
    end

    %% -------------------------
    % X-axis ticks and labels
    %% -------------------------
    xTicks = round(linspace(min(x), max(x), nXTicks));
    xTicks = unique(xTicks);

    xTickLabels = cell(size(xTicks));

    for i = 1:numel(xTicks)

        [~, idxClosest] = min(abs(x - xTicks(i)));

        xTicks(i) = x(idxClosest);

        xTickLabels{i} = sprintf( ...
            '%d / %.0f', ...
            x(idxClosest), ...
            intervalAge(idxClosest));
    end

    set(axLow, ...
        'XTick', xTicks, ...
        'XTickLabel', xTickLabels);

    xlabel(axLow, ...
        'Index / Age [kyr]', ...
        'Color', tickLabelColor, ...
        'FontName', 'Cambria Math', ...
        'FontSize', fontSizeLabel);

    ylabel(axLow, ...
        '\Delta age [kyr]', ...
        'Color', tickLabelColor, ...
        'FontName', 'Cambria Math', ...
        'FontSize', fontSizeLabel);

    if plotType == 2

        set(axTop, ...
            'XTick', [], ...
            'XTickLabel', []);

    elseif plotType == 3

        set(axMid, ...
            'XTick', [], ...
            'XTickLabel', []);

        set(axTop, ...
            'XTick', [], ...
            'XTickLabel', []);
    end

    %% -------------------------
    % Figure title
    %% -------------------------
    if plotType == 1

        title(axLow, displayName, ...
            'Interpreter', 'none', ...
            'FontName', 'Cambria Math', ...
            'FontSize', fontSizeTitle, ...
            'FontWeight', 'bold');

    else

        title(axTop, displayName, ...
            'Interpreter', 'none', ...
            'FontName', 'Cambria Math', ...
            'FontSize', fontSizeTitle, ...
            'FontWeight', 'bold');
    end

    %% -------------------------
    % Axis-break marks
    %% -------------------------
    if plotType == 2

        addYAxisBreakMark( ...
            fig, axLow, 'top', ...
            separatorColor, breakMarkLineW);

        addYAxisBreakMark( ...
            fig, axTop, 'bottom', ...
            separatorColor, breakMarkLineW);

    elseif plotType == 3

        addYAxisBreakMark( ...
            fig, axLow, 'top', ...
            separatorColor, breakMarkLineW);

        addYAxisBreakMark( ...
            fig, axMid, 'bottom', ...
            separatorColor, breakMarkLineW);

        addYAxisBreakMark( ...
            fig, axMid, 'top', ...
            separatorColor, breakMarkLineW);

        addYAxisBreakMark( ...
            fig, axTop, 'bottom', ...
            separatorColor, breakMarkLineW);
    end

    %% -------------------------
    % Display figure and continue
    %% -------------------------
    drawnow;

    fprintf('Displayed %s.\n', displayName);
end

fprintf('All separate figures have been displayed and remain open.\n');

%% ============================================================
% Local functions
%% ============================================================

function ticks = makeUpperTicks(yMin, yMax, maxN)
% Generate a limited number of ticks above the break.

    if yMax <= yMin
        ticks = [];
        return;
    end

    rangeVal = yMax - yMin;

    if rangeVal <= 2
        step = 0.5;

    elseif rangeVal <= 5
        step = 1;

    elseif rangeVal <= 12
        step = 2;

    else
        step = 5;
    end

    ticks = ...
        (ceil(yMin / step) * step):step:yMax;

    ticks = ticks(ticks > yMin);

    if numel(ticks) > maxN
        idx = round(linspace(1, numel(ticks), maxN));
        ticks = ticks(idx);
    end

    ticks = unique(round(ticks, 2));

    if isempty(ticks)
        ticks = round(yMax, 2);
    end
end

function drawSeparatorLine(ax, yVal, lineColor, lineW)

    xl = xlim(ax);

    line(ax, xl, [yVal yVal], ...
        'Color', lineColor, ...
        'LineWidth', lineW, ...
        'Clipping', 'off');
end

function addYAxisBreakMark( ...
    fig, ax, location, lineColor, lineW)

    oldUnits = ax.Units;
    ax.Units = 'normalized';

    pos = ax.Position;

    ax.Units = oldUnits;

    x0 = pos(1) - 0.0015;

    switch lower(location)

        case 'top'
            y0 = pos(2) + pos(4);

        case 'bottom'
            y0 = pos(2);

        otherwise
            error('Location must be top or bottom.');
    end

    dx = 0.0032;
    dy = 0.0025;
    gap = 0.0024;

    annotation(fig, 'line', ...
        [x0 - dx, x0 + dx], ...
        [y0 - dy, y0 + dy], ...
        'Color', lineColor, ...
        'LineWidth', lineW);

    annotation(fig, 'line', ...
        [x0 - dx, x0 + dx], ...
        [y0 - dy - gap, y0 + dy - gap], ...
        'Color', lineColor, ...
        'LineWidth', lineW);
end