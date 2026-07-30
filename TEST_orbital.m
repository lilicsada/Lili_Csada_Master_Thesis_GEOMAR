%% ============================================================
% CURVED-TREND ORBITAL-FORCING SENSITIVITY
%
% ONE TREND ONLY:
% Flat before 600 kyr BP, followed by a gently curved,
% accelerating increase towards the present.
%
% Fixed background:
%   signal:noise = 1:0.2
%   AR(1) phi    = 0.93004 (proxy-like)
%
% Four orbital cases:
%   1. Weak orbital forcing:      signal:orbital = 1:0.2
%   2. Moderate orbital forcing:  signal:orbital = 1:0.5
%   3. Strong orbital forcing:    signal:orbital = 1:1.0
%   4. Phase-shifted moderate orbital forcing:
%      signal:orbital = 1:0.5, 50-kyr temporal shift
%
% Methods:
%   - SiZer
%   - findchangepts
%   - BEAST
%

clear;
close all;
clc;

rng(1);


% 1. BASIC SETTINGS


dt = 0.5;                     % kyr
age_max = 800;                % kyr BP
age = (0:dt:age_max)';        % internal order: 0 -> 800 kyr BP
n = numel(age);

change_age = 600;             % true trend onset, kyr BP

% Proxy-derived parameters
sigma_proxy = 0.88047;
phi_proxy   = 0.93004;

% Total imposed long-term trend change
deltaT = 2 * sigma_proxy;

% Fixed background noise:
% DeltaT : 2*sigma_noise = 1 : 0.2
noise_ratio = 0.2;
sigma_noise = noise_ratio * deltaT / 2;

% Monte Carlo settings
nRep = 1000 ;
burnin = 2000;


% 2. FOUR ORBITAL CASES


orbital_ratios =[0.2 0.5 1.0 1.0];
phase_shift_kyr = [0 0 0 70];
phase_deg_100k = 360 * phase_shift_kyr / 100;

case_names = {
    'Weak orbital: 1:0.2'
    'Moderate orbital: 1:0.5'
    'Strong orbital: 1:1.0'
    'Moderate 1:0.5, shifted 50 kyr'
    };

panel_labels = {
    '1:0.2'
    '1:0.5'
    '1:1.0'
    '1:0.5 + 50 kyr'
    };

row_letters = {'(a)','(b)','(c)','(d)'};

nCase = numel(orbital_ratios);


% 3. ORBITAL PERIODS AND RELATIVE AMPLITUDES


periods = [100 41 23];
relative_amp = [1.0 0.7 0.5];


% 4. CLEAN CURVED TREND


clean_trend = zeros(n,1);
idx_post = age <= change_age;
u = (change_age - age(idx_post)) / change_age;
clean_trend(idx_post) = deltaT * (1 - cos((pi/2) * u));


% 5. SIZER AND FINDCHANGEPTS SETTINGS


opts = struct();
opts.methods = ["sizer" "findchangepts"];

% SiZer
opts.sizer.bandwidth_min = 60;
opts.sizer.bandwidth_max = 100;
opts.sizer.n_bandwidths = 41;
opts.sizer.x_step = 1.0;
opts.sizer.alpha = 0.10;
opts.sizer.simflag = 2;
opts.sizer.iregtdist = 0;
opts.sizer.endpoint_flag = 0;
opts.sizer.boundary_adjustment = 0;
opts.sizer.onset_mode = "all";
opts.sizer.min_run_kyr = 40;
opts.sizer.edge_ignore_kyr = 0;     % no edge clipping
opts.sizer.cluster_gap_kyr = 35;
opts.sizer.min_bw_fraction = 0.80;
opts.sizer.uncertainty_mode = "iqr";

% findchangepts
opts.findchangepts.smooth_kyr = 0;  % no 50-kyr smoothing
opts.findchangepts.max_num_changes = 1;


% 6. BEAST SETTINGS


if exist('beast','file') == 0
    error('BEAST was not found on the MATLAB path.');
end

tcp_minmax = [0 1];
torder_minmax = [0 1];
mcmc_samples = 4000;
mcmc_burnin = 500;
mcmc_chains = 2;
mcmc_thin = 5;


% 7. GENERATE THE SAME AR(1) NOISE REALISATIONS


base_noise = nan(n,nRep);

for r = 1:nRep
    noise_long = zeros(n + burnin, 1);
    for k = 2:(n + burnin)
        noise_long(k) = phi_proxy * noise_long(k-1) + randn;
    end

    noise_raw = noise_long(burnin+1:end);
    base_noise(:,r) = (noise_raw - mean(noise_raw,'omitnan')) ./ std(noise_raw,'omitnan');
end


% 8. BUILD THE FOUR ORBITAL COMPONENTS


orbital_component = nan(n,nCase);
sigma_orbital_case = nan(1,nCase);

for q = 1:nCase
    orbital_ratio = orbital_ratios(q);
    sigma_orbital_case(q) = orbital_ratio * deltaT / 2;
    shift = phase_shift_kyr(q);

    orbital_raw = relative_amp(1) * sin(2*pi*(age + shift) / periods(1)) + ...
                  relative_amp(2) * sin(2*pi*(age + shift) / periods(2)) + ...
                  relative_amp(3) * sin(2*pi*(age + shift) / periods(3));

    orbital_standardised = (orbital_raw - mean(orbital_raw,'omitnan')) ./ std(orbital_raw,'omitnan');
    orbital_component(:,q) = sigma_orbital_case(q) * orbital_standardised;
end


% 9. STORAGE


cp_sizer = nan(nRep,nCase);
cp_find  = nan(nRep,nCase);
cp_beast = nan(nRep,nCase);
beast_ncp_mode = nan(nRep,nCase);
example_series = nan(n,nCase);


% 10. MONTE CARLO LOOP


fprintf('\n');
fprintf('CURVED-TREND ORBITAL-FORCING SENSITIVITY\n');
fprintf('============================================================\n');
fprintf('True onset = %.1f kyr BP\n',change_age);
fprintf('Signal:noise = 1:%.1f\n',noise_ratio);
fprintf('Proxy-like AR(1) phi = %.5f\n',phi_proxy);
fprintf('SiZer edge clipping = NONE\n');
fprintf('findchangepts smoothing = NONE\n');
fprintf('Realisations = %d per orbital case\n',nRep);
fprintf('============================================================\n');

for q = 1:nCase
    fprintf('\nCase %d/%d: %s\n', q, nCase, case_names{q});
    fprintf('  signal:orbital = 1:%.1f\n', orbital_ratios(q));
    fprintf('  temporal phase shift = %.1f kyr\n', phase_shift_kyr(q));

    for r = 1:nRep
        noise = sigma_noise * base_noise(:,r);
        synthetic_raw = clean_trend + orbital_component(:,q) + noise;
        synthetic_z = (synthetic_raw - mean(synthetic_raw,'omitnan')) ./ std(synthetic_raw,'omitnan');

        if r == 1
            example_series(:,q) = synthetic_raw;
        end

        % SiZer and findchangepts
        try
            out = chpttests(age, synthetic_z, opts);
            cp_sizer(r,q) = out.sizer.cp_age;
            cp_find(r,q)  = out.findchangepts.cp_age;
        catch ME
            if r == 1
                fprintf('  SiZer/findchangepts error: %s\n', ME.message);
            end
        end

        % BEAST (no smoothing)
        try
            seed = (q-1) * nRep + r;
            [cp_beast(r,q), beast_ncp_mode(r,q)] = run_beast_once( ...
                synthetic_z, age, dt, age_max, tcp_minmax, torder_minmax, ...
                seed, mcmc_samples, mcmc_burnin, mcmc_chains, mcmc_thin);
        catch ME
            if r == 1
                fprintf('  BEAST error: %s\n', ME.message);
            end
        end

        if mod(r,10) == 0
            fprintf('  %d / %d\n', r, nRep);
        end
    end
end

fprintf('\nAll calculations completed.\n');


% 11. SUMMARY TABLE


method_names = ["SiZer"; "findchangepts"; "BEAST"];

OrbitalRatio = [];
PhaseShift_kyr = [];
Phase_100k_deg = [];
OrbitalCase = strings(0,1);
Method = strings(0,1);
DetectionRate_pct = [];
MedianAge_kyr = [];
Q25_kyr = [];
Q75_kyr = [];
P05_kyr = [];
P95_kyr = [];
MedianError_kyr = [];

for q = 1:nCase
    cp_all = {cp_sizer(:,q); cp_find(:,q); cp_beast(:,q)};

    for m = 1:3
        cp = cp_all{m};
        cp = cp(isfinite(cp));

        OrbitalRatio(end+1,1) = orbital_ratios(q);
        PhaseShift_kyr(end+1,1) = phase_shift_kyr(q);
        Phase_100k_deg(end+1,1) = phase_deg_100k(q);
        OrbitalCase(end+1,1) = string(case_names{q});
        Method(end+1,1) = method_names(m);
        DetectionRate_pct(end+1,1) = 100 * numel(cp) / nRep;

        if isempty(cp)
            MedianAge_kyr(end+1,1) = NaN;
            Q25_kyr(end+1,1) = NaN;
            Q75_kyr(end+1,1) = NaN;
            P05_kyr(end+1,1) = NaN;
            P95_kyr(end+1,1) = NaN;
            MedianError_kyr(end+1,1) = NaN;
        else
            MedianAge_kyr(end+1,1) = median(cp,'omitnan');
            Q25_kyr(end+1,1) = prctile(cp,25);
            Q75_kyr(end+1,1) = prctile(cp,75);
            P05_kyr(end+1,1) = prctile(cp,5);
            P95_kyr(end+1,1) = prctile(cp,95);
            MedianError_kyr(end+1,1) = median(cp - change_age,'omitnan');
        end
    end
end

Summary = table(OrbitalRatio, PhaseShift_kyr, Phase_100k_deg, OrbitalCase, Method, ...
    DetectionRate_pct, MedianAge_kyr, Q25_kyr, Q75_kyr, P05_kyr, P95_kyr, MedianError_kyr);

disp(' ');
disp(Summary);


% 12. BEAST MODEL-SELECTION SUMMARY


fprintf('\n');
fprintf('BEAST MODEL-SELECTION SUMMARY\n');
fprintf('============================================================\n');

for q = 1:nCase
    valid = isfinite(beast_ncp_mode(:,q));

    if ~any(valid)
        fprintf('%s -> no valid BEAST runs\n', case_names{q});
        continue;
    end

    pct0 = 100 * mean(beast_ncp_mode(valid,q) == 0);
    pct1 = 100 * mean(beast_ncp_mode(valid,q) >= 1);

    fprintf('%s -> 0 CP = %.1f%%; 1 CP = %.1f%%\n', case_names{q}, pct0, pct1);
end


% 13. FIGURE SETTINGS


font_name = 'Cambria Math';

col_sizer = [0.05 0.15 0.95];
col_find  = [1.00 0.45 0.00];
col_beast = [0.72 0.00 0.88];

all_error = [cp_sizer(:) - change_age; cp_find(:) - change_age; cp_beast(:) - change_age];
all_error = all_error(isfinite(all_error));

if isempty(all_error)
    error_limit = 100;
else
    error_limit = max(50, 25 * ceil(prctile(abs(all_error),95) / 25));
end

error_limit = min(error_limit, 400);


% 14. FIGURE


fig = figure( ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [50 40 1250 700], ...
    'Renderer', 'painters');

left_margin = 0.075;
right_margin = 0.025;
top_margin = 0.075;
bottom_margin = 0.085;

row_gap = 0.010;
column_gap = 0.014;
column_widths = [0.40 0.18 0.18 0.18];

usable_width = 1 - left_margin - right_margin - 3 * column_gap;
column_widths = column_widths / sum(column_widths) * usable_width;

x0 = nan(1,4);
x0(1) = left_margin;
for c = 2:4
    x0(c) = x0(c-1) + column_widths(c-1) + column_gap;
end

panel_height = (1 - top_margin - bottom_margin - (nCase - 1) * row_gap) / nCase;

all_example = [example_series(:); clean_trend(:)];
all_example = all_example(isfinite(all_example));
y_min = min(all_example);
y_max = max(all_example);
y_pad = 0.06 * (y_max - y_min);
if y_pad == 0
    y_pad = 1;
end
example_ylim = [y_min - y_pad, y_max + y_pad];

for q = 1:nCase
    y0 = 1 - top_margin - q * panel_height - (q - 1) * row_gap;

    %% Synthetic example
    ax1 = axes('Parent', fig, 'Position', [x0(1), y0, column_widths(1), panel_height]);

    plot(ax1, age, example_series(:,q), 'Color', [0.70 0.70 0.70], 'LineWidth', 0.55);
    hold(ax1, 'on');
    plot(ax1, age, clean_trend, 'k-', 'LineWidth', 1.8);
    xline(ax1, change_age, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.0);

    set(ax1, ...
        'XDir', 'reverse', ...
        'XLim', [0 age_max], ...
        'YLim', example_ylim, ...
        'FontName', font_name, ...
        'FontSize', 18, ...
        'LineWidth', 0.7, ...
        'TickDir', 'out', ...
        'Box', 'on', ...
        'Color', [0.95 0.92 0.94], ...
        'XGrid', 'on', ...
        'YGrid', 'on', ...
        'GridColor', [1 1 1], ...
        'GridAlpha', 1.0, ...
        'GridLineStyle', '-', ...
        'Layer', 'bottom');

    y_range = example_ylim(2) - example_ylim(1);
    ax1.YTick = [example_ylim(1) + 0.25 * y_range, example_ylim(1) + 0.75 * y_range];
    ax1.YTickLabel = compose('%.1f', ax1.YTick);
    ax1.XTick = [0 200 400 600 800];

    if q == 1
        title(ax1, 'Synthetic mMCV trend', 'FontWeight', 'normal', 'FontSize', 18, 'FontName', font_name);
    end

    text(ax1, 0.018, 0.82, panel_labels{q}, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontName', font_name, ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'Clipping', 'on');

    text(ax1, -0.24, 0.82, row_letters{q}, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontName', font_name, ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'Clipping', 'off');

    if q < nCase
        ax1.XTickLabel = [];
    else
        xlabel(ax1, 'Age BP (kyr)', 'FontName', font_name, 'FontSize', 18);
    end

    %% SiZer
    ax2 = axes('Parent', fig, 'Position', [x0(2), y0, column_widths(2), panel_height]);
    draw_distribution(ax2, cp_sizer(:,q), change_age, nRep, col_sizer, error_limit, font_name);
    if q == 1
        title(ax2, 'SiZer', 'FontWeight', 'normal', 'FontSize', 18, 'FontName', font_name);
    end
    if q < nCase
        ax2.XTickLabel = [];
    else
        xlabel(ax2, 'Detected age - 600 kyr (kyr)', 'FontName', font_name, 'FontSize', 18);
    end

    %% findchangepts
    ax3 = axes('Parent', fig, 'Position', [x0(3), y0, column_widths(3), panel_height]);
    draw_distribution(ax3, cp_find(:,q), change_age, nRep, col_find, error_limit, font_name);
    if q == 1
        title(ax3, 'findchangepts', 'FontWeight', 'normal', 'FontSize', 18, 'FontName', font_name);
    end
    if q < nCase
        ax3.XTickLabel = [];
    else
        xlabel(ax3, 'Detected age - 600 kyr (kyr)', 'FontName', font_name, 'FontSize', 18);
    end

    %% BEAST
    ax4 = axes('Parent', fig, 'Position', [x0(4), y0, column_widths(4), panel_height]);
    draw_distribution(ax4, cp_beast(:,q), change_age, nRep, col_beast, error_limit, font_name);
    if q == 1
        title(ax4, 'BEAST', 'FontWeight', 'normal', 'FontSize', 18, 'FontName', font_name);
    end
    if q < nCase
        ax4.XTickLabel = [];
    else
        xlabel(ax4, 'Detected age - 600 kyr (kyr)', 'FontName', font_name, 'FontSize', 18);
    end
end

annotation(fig, 'textbox', [0.035 0.42 0.03 0.20], ...
    'String', 'mMCV (z-score)', ...
    'LineStyle', 'none', ...
    'Rotation', 90, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontName', font_name, ...
    'FontSize', 18, ...
    'FontWeight', 'bold');

annotation(fig, 'textbox', [0.060 0.958 0.40 0.03], ...
    'String', 'Changepoint tests: orbital-forcing sensitivity', ...
    'LineStyle', 'none', ...
    'FontName', font_name, ...
    'FontSize', 18, ...
    'FontWeight', 'normal');

fprintf('\nFigure created. No files were saved.\n');


% LOCAL FUNCTION: RUN BEAST ONCE


function [cp_age,ncp_mode] = run_beast_once( ...
    y,age,dt,age_max,tcp_minmax,torder_minmax,seed, ...
    mcmc_samples,mcmc_burnin,mcmc_chains,mcmc_thin)

cp_age = NaN;
ncp_mode = NaN;

out = beast( ...
    y(:), ...
    'start',0, ...
    'deltat',dt, ...
    'season','none', ...
    'tcp.minmax',tcp_minmax, ...
    'torder.minmax',torder_minmax, ...
    'mcmc.seed',seed, ...
    'mcmc.samples',mcmc_samples, ...
    'mcmc.burnin',mcmc_burnin, ...
    'mcmc.chains',mcmc_chains, ...
    'mcmc.thin',mcmc_thin, ...
    'print.progress',false, ...
    'print.param',false, ...
    'print.warning',false, ...
    'quiet',true, ...
    'dump.ci',false);

if ~isfield(out,'trend')
    return
end

if isfield(out.trend,'ncp_mode') && ~isempty(out.trend.ncp_mode)
    z = double(out.trend.ncp_mode(:));
    z = z(isfinite(z));
    if ~isempty(z)
        ncp_mode = z(1);
    end
end

if ~isfinite(ncp_mode) && isfield(out.trend,'ncpPr') && ~isempty(out.trend.ncpPr)
    p = double(out.trend.ncpPr(:));
    [~,ii] = max(p);
    ncp_mode = ii - 1;
end

if ~isfinite(ncp_mode) || ncp_mode < 1
    return
end

if isfield(out.trend,'cp') && ~isempty(out.trend.cp)
    cp = double(out.trend.cp(:));

    if isfield(out.trend,'cpPr') && ~isempty(out.trend.cpPr)
        cpPr = double(out.trend.cpPr(:));
    else
        cpPr = nan(size(cp));
    end

    valid = isfinite(cp) & cp >= 0 & cp <= age_max;

    if any(valid)
        cp = cp(valid);

        if numel(cpPr) == numel(valid)
            cpPr = cpPr(valid);
        else
            cpPr = nan(size(cp));
        end

        if any(isfinite(cpPr))
            [~,ii] = max(cpPr);
            cp_age = cp(ii);
        else
            cp_age = cp(1);
        end
    end
end

if ~isfinite(cp_age) && isfield(out.trend,'cpOccPr') && ~isempty(out.trend.cpOccPr)
    p = double(out.trend.cpOccPr(:));
    if numel(p) == numel(age)
        [~,ii] = max(p);
        cp_age = age(ii);
    end
end

end


% LOCAL FUNCTION: DRAW DETECTION DISTRIBUTION


function draw_distribution(ax,cp_age,true_age,nRep,colour,error_limit,font_name)

hold(ax,'on');

xline(ax,0,'--', 'Color',[0.35 0.35 0.35], 'LineWidth',0.9);

valid = cp_age(isfinite(cp_age));
err = valid - true_age;

if isempty(err)
    text(ax,0.50,0.58,'ND', ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName',font_name, ...
        'FontSize',18, ...
        'FontWeight','bold', ...
        'Color',colour);
else
    p05 = prctile(err,5);
    p95 = prctile(err,95);
    q25 = prctile(err,25);
    q75 = prctile(err,75);
    med = median(err,'omitnan');

    plot(ax,[p05 p95],[0.56 0.56],'-', 'Color',colour,'LineWidth',1.0);
    plot(ax,[q25 q75],[0.56 0.56],'-', 'Color',colour,'LineWidth',3.0);
    plot(ax,med,0.56,'s', 'MarkerSize',8, 'MarkerFaceColor','w', ...
        'MarkerEdgeColor',colour, 'LineWidth',1.5);
end

% Detection percentage in the upper-right corner
text(ax,0.97,0.87, sprintf('%.0f%%',100*numel(valid)/nRep), ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','middle', ...
    'FontName',font_name, ...
    'FontSize',18, ...
    'FontWeight','bold', ...
    'Color',[0.25 0.25 0.25]);

set(ax, ...
    'XLim',[-error_limit error_limit], ...
    'YLim',[0 1], ...
    'YTick',[0.25 0.75], ...
    'YTickLabel',[], ...
    'FontName',font_name, ...
    'FontSize',18, ...
    'LineWidth',0.7, ...
    'TickDir','out', ...
    'TickLength',[0 0], ...
    'XMinorTick','off', ...
    'YMinorTick','off', ...
    'Box','on', ...
    'Color',[0.90 0.94 0.99], ...
    'XGrid','on', ...
    'YGrid','off', ...
    'GridColor',[0.72 0.72 0.72], ...
    'GridAlpha',0.30, ...
    'GridLineStyle','-', ...
    'Layer','bottom');

ax.XTick = [-0.50 * error_limit, 0, 0.50 * error_limit];
ax.XTickLabel = compose('%.0f', ax.XTick);

end
