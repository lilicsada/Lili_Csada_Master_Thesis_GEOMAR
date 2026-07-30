clear; close all; clc;


%  1. SETTINGS
opts.age_min = 0;
opts.age_max = 800; %Only investigate the past 800kyr

opts.interp_method = "pchip"; %PCHIP interpolation is used

% Band-pass filter
opts.period_min = 2;          % kyr
opts.period_max = 10;         % kyr
opts.filter_order = 4;

% Long-term trend
opts.gauss_window_kyr = 100;  % kyr

% Padding and edge cut
opts.pad_kyr = 10;           % kyr mirror padding before filtering and hilbert amplitude
opts.edge_cut_kyr = 0;       % Extra cutting in the edges if needed based on methodological biases

% Plot high-frequency signal in proxy units
opts.zscore_highfreq = false;

% Plot Hilbert amplitude as z-score
opts.zscore_hilbert = true;

% Mark data gaps larger than this as hiatus 
opts.hiatus_threshold_kyr = 5; % kyr

% EPICA Dome C cut
% The first record is EPICA Dome C. Keep only the first 5378 valid
% datapoints. Due to resolution issues
opts.epica_cut_after_point = 5378;

% LR04 MIS boundary file
% To the glacial-interglacial boundaries
opts.mis_file = "LR04_MISboundaries.txt";



%  2. LOAD DATA

files = [
    "Epica_deuterium.txt"
    "Fuji_deuterium.txt"
    "Speleothem_D18O.txt"
    "1094_N. pachydermaD18O.txt"
    "1089_G.bulloidesD18O.txt"
    "983_Nps%.txt"
    "PS75_076_lnZrRb.txt"
    "U1385_logCaTi.txt"
    "U1308_logSiSr.txt"
    
];

labels = [
    "EPICA"
    "Dome Fuji"
    "Speleothem"
    "ODP 1094"
    "ODP 1089"
    "ODP 983"
    "PS75/076"
    "U1385"
    "U1308"
];

% Interpolation number, either 500 or 250 yr
dt_list = [
    0.5
    0.5
    0.5
    0.5
    0.5
    0.5
    0.5
    0.25
    0.25
];


%3. DEFINING GLACIALS-INTERGLACIALS
mis_intervals = read_lr04_mis_intervals(opts.mis_file, opts.age_min, opts.age_max);

%4. MAIN ANALYSIS LINE
R = run_pipeline(files, labels, dt_list, opts);

%5. PLOTTING SECTION
plot_opts = make_plot_options();
plot_stacked_figure(R, mis_intervals, opts, plot_opts, "raw");
plot_stacked_figure(R, mis_intervals, opts, plot_opts, "interp");
plot_stacked_figure(R, mis_intervals, opts, plot_opts, "highfreq");
plot_stacked_figure(R, mis_intervals, opts, plot_opts, "hilbert");


% SAVE OUTPUT FOR SEPARATE SPECTRAL ANALYSIS
save( ...
    "Main_Data_results.mat", ...
    "R", ...
    "opts", ...
    "plot_opts", ...
    "mis_intervals", ...
    "-v7.3");

fprintf("\nPipeline output saved as Main_Datat_results.mat\n\n");

% 6. FUNCTIONS: MAIN PIPELINE
function R = run_pipeline(files, labels, dt_list, opts)
    nrec = numel(files);
    R = struct();
    for i = 1:nrec
        fname = find_existing_file(files(i));

        % Apply EPICA cut only to the first record, before analysis
        % Optional step: I left it in the code in case other datasets also need to be cut.
        if i == 1
            cut_after_point = opts.epica_cut_after_point;
        else
            cut_after_point = NaN;
        end

        [age_raw, proxy_raw, cut_info] = read_proxy_file( ...
            fname, opts.age_min, opts.age_max, cut_after_point);

        if cut_info.applied
            fprintf("\nEPICA Dome C cut applied before analysis:\n");
            fprintf("  Kept first %d valid datapoints.\n", cut_info.point);
            fprintf("  Last kept datapoint age = %.3f kyr BP = %.0f years BP.\n\n", ...
                cut_info.age_kyr, cut_info.age_kyr * 1000);
        end

        
        % Detect hiatus 
        raw_gaps = diff(age_raw);
        hiatus_idx = find(raw_gaps > opts.hiatus_threshold_kyr);
        hiatus_start = age_raw(hiatus_idx);
        hiatus_end   = age_raw(hiatus_idx + 1);
        dt = dt_list(i);
        age_reg = (ceil(min(age_raw)/dt)*dt : dt : floor(max(age_raw)/dt)*dt)';
        proxy_interp = interp1(age_raw, proxy_raw, age_reg, opts.interp_method);
        proxy_trend = gaussian_smooth_nan(proxy_interp, dt, opts.gauss_window_kyr);

        
        % Full band-pass filtered signal and Hilbert amplitude
        % Padding is added before filtering, kept during Hilbert amplitude
        % calculation, and removed only after Hilbert amplitude calculation.
     
        [proxy_bp_full, hilbert_amp_full] = bandpass_and_hilbert_with_padding( ...
            proxy_interp, dt, opts.period_min, opts.period_max, ...
            opts.filter_order, opts.pad_kyr);
     
        % Full high-frequency signal for plotting
        % The artificial padding has already been removed after Hilbert amplitude calculation.
        
        if opts.zscore_highfreq
            proxy_bp_plot_full = zscore_nan(proxy_bp_full);
        else
            proxy_bp_plot_full = proxy_bp_full;
        end

       
        % Full Hilbert amplitude envelope
        % Calculated from the padded 2–10 kyr band-pass filtered signal.
        % The artificial padding has already been removed here.
       
        % Z-score Hilbert amplitude
        if opts.zscore_hilbert
            hilbert_amp_plot_full = zscore_nan(hilbert_amp_full);
        else
            hilbert_amp_plot_full = hilbert_amp_full;
        end

        % Smooth the plotted Hilbert amplitude version
        % Here this means the 100 kyr smooth is applied to the z-scored Hilbert amplitude.
        hilbert_amp_smooth_full = gaussian_smooth_nan( ...
            hilbert_amp_plot_full, dt, opts.gauss_window_kyr);

        
        % Define edge cut only for plotting
        % Since opts.edge_cut_kyr = 0, no additional real-data edge is removed.
         plot_edge_keep = age_reg >= min(age_reg) + opts.edge_cut_kyr & ...
                         age_reg <= max(age_reg) - opts.edge_cut_kyr;

        
        % Save everything into record structure
       
        R(i).label = labels(i);
        R(i).file = fname;
        R(i).dt = dt;

        R(i).age_raw = age_raw;
        R(i).proxy_raw = proxy_raw;

        R(i).age_reg = age_reg;
        R(i).proxy_interp = proxy_interp;
        R(i).proxy_trend = proxy_trend;

        R(i).age_bp = age_reg;
        R(i).proxy_bp = proxy_bp_plot_full;

        R(i).age_hilbert = age_reg;
        R(i).hilbert_amp_z = hilbert_amp_plot_full;
        R(i).hilbert_amp_z_smooth = hilbert_amp_smooth_full;

        R(i).plot_edge_keep = plot_edge_keep;

        R(i).hiatus_start = hiatus_start;
        R(i).hiatus_end = hiatus_end;

        fprintf("%s | raw points: %d | dt: %.2f kyr | interpolated points: %d | hiatus gaps > %.1f kyr: %d\n", ...
            char(labels(i)), numel(age_raw), dt, numel(age_reg), ...
            opts.hiatus_threshold_kyr, numel(hiatus_start));

    end

end


function [age_raw, proxy_raw, cut_info] = read_proxy_file(filename, age_min, age_max, cut_after_point)
    data = readmatrix(filename);
    if size(data,2) < 2
        error("File %s must contain at least two columns: age and proxy.", filename);
    end

    age_raw = data(:,1);
    proxy_raw = data(:,2);

    % Convert years BP to kyr BP if necessary
    if max(age_raw, [], "omitnan") > 5000
        age_raw = age_raw / 1000;
    end

    valid = ~isnan(age_raw) & ~isnan(proxy_raw);
    age_raw = age_raw(valid);
    proxy_raw = proxy_raw(valid);

    cut_info.applied = false;
    cut_info.point = NaN;
    cut_info.age_kyr = NaN;

   
    % Optional cut before sorting/interpolation/filtering
    % This keeps the first N valid datapoints in the original file order.
    
    if nargin >= 4 && ~isnan(cut_after_point)

        if numel(age_raw) < cut_after_point
            error("Cannot cut %s after datapoint %d because the file only has %d valid datapoints.", ...
                filename, cut_after_point, numel(age_raw));
        end

        cut_info.applied = true;
        cut_info.point = cut_after_point;
        cut_info.age_kyr = age_raw(cut_after_point);

        age_raw = age_raw(1:cut_after_point);
        proxy_raw = proxy_raw(1:cut_after_point);

    end

    [age_raw, idx] = sort(age_raw);
    proxy_raw = proxy_raw(idx);

    [age_raw, unique_idx] = unique(age_raw, "stable");
    proxy_raw = proxy_raw(unique_idx);

    keep = age_raw >= age_min & age_raw <= age_max;
    age_raw = age_raw(keep);
    proxy_raw = proxy_raw(keep);

    if isempty(age_raw)
        error("No data left in %s after cutting to %.0f–%.0f kyr.", ...
            filename, age_min, age_max);
    end

end


function [y_bp, hilbert_amp] = bandpass_and_hilbert_with_padding( ...
    y, dt, period_min, period_max, filter_order, pad_kyr)
    y = y(:);
    fs = 1 / dt;
    f_low  = 1 / period_max;
    f_high = 1 / period_min;
    wn = [f_low f_high] / (fs/2);
    [b, a] = butter(filter_order, wn, "bandpass");
    npad = round(pad_kyr / dt);

    % Minimum padding required for filtfilt stability
    minpad = 3 * max(length(a), length(b));
    npad = max(npad, minpad);
    if npad >= length(y)
        error("Padding is longer than the data series. Reduce pad_kyr.");
    end

    % Add mirror padding before filtering
    y_pad = mirror_pad(y, npad);

    % Apply band-pass filter to the padded signal
    y_bp_pad = filtfilt(b, a, y_pad);

    % Calculate Hilbert amplitude while the padding is still present
    hilbert_amp_pad = abs(hilbert(y_bp_pad));

    % Remove artificial mirror padding only after Hilbert amplitude calculation
    y_bp = y_bp_pad(npad+1 : npad+length(y));
    hilbert_amp = hilbert_amp_pad(npad+1 : npad+length(y));

end


function y_pad = mirror_pad(y, npad)
    y = y(:);
    n = length(y);
    left_pad = y(npad+1:-1:2);
    right_pad = y(n-1:-1:n-npad);
    y_pad = [left_pad; y; right_pad];

end


function y_smooth = gaussian_smooth_nan(y, dt, window_kyr)
    y = y(:);
    window_points = round(window_kyr / dt);
    if mod(window_points, 2) == 0
        window_points = window_points + 1;
    end

    half_window = floor(window_points / 2);
    x = (-half_window:half_window)';
    sigma_points = window_points / 6;
    kernel = exp(-0.5 * (x / sigma_points).^2);
    kernel = kernel / sum(kernel);
    valid = ~isnan(y);
    y0 = y;
    y0(~valid) = 0;
    y_conv = conv(y0, kernel, 'same');
    w_conv = conv(double(valid), kernel, 'same');
    y_smooth = y_conv ./ w_conv;
    y_smooth(w_conv == 0) = NaN;

end


function z = zscore_nan(x)

    mu = mean(x, "omitnan");
    sig = std(x, 0, "omitnan");

    if sig == 0 || isnan(sig)
        z = x - mu;
    else
        z = (x - mu) ./ sig;
    end

end

%7. PLOTTING

function plot_opts = make_plot_options()

    plot_opts.fig_position = [60 20 600 1350];

    % Font
    plot_opts.font_name = "Cambria Math";

    % Y-axis labels for plots shown in proxy units
    % Used for the interpolated and high-frequency figures.
    plot_opts.proxy_ylabels = [
        "EPICA \delta D (‰)"
        "Dome Fuji \delta D (‰)"
        "Speleothem \delta^{18}O (‰)"
        "ODP 1094 \delta^{18}O (‰)"
        "ODP 1089 \delta^{18}O (‰)"
        "ODP 983 N. pachyderma (s.) (%)"
        "PS75/076 ln(Zr/Rb)"
        "U1385 log(Ca/Ti)"
        "U1308 log(Si/Sr)"
    ];

    % The Hilbert-amplitude series are z-standardised and therefore unitless.
    plot_opts.hilbert_ylabels = [
        "EPICA (z-score)"
        "Dome Fuji (z-score)"
        "Speleothem (z-score)"
        "ODP 1094 (z-score)"
        "ODP 1089 (z-score)"
        "ODP 983 (z-score)"
        "PS75/076 (z-score)"
        "U1385 (z-score)"
        "U1308 (z-score)"
    ];

    % Axis background is white, so the MIS shading provides the background colour
    plot_opts.ax_bg = [1 1 1];

    % Custom grid is drawn above MIS shading but below hiatus and data
    plot_opts.grid_col = [1 1 1];
    plot_opts.grid_line_width = 0.65;

    % Signals are black
    plot_opts.raw_col = [0 0 0];
    plot_opts.interp_col = [0 0 0];
    plot_opts.hf_col = [0 0 0];
    plot_opts.hilbert_col = [0 0 0];

    % Long-term trend / smoothed amplitude: softer and thinner red
    plot_opts.trend_col = [0.95 0.60 0.60];
    plot_opts.trend_line_width = 1.0;

    % MIS shading colours
    plot_opts.interglacial_col = [0.975 0.935 0.945];  % pink
    plot_opts.glacial_col      = [0.88 0.93 1.00];     % baby blue

    % MIS label colours
    plot_opts.mis_interglacial_text_col = [0.55 0.10 0.22];  % darker pink/red
    plot_opts.mis_glacial_text_col      = [0.10 0.25 0.55];  % darker blue
    plot_opts.mis_label_font_size = 6;

    % Hiatus / large raw-data gap shading
    plot_opts.hiatus_col = [0.35 0.35 0.35];
    plot_opts.hiatus_alpha = 0.22;

    % Font sizes
    plot_opts.axis_font_size = 9;
    plot_opts.xlabel_font_size = 10;
    plot_opts.ylabel_font_size = 9;
    plot_opts.panel_label_font_size = 11;
    plot_opts.inside_text_font_size = 9;
    plot_opts.title_font_size = 11;
    plot_opts.legend_font_size = 8;

    % Line and marker sizes
    plot_opts.raw_marker_size = 3;
    plot_opts.interp_line_width = 0.65;
    plot_opts.hf_line_width = 0.65;
    plot_opts.hilbert_line_width = 0.65;

    % Same number of y-axis divisions in each panel
    plot_opts.n_yticks = 5;

    % Vertical grid spacing in kyr
    plot_opts.xgrid_spacing = 100;

end


function fig = plot_stacked_figure(R, mis_intervals, opts, plot_opts, mode)
    nrec = numel(R);
    fig = figure('Color', 'w', 'Position', plot_opts.fig_position);

    % Manual stacked-axis layout
    left_margin   = 0.075;
    right_margin  = 0.020;
    top_margin    = 0.080;
    bottom_margin = 0.055;

    gap = 0.002;
    ax_width  = 1 - left_margin - right_margin;
    ax_height = (1 - top_margin - bottom_margin - (nrec-1)*gap) / nrec;
    ax = gobjects(nrec,1);
    for i = 1:nrec
        y0 = 1 - top_margin - i*ax_height - (i-1)*gap;
        ax(i) = axes('Parent', fig, ...
            'Position', [left_margin y0 ax_width ax_height]);

        hold(ax(i), 'on');
        switch mode
            case "raw"
                xdata = R(i).age_raw;
                ydata = R(i).proxy_raw;
                xtrend = [];
                ytrend = [];
                inside_txt = sprintf('%s raw data', R(i).label);
                inside_col = [0.55 0.00 0.10];

            case "interp"
                xdata = R(i).age_reg;
                ydata = R(i).proxy_interp;
                xtrend = R(i).age_reg;
                ytrend = R(i).proxy_trend;
                inside_txt = sprintf('%s + 100 kyr trend', R(i).label);
                inside_col = [0.70 0.25 0.25];

            case "highfreq"

                idx_plot = R(i).plot_edge_keep;
                xdata = R(i).age_bp(idx_plot);
                ydata = R(i).proxy_bp(idx_plot);
                xtrend = [];
                ytrend = [];
                inside_txt = sprintf('%s: 2–10 kyr signal', R(i).label);
                inside_col = [0.05 0.20 0.60];

            case "hilbert"

                idx_plot = R(i).plot_edge_keep;
                xdata = R(i).age_hilbert(idx_plot);
                ydata = R(i).hilbert_amp_z(idx_plot);
                xtrend = R(i).age_hilbert(idx_plot);
                ytrend = R(i).hilbert_amp_z_smooth(idx_plot);
                inside_txt = sprintf('%s: Hilbert amplitude', R(i).label);
                inside_col = [0.70 0.25 0.25];

        end

        xlim(ax(i), [opts.age_min opts.age_max]);

        % Set y limits before drawing shading/grid/data
        if (mode == "interp" || mode == "hilbert") && ~isempty(ytrend)
            all_y = [ydata(:); ytrend(:)];
        else
            all_y = ydata(:);
        end

        set_nice_ylim(ax(i), all_y, mode);
        style_axis(ax(i), plot_opts);

        % Same number of y-tick divisions, with first/last labels hidden
        format_y_axis_same_divisions(ax(i), plot_opts.n_yticks, mode);

        add_mis_shading(ax(i), mis_intervals, plot_opts);
        add_custom_grid(ax(i), opts, plot_opts);
        add_hiatus_markers(ax(i), R(i), plot_opts);

        switch mode
            case "raw"
                plot(ax(i), xdata, ydata, '.', ...
                    'Color', plot_opts.raw_col, ...
                    'MarkerSize', plot_opts.raw_marker_size);
            case "interp"
                h1 = plot(ax(i), xdata, ydata, ...
                    'Color', plot_opts.interp_col, ...
                    'LineWidth', plot_opts.interp_line_width);
                h2 = plot(ax(i), xtrend, ytrend, ...
                    'Color', plot_opts.trend_col, ...
                    'LineWidth', plot_opts.trend_line_width);
                if i == 1
                    h_hiatus = make_hiatus_legend_handle(ax(i), plot_opts);
                    legend(ax(i), [h1 h2 h_hiatus], ...
                        {'Interpolated record', '100 kyr trend', ...
                         sprintf('Hiatus (>%.0f kyr gap)', opts.hiatus_threshold_kyr)}, ...
                        'Location', 'northeast', ...
                        'Box', 'off', ...
                        'FontSize', plot_opts.legend_font_size, ...
                        'FontName', plot_opts.font_name);
                end

            case "highfreq"

                h1 = plot(ax(i), xdata, ydata, ...
                    'Color', plot_opts.hf_col, ...
                    'LineWidth', plot_opts.hf_line_width);

                if i == 1
                    h_hiatus = make_hiatus_legend_handle(ax(i), plot_opts);

                    legend(ax(i), [h1 h_hiatus], ...
                        {'2–10 kyr signal', ...
                         sprintf('Hiatus (>%.0f kyr gap)', opts.hiatus_threshold_kyr)}, ...
                        'Location', 'northeast', ...
                        'Box', 'off', ...
                        'FontSize', plot_opts.legend_font_size, ...
                        'FontName', plot_opts.font_name);
                end

            case "hilbert"

                h1 = plot(ax(i), xdata, ydata, ...
                    'Color', plot_opts.hilbert_col, ...
                    'LineWidth', plot_opts.hilbert_line_width);

                h2 = plot(ax(i), xtrend, ytrend, ...
                    'Color', plot_opts.trend_col, ...
                    'LineWidth', plot_opts.trend_line_width);

                if i == 1
                    h_hiatus = make_hiatus_legend_handle(ax(i), plot_opts);

                    legend(ax(i), [h1 h2 h_hiatus], ...
                        {'Hilbert amplitude (z-score)', '100 kyr smooth', ...
                         sprintf('Hiatus (>%.0f kyr gap)', opts.hiatus_threshold_kyr)}, ...
                        'Location', 'northeast', ...
                        'Box', 'off', ...
                        'FontSize', plot_opts.legend_font_size, ...
                        'FontName', plot_opts.font_name);
                end

        end

        % Add MIS labels only to the top panel to avoid clutter
        if i == 1
            add_mis_labels(ax(i), mis_intervals, plot_opts);
        end

        % Use proxy variables and units only for the processed records
        % shown in physical proxy units. Hilbert amplitudes are unitless
        % because they are z-standardised. Raw data retain short record names.
        switch mode
            case {"interp", "highfreq"}
                y_label_text = plot_opts.proxy_ylabels(i);
                y_label_interpreter = 'tex';

            case "hilbert"
                y_label_text = plot_opts.hilbert_ylabels(i);
                y_label_interpreter = 'tex';

            otherwise
                y_label_text = R(i).label;
                y_label_interpreter = 'none';
        end

        ylabel(ax(i), y_label_text, ...
            'Interpreter', y_label_interpreter, ...
            'FontSize', plot_opts.ylabel_font_size, ...
            'FontWeight', 'bold', ...
            'FontName', plot_opts.font_name);

        add_panel_label(ax(i), i, plot_opts);

        % EPICA label is shifted lower to avoid overlap with the top axis
        add_inside_text(ax(i), inside_txt, inside_col, plot_opts, i);

        % Only show time axis at top and bottom
        if i == 1

            ax(i).XAxisLocation = 'top';
            ax(i).XAxis.Visible = 'on';
            xlabel(ax(i), 'Age BP(kyr)', ...
                'FontSize', plot_opts.xlabel_font_size, ...
                'FontWeight', 'bold', ...
                'FontName', plot_opts.font_name);

        elseif i == nrec

            ax(i).XAxisLocation = 'bottom';
            ax(i).XAxis.Visible = 'on';
            xlabel(ax(i), 'Age BP(kyr)', ...
                'FontSize', plot_opts.xlabel_font_size, ...
                'FontWeight', 'bold', ...
                'FontName', plot_opts.font_name);

        else

            ax(i).XAxis.Visible = 'off';

        end

    end

    linkaxes(ax, 'x');

    switch mode

        case "raw"
            main_title = "Raw proxy data points";

        case "interp"
            main_title = "Interpolated proxy records";

        case "highfreq"

            if opts.zscore_highfreq
                unit_txt = "z-scored";
            else
                unit_txt = "proxy units";
            end

            main_title = sprintf("High-frequency signal: 2–10 kyr band-pass filtered records (%s)", unit_txt);

        case "hilbert"

            if opts.zscore_hilbert
                unit_txt = "z-scored";
            else
                unit_txt = "proxy units";
            end

            main_title = sprintf("Hilbert amplitude envelope of 2–10 kyr variability (%s)", unit_txt);

    end

    sgtitle(fig, main_title, ...
        'FontWeight', 'bold', ...
        'FontSize', plot_opts.title_font_size, ...
        'FontName', plot_opts.font_name);

    % Force all remaining text/axis objects to Cambria Math
    set(findall(fig, '-property', 'FontName'), 'FontName', plot_opts.font_name);

end


function set_nice_ylim(ax, y, mode)
    y = y(~isnan(y));
    if isempty(y)
        ylim(ax, [-1 1]);
        return
    end

    ymin = min(y);
    ymax = max(y);
    yr = ymax - ymin;

    if yr == 0
        yr = max(abs(ymax), 1);
        ymin = ymin - 0.5 * yr;
        ymax = ymax + 0.5 * yr;
    end

    pad = 0.08 * yr;

    ylim(ax, [ymin - pad, ymax + pad]);

end


function format_y_axis_same_divisions(ax, nticks, mode)
    yl = ylim(ax);
    yt = linspace(yl(1), yl(2), nticks);
    ax.YTick = yt;
    tick_labels = strings(size(yt));

    for k = 1:numel(yt)

        if abs(yt(k) - round(yt(k))) < 1e-8
            tick_labels(k) = sprintf('%.0f', yt(k));
        else
            tick_labels(k) = sprintf('%.1f', yt(k));
        end

    end

    % Hide first and last y-axis tick labels
    tick_labels(1) = "";
    tick_labels(end) = "";

    ax.YTickLabel = tick_labels;

end


function add_mis_shading(ax, mis_intervals, plot_opts)

    yl = ylim(ax);

    for j = 1:height(mis_intervals)

        x1 = mis_intervals.start_age(j);
        x2 = mis_intervals.end_age(j);

        if mis_intervals.is_interglacial(j)
            col = plot_opts.interglacial_col;
        else
            col = plot_opts.glacial_col;
        end

        patch(ax, ...
            [x1 x2 x2 x1], ...
            [yl(1) yl(1) yl(2) yl(2)], ...
            col, ...
            'FaceAlpha', 1.0, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

    end

    ylim(ax, yl);

end


function add_mis_labels(ax, mis_intervals, plot_opts)

    yl = ylim(ax);

    % Labels near the top of the panel
    y_pos = yl(2) - 0.12 * range(yl);

    for j = 1:height(mis_intervals)

        x1 = mis_intervals.start_age(j);
        x2 = mis_intervals.end_age(j);
        x_mid = (x1 + x2) / 2;

        mis_number = mis_intervals.MIS(j);

        if mis_intervals.is_interglacial(j)
            txt_col = plot_opts.mis_interglacial_text_col;
        else
            txt_col = plot_opts.mis_glacial_text_col;
        end

        text(ax, x_mid, y_pos, sprintf('%g', mis_number), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontSize', plot_opts.mis_label_font_size, ...
            'FontWeight', 'bold', ...
            'FontName', plot_opts.font_name, ...
            'Color', txt_col, ...
            'Clipping', 'on');

    end

    ylim(ax, yl);

end

function add_hiatus_markers(ax, record, plot_opts)

    yl = ylim(ax);

    if ~isfield(record, 'hiatus_start') || isempty(record.hiatus_start)
        return
    end

    for h = 1:numel(record.hiatus_start)

        x1 = record.hiatus_start(h);
        x2 = record.hiatus_end(h);

        patch(ax, ...
            [x1 x2 x2 x1], ...
            [yl(1) yl(1) yl(2) yl(2)], ...
            plot_opts.hiatus_col, ...
            'FaceAlpha', plot_opts.hiatus_alpha, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

    end

    ylim(ax, yl);

end


function h = make_hiatus_legend_handle(ax, plot_opts)

    % Dummy patch used only to explain the grey hiatus shading in the legend.
    h = patch(ax, ...
        [NaN NaN NaN NaN], ...
        [NaN NaN NaN NaN], ...
        plot_opts.hiatus_col, ...
        'FaceAlpha', plot_opts.hiatus_alpha, ...
        'EdgeColor', 'none');

end


function add_custom_grid(ax, opts, plot_opts)

    xl = xlim(ax);
    yl = ylim(ax);

    % Vertical grid every 100 kyr
    xgrid = opts.age_min:plot_opts.xgrid_spacing:opts.age_max;

    for k = 1:numel(xgrid)

        plot(ax, [xgrid(k) xgrid(k)], yl, ...
            'Color', plot_opts.grid_col, ...
            'LineWidth', plot_opts.grid_line_width, ...
            'HandleVisibility', 'off');

    end

    % Horizontal grid at y-ticks
    yt = ax.YTick;

    for k = 1:numel(yt)

        plot(ax, xl, [yt(k) yt(k)], ...
            'Color', plot_opts.grid_col, ...
            'LineWidth', plot_opts.grid_line_width, ...
            'HandleVisibility', 'off');

    end

    xlim(ax, xl);
    ylim(ax, yl);

end


function style_axis(ax, plot_opts)

    set(ax, ...
        'Color', plot_opts.ax_bg, ...
        'FontName', plot_opts.font_name, ...
        'FontSize', plot_opts.axis_font_size, ...
        'LineWidth', 0.5, ...
        'TickDir', 'out', ...
        'Layer', 'top', ...
        'XGrid', 'off', ...
        'YGrid', 'off', ...
        'Box', 'off', ...
        'XColor', [0.35 0.35 0.35], ...
        'YColor', [0.35 0.35 0.35]);

end


function add_panel_label(ax, i, plot_opts)

    letters = 'abcdefghijklmnopqrstuvwxyz';
    txt = sprintf('(%c)', letters(i));

    text(ax, -0.12, 0.84, txt, ...
        'Units', 'normalized', ...
        'FontSize', plot_opts.panel_label_font_size, ...
        'FontWeight', 'bold', ...
        'FontName', plot_opts.font_name, ...
        'Color', [0 0 0], ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle');

end


function add_inside_text(ax, txt, col, plot_opts, i)

    xl = xlim(ax);
    yl = ylim(ax);

    x_pos = xl(1) + 0.02 * range(xl);

    % Shift only the first panel label lower
    if i == 1
        y_pos = yl(2) - 0.22 * range(yl);
    else
        y_pos = yl(2) - 0.12 * range(yl);
    end

    text(ax, x_pos, y_pos, txt, ...
        'FontSize', plot_opts.inside_text_font_size, ...
        'FontWeight', 'bold', ...
        'FontName', plot_opts.font_name, ...
        'Color', col, ...
        'Interpreter', 'none', ...
        'VerticalAlignment', 'top');

end




% 7. FUNCTIONS: LR04 MIS BOUNDARIES

function mis_intervals = read_lr04_mis_intervals(filename, age_min, age_max)

    filename = find_existing_file(filename);

    fid = fopen(filename, 'r');

    if fid == -1
        error("Could not open MIS file: %s", filename);
    end

    lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);

    lines = string(lines{1});

    boundary_age = [];
    stage_left = [];
    stage_right = [];

    for i = 1:numel(lines)

        line = strtrim(lines(i));

        tok = regexp(line, '^(\d+)/(\d+)\s+([0-9]+\.?[0-9]*)', ...
            'tokens', 'once');

        if ~isempty(tok)

            left_stage = str2double(tok{1});
            right_stage = str2double(tok{2});
            age = str2double(tok{3});

            if right_stage == left_stage + 1
                stage_left(end+1,1) = left_stage;
                stage_right(end+1,1) = right_stage;
                boundary_age(end+1,1) = age;
            end

        end

    end

    [boundary_age, idx] = sort(boundary_age);
    stage_left = stage_left(idx);
    stage_right = stage_right(idx);

    keep = boundary_age > age_min & boundary_age < age_max;
    boundary_age = boundary_age(keep);
    stage_left = stage_left(keep);
    stage_right = stage_right(keep);

    if isempty(boundary_age)
        mis_intervals = table(age_min, age_max, 1, true, ...
            'VariableNames', {'start_age','end_age','MIS','is_interglacial'});
        return
    end

    starts = [];
    stops = [];
    stages = [];

    starts(end+1,1) = age_min;
    stops(end+1,1) = boundary_age(1);
    stages(end+1,1) = stage_left(1);

    for j = 1:length(boundary_age)

        starts(end+1,1) = boundary_age(j);

        if j < length(boundary_age)
            stops(end+1,1) = boundary_age(j+1);
        else
            stops(end+1,1) = age_max;
        end

        stages(end+1,1) = stage_right(j);

    end

    is_interglacial = mod(stages,2) == 1;

    mis_intervals = table(starts, stops, stages, is_interglacial, ...
        'VariableNames', {'start_age','end_age','MIS','is_interglacial'});

end



%  FUNCTION: FILE HANDLING


function fname_out = find_existing_file(fname_in)

    fname_in = string(fname_in);

    if isfile(fname_in)
        fname_out = fname_in;
        return
    end

    d = dir;
    all_names = string({d.name});

    idx = find(strcmpi(all_names, fname_in), 1);

    if isempty(idx)
        error("File not found: %s", fname_in);
    end

    fname_out = all_names(idx);

end