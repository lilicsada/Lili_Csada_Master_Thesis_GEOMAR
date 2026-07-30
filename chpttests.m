function results = chpttests(age, data, opts)
%CHPTTESTS Run SiZer, findchangepts, and BEAST

    if nargin < 2
        error([ ...
            'chpttests requires an age vector and a data vector.' newline ...
            'Example: results = chpttests(age, data);']);
    end

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    opts = apply_defaults(opts);

    %% Prepare input
    age = age(:);
    data = data(:);

    if numel(age) ~= numel(data)
        error('age and data must contain the same number of values.');
    end

    valid = isfinite(age) & isfinite(data);
    age = age(valid);
    data = data(valid);

    [age, order] = sort(age, 'ascend');
    data = data(order);

    [age, unique_index] = unique(age, 'stable');
    data = data(unique_index);

    if numel(age) < 100
        error('At least 100 valid points are required.');
    end

    dt = median(diff(age));

    if ~isfinite(dt) || dt <= 0
        error('The age vector must increase with a positive time step.');
    end

    tolerance = max(1e-8, 1e-3 * dt);

    if any(abs(diff(age) - dt) > tolerance)
        error('The input age vector must be regularly spaced.');
    end

    data_z = safe_zscore(data);

    results = struct();
    results.age = age;
    results.data_z = data_z;
    results.dt = dt;
    results.settings = opts;

    methods = lower(string(opts.methods));

    summary_method = strings(0,1);
    summary_detected = false(0,1);
    summary_cp_age = nan(0,1);
    summary_success = false(0,1);



    %% SiZer

    if any(methods == "sizer")

        try
            results.sizer = run_sizer_method(age, data_z, opts.sizer);
            results.sizer.run_successful = true;
            results.sizer.error_message = "";

        catch ME
            results.sizer = empty_sizer_result();
            results.sizer.run_successful = false;
            results.sizer.error_message = string(ME.message);
        end

        summary_method(end+1,1) = "SiZer";
        summary_detected(end+1,1) = isfinite(results.sizer.cp_age);
        summary_cp_age(end+1,1) = results.sizer.cp_age;
        summary_success(end+1,1) = results.sizer.run_successful;

    end



    %% findchangepts

    if any(methods == "findchangepts")

        try
            results.findchangepts = ...
                run_findchangepts_method(age, data_z, dt, opts.findchangepts);

            results.findchangepts.run_successful = true;
            results.findchangepts.error_message = "";

        catch ME
            results.findchangepts = empty_findchangepts_result();
            results.findchangepts.run_successful = false;
            results.findchangepts.error_message = string(ME.message);
        end

        summary_method(end+1,1) = "findchangepts";
        summary_detected(end+1,1) = isfinite(results.findchangepts.cp_age);
        summary_cp_age(end+1,1) = results.findchangepts.cp_age;
        summary_success(end+1,1) = results.findchangepts.run_successful;

    end



    %% BEAST

    if any(methods == "beast")

        try
            results.beast = run_beast_method(age, data_z, dt, opts.beast);
            results.beast.run_successful = true;
            results.beast.error_message = "";

        catch ME
            results.beast = empty_beast_result();
            results.beast.run_successful = false;
            results.beast.error_message = string(ME.message);
        end

        summary_method(end+1,1) = "BEAST";
        summary_detected(end+1,1) = isfinite(results.beast.cp_age);
        summary_cp_age(end+1,1) = results.beast.cp_age;
        summary_success(end+1,1) = results.beast.run_successful;

    end



    %% Summary

    results.summary = table( ...
        summary_method, ...
        summary_detected, ...
        summary_cp_age, ...
        summary_success, ...
        'VariableNames', ...
        {'Method','Detected','CP_age_kyr','RunSuccessful'});

end



% DEFAULT SETTINGS

function opts = apply_defaults(opts)

    if ~isfield(opts, 'methods')
        opts.methods = ["sizer", "findchangepts", "beast"];
    end

    %% SiZer: matches the final Abram/Marron setup

    if ~isfield(opts, 'sizer')
        opts.sizer = struct();
    end

    if ~isfield(opts.sizer, 'folder')
        opts.sizer.folder = "";
    end

    if ~isfield(opts.sizer, 'bandwidth_min')
        opts.sizer.bandwidth_min = 60;
    end

    if ~isfield(opts.sizer, 'bandwidth_max')
        opts.sizer.bandwidth_max = 100;
    end

    if ~isfield(opts.sizer, 'n_bandwidths')
        opts.sizer.n_bandwidths = 41;
    end

    if ~isfield(opts.sizer, 'x_step')
        opts.sizer.x_step = 1.0;
    end

    if ~isfield(opts.sizer, 'alpha')
        opts.sizer.alpha = 0.10;
    end

    if ~isfield(opts.sizer, 'simflag')
        opts.sizer.simflag = 2;
    end

    if ~isfield(opts.sizer, 'iregtdist')
        opts.sizer.iregtdist = 0;
    end

    if ~isfield(opts.sizer, 'endpoint_flag')
        opts.sizer.endpoint_flag = 0;
    end

    if ~isfield(opts.sizer, 'boundary_adjustment')
        opts.sizer.boundary_adjustment = 0;
    end

    if ~isfield(opts.sizer, 'onset_mode')
        opts.sizer.onset_mode = "all";
    end

    if ~isfield(opts.sizer, 'min_run_kyr')
        opts.sizer.min_run_kyr = 40;
    end

    if ~isfield(opts.sizer, 'edge_ignore_kyr')
        opts.sizer.edge_ignore_kyr = 0;
    end

    if ~isfield(opts.sizer, 'cluster_gap_kyr')
        opts.sizer.cluster_gap_kyr = 35;
    end

    if ~isfield(opts.sizer, 'min_bw_fraction')
        opts.sizer.min_bw_fraction = 0.80;
    end

    if ~isfield(opts.sizer, 'uncertainty_mode')
        opts.sizer.uncertainty_mode = "iqr";
    end

    %% findchangepts

    if ~isfield(opts, 'findchangepts')
        opts.findchangepts = struct();
    end

    if ~isfield(opts.findchangepts, 'smooth_kyr')
        opts.findchangepts.smooth_kyr = 0;
    end

    if ~isfield(opts.findchangepts, 'max_num_changes')
        opts.findchangepts.max_num_changes = 1;
    end

    %% BEAST: matches the working standalone BEAST code

    if ~isfield(opts, 'beast')
        opts.beast = struct();
    end

    if ~isfield(opts.beast, 'tcp_minmax')
        opts.beast.tcp_minmax = [0 1];
    end

    if ~isfield(opts.beast, 'torder_minmax')
        opts.beast.torder_minmax = [0 1];
    end

    if ~isfield(opts.beast, 'samples')
        opts.beast.samples = 4000;
    end

    if ~isfield(opts.beast, 'burnin')
        opts.beast.burnin = 1000;
    end

    if ~isfield(opts.beast, 'chains')
        opts.beast.chains = 2;
    end

    if ~isfield(opts.beast, 'thin')
        opts.beast.thin = 5;
    end

    if ~isfield(opts.beast, 'seed')
        opts.beast.seed = 1;
    end

end



% SIZER


function out = run_sizer_method(age, data_z, settings)

    sizer_folder = resolve_sizer_folder(settings.folder);

    if exist('norminv', 'file') ~= 2
        error(['The Abram/Marron SiZer code requires norminv. ' ...
               'Enable the Statistics and Machine Learning Toolbox.']);
    end

    addpath(sizer_folder);

    if exist('sz2SM_edit', 'file') ~= 2
        error('sz2SM_edit.m could not be found after adding the SiZer folder.');
    end

    %% Run the original Abram/Marron SiZer directly

    [class_abram, class_plot, x_grid, h_grid, map90, map0] = ...
        run_original_sizer(age, data_z, settings);

    %% Detect sustained significant onsets

    [onset_age, onset_bw, onset_sign] = detect_onsets( ...
        x_grid, ...
        h_grid, ...
        class_plot, ...
        settings.onset_mode, ...
        settings.min_run_kyr, ...
        settings.edge_ignore_kyr);

    %% Cluster onsets across bandwidths

    min_bw_count = ceil( ...
        settings.min_bw_fraction * settings.n_bandwidths);

    clusters = cluster_onsets( ...
        onset_age, ...
        onset_bw, ...
        onset_sign, ...
        settings.cluster_gap_kyr, ...
        min_bw_count, ...
        settings.uncertainty_mode);

    %% Select one primary transition for Monte Carlo summaries
    % The primary cluster is the one supported by the largest number of
    % bandwidths. Ties are resolved by the number of onset points and then
    % by the older median age. All robust clusters are still returned.

    cp_age = NaN;
    q10 = NaN;
    q90 = NaN;
    primary_cluster_index = NaN;

    if ~isempty(clusters)

        support = [clusters.n_bandwidths]';
        n_points = [clusters.n_points]';
        median_age = [clusters.median_age]';

        ranking = table( ...
            (1:numel(clusters))', ...
            support, ...
            n_points, ...
            median_age, ...
            'VariableNames', ...
            {'Index','Support','NPoints','MedianAge'});

        ranking = sortrows( ...
            ranking, ...
            {'Support','NPoints','MedianAge'}, ...
            {'descend','descend','descend'});

        primary_cluster_index = ranking.Index(1);
        primary = clusters(primary_cluster_index);

        cp_age = primary.median_age;
        q10 = primary.q10_age;
        q90 = primary.q90_age;

    end

    out = struct();
    out.cp_age = cp_age;
    out.q10 = q10;
    out.q90 = q90;
    out.primary_cluster_index = primary_cluster_index;
    out.clusters = clusters;
    out.all_cp_ages = [clusters.median_age]';
    out.all_onset_ages = onset_age;
    out.all_onset_bandwidths = onset_bw;
    out.all_onset_signs = onset_sign;
    out.class_abram = class_abram;
    out.class_plot = class_plot;
    out.map_90 = map90;
    out.map_sign_only = map0;
    out.age = x_grid;
    out.bandwidths = h_grid;
    out.sizer_folder = string(sizer_folder);

end


function folder = resolve_sizer_folder(requested_folder)

    if strlength(string(requested_folder)) > 0

        requested_folder = char(string(requested_folder));

        if isfile(fullfile(requested_folder, 'sz2SM_edit.m'))
            folder = requested_folder;
            return
        end

    end

    root_dir = pwd;

    candidates = { ...
        fullfile(root_dir, 'sizer'), ...
        fullfile(root_dir, 'AbrametalPAGES2k2016_Supp3_code', 'sizer') ...
        };

    for i = 1:numel(candidates)

        if isfile(fullfile(candidates{i}, 'sz2SM_edit.m'))
            folder = candidates{i};
            return
        end

    end

    hit = dir(fullfile(root_dir, '**', 'sz2SM_edit.m'));

    if ~isempty(hit)
        folder = hit(1).folder;
        return
    end

    error(['Could not find sz2SM_edit.m. Keep the complete Abram/Marron ' ...
           'SiZer folder beside this function or provide opts.sizer.folder.']);

end


function [class_abram, class_plot, x_grid, h_grid, map90, map0] = ...
    run_original_sizer(age, y, settings)

    age = age(:);
    y = y(:);

    xmin = min(age);
    xmax = max(age);

    nbin = max( ...
        round((xmax - xmin) / settings.x_step) + 1, ...
        101);

    input_data = [age, y];

    par = struct( ...
        'vxgp', [xmin xmax nbin], ...
        'vhgp', [ ...
            settings.bandwidth_min ...
            settings.bandwidth_max ...
            settings.n_bandwidths], ...
        'eptflag', settings.endpoint_flag, ...
        'ibdryadj', settings.boundary_adjustment, ...
        'iregtdist', settings.iregtdist, ...
        'simflag', settings.simflag, ...
        'icolor', 1, ...
        'iplot', 0);

    par.alpha = settings.alpha;
    [map90, x_grid] = sz2SM_edit(input_data, par);

    par.alpha = 1;
    [map0, x_grid0] = sz2SM_edit(input_data, par);

    x_grid = x_grid(:)';
    x_grid0 = x_grid0(:)';

    if numel(x_grid) ~= numel(x_grid0) || ...
            any(abs(x_grid - x_grid0) > 1e-10)

        error('The alpha=0.10 and alpha=1 SiZer age grids differ.');

    end

    h_grid = linspace( ...
        settings.bandwidth_min, ...
        settings.bandwidth_max, ...
        settings.n_bandwidths);

    if size(map90,1) ~= numel(h_grid) && ...
            size(map90,2) == numel(h_grid)
        map90 = map90';
    end

    if size(map0,1) ~= numel(h_grid) && ...
            size(map0,2) == numel(h_grid)
        map0 = map0';
    end

    if size(map90,1) ~= numel(h_grid) || ...
            size(map90,2) ~= numel(x_grid)

        error('Unexpected output dimensions from sz2SM_edit.m.');

    end

    if ~isequal(size(map90), size(map0))
        error('The alpha=0.10 and alpha=1 SiZer maps differ in size.');
    end

    %% Original Abram combination
    % 1 = significant positive
    % 2 = non-significant positive
    % 3 = non-significant negative
    % 4 = significant negative
    % -999 = insufficient information

    z90 = double(map90);
    z0 = double(map0);

    z90(z90 == 2) = -999;
    z0(z0 == 2) = -999;

    class_abram = (z0 + z90) / 2;
    class_abram(abs(class_abram - 3.5) < 1e-12) = 3;
    class_abram(class_abram < 0) = -999;

    %% Remap to final plotting key
    % 1 = significant negative
    % 2 = non-significant negative
    % 3 = non-significant positive
    % 4 = significant positive
    % 5 = insufficient information

    class_plot = 5 * ones(size(class_abram));
    ok = class_abram >= 1 & class_abram <= 4;
    class_plot(ok) = 5 - class_abram(ok);

end


function [ages, bws, signs] = detect_onsets( ...
    x, h, C, mode, min_run, edge_ignore)

    ages = [];
    bws = [];
    signs = [];

    young_limit = min(x) + edge_ignore;
    old_limit = max(x) - edge_ignore;

    for ih = 1:numel(h)

        row = C(ih,:);

        neg = get_runs( ...
            x, row == 1, min_run, young_limit, old_limit);

        pos = get_runs( ...
            x, row == 4, min_run, young_limit, old_limit);

        switch lower(string(mode))

            case "negative_only"
                [ages, bws, signs] = append_runs( ...
                    ages, bws, signs, neg, h(ih), -1);

            case "positive_only"
                [ages, bws, signs] = append_runs( ...
                    ages, bws, signs, pos, h(ih), 1);

            case "all"
                [ages, bws, signs] = append_runs( ...
                    ages, bws, signs, neg, h(ih), -1);

                [ages, bws, signs] = append_runs( ...
                    ages, bws, signs, pos, h(ih), 1);

            otherwise
                error('Unknown SiZer onset mode: %s', char(string(mode)));

        end

    end

end


function [ages, bws, signs] = append_runs( ...
    ages, bws, signs, runs, bandwidth, sign_value)

    for k = 1:numel(runs)

        ages(end+1,1) = runs(k).old_edge_age;
        bws(end+1,1) = bandwidth;
        signs(end+1,1) = sign_value;

    end

end


function runs = get_runs( ...
    age, mask, min_run, young_limit, old_limit)

    age = age(:);
    mask = logical(mask(:)');

    runs = struct( ...
        'young_edge_age', {}, ...
        'old_edge_age', {}, ...
        'duration', {});

    if isempty(mask) || ~any(mask)
        return
    end

    difference = diff([0 mask 0]);
    start_indices = find(difference == 1);
    end_indices = find(difference == -1) - 1;

    for k = 1:numel(start_indices)

        young = age(start_indices(k));
        old = age(end_indices(k));
        duration = old - young;

        if old <= young_limit || old >= old_limit
            continue
        end

        if duration >= min_run

            runs(end+1).young_edge_age = young; 
            runs(end).old_edge_age = old;
            runs(end).duration = duration;

        end

    end

end


function clusters = cluster_onsets( ...
    ages, bws, signs, gap, min_bws, uncertainty_mode)

    clusters = struct( ...
        'median_age', {}, ...
        'q10_age', {}, ...
        'q90_age', {}, ...
        'unc_low_age', {}, ...
        'unc_high_age', {}, ...
        'sign', {}, ...
        'n_points', {}, ...
        'n_bandwidths', {}, ...
        'ages', {}, ...
        'bandwidths', {});

    if isempty(ages)
        return
    end

    sign_values = unique(signs(isfinite(signs)));

    for ss = 1:numel(sign_values)

        this_sign = sign_values(ss);

        use = signs == this_sign & isfinite(ages) & isfinite(bws);

        a = ages(use);
        b = bws(use);

        if isempty(a)
            continue
        end

        [a, order] = sort(a);
        b = b(order);

        current_a = a(1);
        current_b = b(1);

        for i = 2:numel(a)

            if a(i) - a(i-1) <= gap

                current_a(end+1,1) = a(i); 
                current_b(end+1,1) = b(i); 

            else

                clusters = add_cluster( ...
                    clusters, ...
                    current_a, ...
                    current_b, ...
                    this_sign, ...
                    min_bws, ...
                    uncertainty_mode);

                current_a = a(i);
                current_b = b(i);

            end

        end

        clusters = add_cluster( ...
            clusters, ...
            current_a, ...
            current_b, ...
            this_sign, ...
            min_bws, ...
            uncertainty_mode);

    end

    if ~isempty(clusters)

        [~, order] = sort([clusters.median_age], 'descend');
        clusters = clusters(order);

    end

end


function clusters = add_cluster( ...
    clusters, ages, bws, sign_value, min_bws, uncertainty_mode)

    unique_bws = unique(bws);

    if numel(unique_bws) < min_bws
        return
    end

    q10 = local_percentile(ages, 10);
    q90 = local_percentile(ages, 90);

    switch lower(string(uncertainty_mode))

        case "iqr"
            low_age = q10;
            high_age = q90;

        case "range"
            low_age = min(ages, [], 'omitnan');
            high_age = max(ages, [], 'omitnan');

        otherwise
            error('Unknown uncertainty mode: %s', ...
                char(string(uncertainty_mode)));

    end

    clusters(end+1).median_age = median(ages, 'omitnan');
    clusters(end).q10_age = q10;
    clusters(end).q90_age = q90;
    clusters(end).unc_low_age = low_age;
    clusters(end).unc_high_age = high_age;
    clusters(end).sign = sign_value;
    clusters(end).n_points = numel(ages);
    clusters(end).n_bandwidths = numel(unique_bws);
    clusters(end).ages = ages;
    clusters(end).bandwidths = bws;

end


function value = local_percentile(x, percentile)

    x = sort(x(isfinite(x)));

    if isempty(x)
        value = NaN;
        return
    end

    x = x(:);
    n = numel(x);

    if n == 1
        value = x;
        return
    end

    position = 1 + (percentile / 100) * (n - 1);
    low_index = floor(position);
    high_index = ceil(position);

    if low_index == high_index
        value = x(low_index);
    else
        weight = position - low_index;
        value = (1 - weight) * x(low_index) + weight * x(high_index);
    end

end



% FINDCHANGEPTS


function out = run_findchangepts_method(age, data_z, dt, settings)

    if exist('findchangepts', 'file') ~= 2
        error(['findchangepts is unavailable. ' ...
               'The Signal Processing Toolbox is required.']);
    end

    smooth_points = round(settings.smooth_kyr / dt);
    smooth_points = max(smooth_points, 1);

    if mod(smooth_points, 2) == 0
        smooth_points = smooth_points + 1;
    end

    detection_input = movmean( ...
        data_z, ...
        smooth_points, ...
        'omitnan', ...
        'Endpoints', ...
        'shrink');

    change_index = findchangepts( ...
        detection_input, ...
        'Statistic', ...
        'linear', ...
        'MaxNumChanges', ...
        settings.max_num_changes);

    if isempty(change_index)
        cp_age = NaN;
    else
        cp_age = age(change_index(1));
    end

    out = struct();
    out.cp_age = cp_age;
    out.change_index = change_index;
    out.detection_input = detection_input;
    out.smooth_kyr = settings.smooth_kyr;
    out.smooth_points = smooth_points;

end



% BEAST


function result = run_beast_method(age, data_z, dt, settings)

    if exist('beast', 'file') ~= 2
        error('The BEAST function was not found on the MATLAB path.');
    end

    beast_output = beast( ...
        data_z(:), ...
        'start', age(1), ...
        'deltat', dt, ...
        'season', 'none', ...
        'tcp.minmax', settings.tcp_minmax, ...
        'torder.minmax', settings.torder_minmax, ...
        'mcmc.seed', settings.seed, ...
        'mcmc.samples', settings.samples, ...
        'mcmc.burnin', settings.burnin, ...
        'mcmc.chains', settings.chains, ...
        'mcmc.thin', settings.thin, ...
        'print.progress', false, ...
        'print.options', false);

    ncp_mode = NaN;
    probability_0cp = NaN;
    probability_1cp = NaN;
    cp_age = NaN;
    cp_probability = NaN;
    trend = nan(size(data_z));
    cp_occurrence_probability = nan(size(data_z));

    if isfield(beast_output, 'trend')

        beast_trend = beast_output.trend;

        %% Most probable number of trend changepoints

        if isfield(beast_trend, 'ncp_mode') && ...
                ~isempty(beast_trend.ncp_mode)

            ncp_mode = beast_trend.ncp_mode(1);

        end

        %% Probability of 0 and 1 changepoints

        if isfield(beast_trend, 'ncpPr') && ...
                ~isempty(beast_trend.ncpPr)

            ncp_probability = beast_trend.ncpPr(:);

            if numel(ncp_probability) >= 1
                probability_0cp = ncp_probability(1);
            end

            if numel(ncp_probability) >= 2
                probability_1cp = ncp_probability(2);
            end

            % Important fallback used in the working standalone code.
            if ~isfinite(ncp_mode)
                [~, mode_index] = max(ncp_probability);
                ncp_mode = mode_index - 1;
            end

        end

        %% Extract CP only when BEAST selected at least one CP

        if isfinite(ncp_mode) && ncp_mode >= 1

            if isfield(beast_trend, 'cp') && ...
                    ~isempty(beast_trend.cp)

                cp_candidates = beast_trend.cp(:);

                if isfield(beast_trend, 'cpPr') && ...
                        ~isempty(beast_trend.cpPr)

                    cp_probabilities = beast_trend.cpPr(:);

                else

                    cp_probabilities = nan(size(cp_candidates));

                end

                valid_cp = ...
                    isfinite(cp_candidates) & ...
                    cp_candidates >= min(age) & ...
                    cp_candidates <= max(age);

                cp_candidates = cp_candidates(valid_cp);
                cp_probabilities = cp_probabilities(valid_cp);

                if ~isempty(cp_candidates)

                    if any(isfinite(cp_probabilities))

                        finite_probability = isfinite(cp_probabilities);

                        finite_candidates = ...
                            cp_candidates(finite_probability);

                        finite_cp_probability = ...
                            cp_probabilities(finite_probability);

                        [cp_probability, best_index] = ...
                            max(finite_cp_probability);

                        cp_age = finite_candidates(best_index);

                    else

                        cp_age = cp_candidates(1);

                    end

                end

            end

            %% Fallback: maximum CP occurrence probability

            if ~isfinite(cp_age) && ...
                    isfield(beast_trend, 'cpOccPr') && ...
                    ~isempty(beast_trend.cpOccPr)

                cp_occurrence = beast_trend.cpOccPr(:);

                if numel(cp_occurrence) == numel(age)
                    [cp_probability, best_index] = max(cp_occurrence);
                    cp_age = age(best_index);
                end

            end

        end

        %% Store fitted BEAST trend

        if isfield(beast_trend, 'Y') && ...
                ~isempty(beast_trend.Y)

            candidate_trend = beast_trend.Y(:);

            if numel(candidate_trend) == numel(data_z)
                trend = candidate_trend;
            end

        end

        %% Store CP occurrence probability

        if isfield(beast_trend, 'cpOccPr') && ...
                ~isempty(beast_trend.cpOccPr)

            candidate_probability = beast_trend.cpOccPr(:);

            if numel(candidate_probability) == numel(data_z)
                cp_occurrence_probability = candidate_probability;
            end

        end

    end

    result = struct();
    result.cp_age = cp_age;
    result.cp_probability = cp_probability;
    result.ncp_mode = ncp_mode;
    result.probability_0cp = probability_0cp;
    result.probability_1cp = probability_1cp;
    result.trend = trend;
    result.cp_occurrence_probability = cp_occurrence_probability;

end



% EMPTY RESULTS


function out = empty_sizer_result()

    out = struct( ...
        'cp_age', NaN, ...
        'q10', NaN, ...
        'q90', NaN, ...
        'primary_cluster_index', NaN, ...
        'clusters', [], ...
        'all_cp_ages', [], ...
        'all_onset_ages', [], ...
        'all_onset_bandwidths', [], ...
        'all_onset_signs', [], ...
        'class_abram', [], ...
        'class_plot', [], ...
        'map_90', [], ...
        'map_sign_only', [], ...
        'age', [], ...
        'bandwidths', [], ...
        'sizer_folder', "");

end


function out = empty_findchangepts_result()

    out = struct( ...
        'cp_age', NaN, ...
        'change_index', [], ...
        'detection_input', [], ...
        'smooth_kyr', NaN, ...
        'smooth_points', NaN);

end


function out = empty_beast_result()

    out = struct( ...
        'cp_age', NaN, ...
        'cp_probability', NaN, ...
        'ncp_mode', NaN, ...
        'probability_0cp', NaN, ...
        'probability_1cp', NaN, ...
        'trend', [], ...
        'cp_occurrence_probability', []);

end



% SAFE Z-SCORE


function z = safe_zscore(x)

    x_mean = mean(x, 'omitnan');
    x_std = std(x, 0, 'omitnan');

    if ~isfinite(x_std) || x_std == 0
        z = x - x_mean;
    else
        z = (x - x_mean) ./ x_std;
    end

end
