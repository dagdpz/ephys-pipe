function epp_plot_psth(cfg, unit_row)
%EPP_PLOT_PSTH  Stage 4: condition-split PSTH and binned raster for one unit.
%
%   epp_plot_psth(cfg, unit_row)
%
%   Inputs:
%     cfg - struct from epp_load_cfg.
%     unit_row - one element of unit_info. Reads
%           cfg.roots.raster/<Neuron_ID>_raster.mat. Returns without writing
%           if that file is missing.
%           Uses cfg.CONDITIONS (name, color, parameters) and cfg.psth
%           (bin_size_s, smoothing_kernel 'gaussian'|'box', smoothing_width_s).
%           Condition parameters are trial-field names; a trial matches a
%           field if its value is in the list, and all fields are ANDed.
%           The same trial may match more than one condition.
%
%   Outputs:
%     none. For each unit that has a raster file, writes
%       cfg.roots.psth/<Neuron_ID>_psth.png
%       cfg.roots.psth/<Neuron_ID>_psth.mat  (variable psth_data)
%     psth_data.Neuron_ID, psth_data.windows (raster windows plus .conditions).
%     conditions(ci): name, n_trials, bin_centers_s, mean_rate_hz, sem_rate_hz.
%     Windows are concatenated on one axes (0.1 s gap). X tick labels are
%     time within each window. PSTHs sit above 0, rasters below, then eye
%     (red) and hand traces (left blue, right green). Raster squares are one
%     per occupied bin, colored by the first matching condition (gray if none).
%
% Stage 4: load rasters, condition-split, and plot PSTH + rasters.

Neuron_ID = unit_row.Neuron_ID;
raster_file = fullfile(cfg.roots.raster, sprintf('%s_raster.mat', Neuron_ID));
if ~isfile(raster_file)
    return
end
loaded = load(raster_file, 'raster_data');
raster_data = loaded.raster_data;

    psth_data = struct();
    psth_data.Neuron_ID = Neuron_ID;
    psth_data.windows = raster_data.windows;

    n_windows = numel(raster_data.windows);
    if n_windows == 0
        return
    end
    epp_ensure_paths();
    fig = figure('units','normalized','outerposition',[0 0 1 1],'color','w');
    ax = axes(fig);
    hold(ax, 'on');
    raster_trial_height = 0.1;
    window_gap_s = 0.1;
    state_separator = 0;
    y_max = 1;
    y_min = 0;
    align_x = zeros(n_windows, 1);
    label_x = zeros(n_windows, 1);
    label_name = cell(n_windows, 1);
    legend_h = gobjects(numel(cfg.CONDITIONS), 1);

    for wi = 1:n_windows
        w = raster_data.windows(wi);
        n_events = size(w.raster, 1);
        psth_edges = w.t_start_s:cfg.psth.bin_size_s:w.t_end_s;
        if psth_edges(end) < w.t_end_s
            psth_edges(end+1) = w.t_end_s; %#ok<AGROW>
        end
        psth_centers = psth_edges(1:end-1).' + diff(psth_edges(1:2))/2;
        raster_to_psth_bin = discretize(w.bin_centers, psth_edges);
        n_psth_bins = numel(psth_centers);
        raster_bin_size_s = w.bin_size_s;

        kernel_name = lower(char(cfg.psth.smoothing_kernel));
        switch kernel_name
            case 'gaussian'
                sigma_bins = max(cfg.psth.smoothing_width_s / raster_bin_size_s, eps);
                radius_bins = max(1, ceil(4 * sigma_bins));
                x = -radius_bins:radius_bins;
                conv_kernel = exp(-0.5 * (x / sigma_bins).^2);
                conv_kernel = conv_kernel ./ sum(conv_kernel);
            case 'box'
                box_bins = max(1, round(cfg.psth.smoothing_width_s / raster_bin_size_s));
                conv_kernel = ones(1, box_bins) / box_bins;
        end

        conditions_out = repmat(struct('name', '', 'n_trials', 0, 'bin_centers_s', [], 'mean_rate_hz', [], 'sem_rate_hz', []), numel(cfg.CONDITIONS), 1);
        cond_masks = false(n_events, numel(cfg.CONDITIONS));
        cond_colors = zeros(numel(cfg.CONDITIONS), 3);
        assigned_mask = false(n_events,1);

        % Same seconds-per-pixel in every window. 0.1 s gap is the separator.
        state_shift = state_separator - w.t_start_s;
        state_separator = state_shift + w.t_end_s + window_gap_s;
        align_x(wi) = state_shift;
        label_x(wi) = state_shift + (w.t_start_s + w.t_end_s)/2;
        label_name{wi} = w.name;
        x_win = state_shift + [w.t_start_s w.t_end_s];

        % No align events: do not read trial fields. Empty rates are skipped
        % by the population mean (isempty(mean_rate_hz)).
        if n_events == 0 || isempty(w.trial)
            for ci = 1:numel(cfg.CONDITIONS)
                conditions_out(ci).name = cfg.CONDITIONS(ci).name;
                conditions_out(ci).n_trials = 0;
                conditions_out(ci).bin_centers_s = psth_centers;
                conditions_out(ci).mean_rate_hz = [];
                conditions_out(ci).sem_rate_hz = [];
            end
            psth_data.windows(wi).conditions = conditions_out;
            continue;
        end

        for ci = 1:numel(cfg.CONDITIONS)
            cond_def = cfg.CONDITIONS(ci);
            keep_mask = true(n_events, 1);
            param_names = fieldnames(cond_def.parameters);
            for pi = 1:numel(param_names)
                pn = param_names{pi};
                target_vals = cond_def.parameters.(pn);
                trial_vals = [w.trial.(pn)].';
                keep_mask = keep_mask & ismember(trial_vals, target_vals(:));
            end
            cond_masks(:, ci) = keep_mask;

            cond_counts = w.raster(keep_mask, :);
            cond_rate_raster = cond_counts ./ raster_bin_size_s;
            cond_rate_raster = conv2(cond_rate_raster, conv_kernel, 'same');
            cond_rate_psth = NaN(size(cond_counts, 1), n_psth_bins);
            for bi = 1:n_psth_bins
                src_cols = raster_to_psth_bin == bi;
                if any(src_cols)
                    cond_rate_psth(:, bi) = mean(cond_rate_raster(:, src_cols), 2, 'omitnan');
                end
            end
            cond_rate = cond_rate_psth;
            mean_rate = mean(cond_rate, 1, 'omitnan');
            sem_rate = std(cond_rate, 0, 1, 'omitnan') ./ sqrt(max(size(cond_rate,1), 1));
            conditions_out(ci).name = cond_def.name;
            conditions_out(ci).n_trials = size(cond_rate, 1);
            conditions_out(ci).bin_centers_s = psth_centers;
            conditions_out(ci).mean_rate_hz = mean_rate;
            conditions_out(ci).sem_rate_hz = sem_rate;

            color = cond_def.color;
            if max(color) > 1, color = color / 255; end
            cond_colors(ci,:) = color;
            if isempty(mean_rate) || ~any(isfinite(mean_rate))
                continue
            end
            peak = max(mean_rate + sem_rate, [], 'omitnan');
            if isfinite(peak)
                y_max = max(y_max, peak);
            end
            H = shadedErrorBar(psth_centers(:).' + state_shift, mean_rate(:).', sem_rate(:).', {'color', color, 'linewidth', 1.5}, 1);
            if ~isgraphics(legend_h(ci))
                if isstruct(H) && isfield(H, 'mainLine')
                    legend_h(ci) = H.mainLine;
                elseif isgraphics(H)
                    legend_h(ci) = H;
                end
            end
        end

        raster_order = [];
        raster_group = [];
        for ci = 1:numel(cfg.CONDITIONS)
            idx_ci = find(cond_masks(:,ci) & ~assigned_mask);
            raster_order = [raster_order; idx_ci(:)]; %#ok<AGROW>
            raster_group = [raster_group; ci*ones(numel(idx_ci),1)]; %#ok<AGROW>
            assigned_mask(idx_ci) = true;
        end
        idx_unassigned = find(~assigned_mask);
        raster_order = [raster_order; idx_unassigned(:)];
        raster_group = [raster_group; zeros(numel(idx_unassigned),1)];

        raster_y = 0;
        for ci = [1:numel(cfg.CONDITIONS), 0]
            rows = find(raster_group == ci);
            n_rows = numel(rows);
            if n_rows == 0
                continue
            end
            if ci > 0
                plot_color = cond_colors(ci,:);
            else
                plot_color = [0.5 0.5 0.5];
            end
            line(x_win, [raster_y raster_y], 'Color', plot_color, 'LineWidth', 0.5);
            xs = [];
            ys = [];
            for r = 1:n_rows
                spike_cols = find(w.raster(raster_order(rows(r)),:) > 0);
                if isempty(spike_cols)
                    continue
                end
                xs = [xs, w.bin_centers(spike_cols(:)).' + state_shift]; %#ok<AGROW>
                ys = [ys, (raster_y - r*raster_trial_height) * ones(1, numel(spike_cols))]; %#ok<AGROW>
            end
            if ~isempty(xs)
                scatter(ax, xs, ys, 10, plot_color, 's', 'filled');
            end
            raster_y = raster_y - n_rows*raster_trial_height;
            line(x_win, [raster_y raster_y], 'Color', plot_color, 'LineWidth', 0.5);
            if wi == 1
                text(x_win(1), raster_y + n_rows*raster_trial_height, sprintf('N=%d', n_rows), 'Color', plot_color, 'FontSize', 5, 'VerticalAlignment', 'top');
            end
        end
        y_min = min(y_min, raster_y);

        psth_data.windows(wi).conditions = conditions_out;
    end

    % Eye, then hand, under the lowest raster. x and y share a color.
    trace_scale = 0.5;
    band = 20 * trace_scale;
    gap = 1;
    eye_center = y_min - gap - band/2;
    hand_center = eye_center - gap - band;
    y_min = hand_center - band/2;
    for wi = 1:n_windows
        w = raster_data.windows(wi);
        if ~isfield(w, 'behavior') || isempty(w.behavior) || isempty(w.trial)
            continue
        end
        xx = w.bin_edges(:) + align_x(wi);
        for k = 1:numel(w.behavior)
            line(ax, xx, w.behavior(k).x_eye(:) * trace_scale + eye_center, 'Color', [1 0 0], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            line(ax, xx, w.behavior(k).y_eye(:) * trace_scale + eye_center, 'Color', [1 0 0], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            if ~isfield(w.trial, 'reach_hand')
                continue
            end
            if w.trial(k).reach_hand == 1
                hcol = [0 0 1];
            elseif w.trial(k).reach_hand == 2
                hcol = [0 1 0];
            else
                continue
            end
            line(ax, xx, w.behavior(k).x_hnd(:) * trace_scale + hand_center, 'Color', hcol, 'LineWidth', 0.5, 'HandleVisibility', 'off');
            line(ax, xx, w.behavior(k).y_hnd(:) * trace_scale + hand_center, 'Color', hcol, 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end

    xlim(ax, [0, state_separator - window_gap_s]);
    ylim(ax, [y_min, y_max]);
    xtick = [];
    xticklabel = {};
    tick_step = 0.2;
    for wi = 1:n_windows
        w = raster_data.windows(wi);
        rel = ceil(w.t_start_s/tick_step)*tick_step : tick_step : floor(w.t_end_s/tick_step)*tick_step;
        rel = round(rel/tick_step)*tick_step;
        if w.t_start_s <= 0 && w.t_end_s >= 0
            rel = unique([rel, 0]);
        end
        xtick = [xtick, align_x(wi) + rel]; %#ok<AGROW>
        xticklabel = [xticklabel, arrayfun(@(t) sprintf('%g', t), rel, 'UniformOutput', false)]; %#ok<AGROW>
    end
    set(ax, 'XTick', xtick, 'XTickLabel', xticklabel);
    for wi = 1:n_windows
        line(ax, [align_x(wi) align_x(wi)], [y_min y_max], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
        text(ax, label_x(wi), y_max - 0.05*y_max, label_name{wi}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Interpreter', 'none');
    end
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Rate (Hz)');
    title(ax, Neuron_ID, 'Interpreter', 'none');
    if any(isgraphics(legend_h))
        legend(ax, legend_h(isgraphics(legend_h)), {cfg.CONDITIONS(isgraphics(legend_h)).name}, 'Interpreter', 'none', 'Location', 'northeast');
    end
    hold(ax, 'off');
    saveas(fig, fullfile(cfg.roots.psth, sprintf('%s_psth.png', Neuron_ID)));
    close(fig);
    save(fullfile(cfg.roots.psth, sprintf('%s_psth.mat', Neuron_ID)), 'psth_data');
end
