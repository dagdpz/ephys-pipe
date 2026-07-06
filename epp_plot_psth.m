function epp_plot_psth(cfg)
% Stage 3: load rasters, condition-split, and plot PSTH + rasters.

load(fullfile(cfg.roots.processed_trials, 'unit_info.mat'), 'unit_info');

raster_files = dir(fullfile(cfg.roots.raster, '*_raster.mat'));
all_rasters = struct([]);
for i = 1:numel(raster_files)
    loaded = load(fullfile(cfg.roots.raster, raster_files(i).name), 'raster_data');
    all_rasters = [all_rasters; loaded.raster_data]; %#ok<AGROW>
end

for u = 1:numel(unit_info)
    Neuron_ID = unit_info(u).Neuron_ID;
    idx_raster = find(strcmp({all_rasters.Neuron_ID}, Neuron_ID), 1, 'first');
    if isempty(idx_raster), continue; end
    raster_data = all_rasters(idx_raster);

    psth_data = struct();
    psth_data.Neuron_ID = Neuron_ID;
    psth_data.windows = raster_data.windows;

    n_windows = numel(raster_data.windows);
    if n_windows == 0, continue; end
    fig = figure('Color', 'w');
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

        ax_raster = subplot(n_windows, 2, (wi-1)*2 + 1);
        hold(ax_raster, 'on');
        ax = subplot(n_windows, 2, (wi-1)*2 + 2);
        hold(ax, 'on');
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
            plot(ax, psth_centers, mean_rate, 'LineWidth', 1.5, 'Color', color);
        end

        xlabel(ax, 'Time (s)');
        ylabel(ax, 'Rate (Hz)');
        title(ax, sprintf('%s (%d)', w.name, w.align_state));
        legend(ax, {cfg.CONDITIONS.name}, 'Interpreter', 'none', 'Location', 'best');
        hold(ax, 'off');

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

        for rr = 1:numel(raster_order)
            src_idx = raster_order(rr);
            spike_cols = find(w.raster(src_idx,:) > 0);
            if isempty(spike_cols), continue; end
            x_vals = w.bin_centers(spike_cols);
            y_vals = rr * ones(numel(spike_cols),1);
            ci = raster_group(rr);
            if ci > 0
                plot_color = cond_colors(ci,:);
            else
                plot_color = [0.5 0.5 0.5];
            end
            plot(ax_raster, x_vals, y_vals, '.', 'Color', plot_color, 'MarkerSize', 6);
        end
        set(ax_raster, 'YDir', 'reverse');
        xlabel(ax_raster, 'Time (s)');
        ylabel(ax_raster, 'Trials');
        title(ax_raster, sprintf('%s raster', w.name));
        xlim(ax_raster, [w.t_start_s w.t_end_s]);
        hold(ax_raster, 'off');

        psth_data.windows(wi).conditions = conditions_out;
    end

    %sgtitle(fig, sprintf('Unit %s', Neuron_ID), 'Interpreter', 'none');
    saveas(fig, fullfile(cfg.roots.psth, sprintf('%s_psth.png', Neuron_ID)));
    close(fig);
    save(fullfile(cfg.roots.psth, sprintf('%s_psth.mat', Neuron_ID)), 'psth_data');
end
end
