function epp_plot_population_psth(cfg)
%EPP_PLOT_POPULATION_PSTH  Population mean PSTH from unit rasters.
%
%   epp_plot_population_psth(cfg)
%
%   Inputs:
%     cfg - from epp_load_cfg. Loads every *_raster.mat under roots.raster
%           (Neuron_ID, windows). Builds condition means from raster +
%           cond_mask (same pad/smooth/crop as epp_plot_psth).
%
%   Outputs:
%     cfg.roots.population_psth/population_psth.png
%     population_psth.mat (pop_data): n_units, Neuron_ID, windows(wi).conditions
%     with mean/sem across units (not trials).

epp_ensure_paths();
unit_rasters = load_unit_raster_windows(cfg.roots.raster);
if isempty(unit_rasters)
    warning('epp_plot_population_psth:NoData', 'No unit raster files found in %s.', cfg.roots.raster);
    return;
end

n_windows = numel(unit_rasters(1).windows);
if n_windows == 0
    return;
end

n_conditions = numel(cfg.CONDITIONS);
conv_kernel = psth_conv_kernel(cfg);
pop_data = struct();
pop_data.n_units = numel(unit_rasters);
pop_data.Neuron_ID = {unit_rasters.Neuron_ID};
pop_data.windows = repmat(struct( ...
    'name', '', 'align_state', [], 'conditions', []), n_windows, 1);

fig = figure('Color', 'w');
for wi = 1:n_windows
    w_ref = unit_rasters(1).windows(wi);
    pop_data.windows(wi).name = w_ref.name;
    pop_data.windows(wi).align_state = w_ref.align_state;
    pop_data.windows(wi).conditions = repmat(struct( ...
        'name', '', 'n_units', 0, 'n_trials_total', 0, ...
        'bin_centers_s', [], 'mean_rate_hz', [], 'sem_rate_hz', []), n_conditions, 1);

    ax = subplot(n_windows, 1, wi);
    hold(ax, 'on');

    for ci = 1:n_conditions
        [unit_means, bin_centers, n_trials_vec] = collect_condition_means( ...
            unit_rasters, wi, ci, cfg, conv_kernel);
        if isempty(unit_means)
            continue;
        end

        pop_mean = mean(unit_means, 1, 'omitnan');
        n_units_ci = size(unit_means, 1);
        pop_sem = std(unit_means, 0, 1, 'omitnan') ./ sqrt(max(n_units_ci, 1));

        pop_data.windows(wi).conditions(ci).name = cfg.CONDITIONS(ci).name;
        pop_data.windows(wi).conditions(ci).n_units = n_units_ci;
        pop_data.windows(wi).conditions(ci).n_trials_total = sum(n_trials_vec);
        pop_data.windows(wi).conditions(ci).bin_centers_s = bin_centers;
        pop_data.windows(wi).conditions(ci).mean_rate_hz = pop_mean;
        pop_data.windows(wi).conditions(ci).sem_rate_hz = pop_sem;

        color = normalize_plot_color(cfg.CONDITIONS(ci).color);
        lineProps = {'LineWidth', 1.5, 'Color', color};
        shadedErrorBar(bin_centers(:).', pop_mean(:).', pop_sem(:).', lineProps, 0);
    end

    hold(ax, 'off');
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Rate (Hz)');
    title(ax, sprintf('%s (%d)', w_ref.name, w_ref.align_state));
    legend(ax, {cfg.CONDITIONS.name}, 'Interpreter', 'none', 'Location', 'best');
end

saveas(fig, fullfile(cfg.roots.population_psth, 'population_psth.png'));
close(fig);
save(fullfile(cfg.roots.population_psth, 'population_psth.mat'), 'pop_data');
end

function unit_rasters = load_unit_raster_windows(raster_root)
unit_rasters = struct([]);
raster_files = dir(fullfile(raster_root, '*_raster.mat'));
for i = 1:numel(raster_files)
    loaded = load(fullfile(raster_root, raster_files(i).name), 'Neuron_ID', 'windows');
    if ~isfield(loaded, 'windows') || isempty(loaded.windows)
        continue
    end
    entry = struct('Neuron_ID', loaded.Neuron_ID, 'windows', loaded.windows);
    unit_rasters = [unit_rasters; entry]; %#ok<AGROW>
end
end

function [unit_means, bin_centers, n_trials_vec] = collect_condition_means( ...
    unit_rasters, wi, ci, cfg, conv_kernel)
unit_means = [];
bin_centers = [];
n_trials_vec = [];

for u = 1:numel(unit_rasters)
    if numel(unit_rasters(u).windows) < wi
        continue
    end
    w = unit_rasters(u).windows(wi);
    [mean_rate, centers, n_trials] = unit_condition_psth(w, ci, cfg, conv_kernel);
    if isempty(mean_rate)
        continue
    end
    if isempty(bin_centers)
        bin_centers = centers(:).';
    elseif numel(centers) ~= numel(bin_centers) || any(centers(:) ~= bin_centers(:))
        continue
    end
    unit_means = [unit_means; mean_rate(:).']; %#ok<AGROW>
    n_trials_vec = [n_trials_vec; n_trials]; %#ok<AGROW>
end
end

function [mean_rate, psth_centers, n_trials] = unit_condition_psth(w, ci, cfg, conv_kernel)
% Same pad -> smooth -> crop -> psth-bin path as epp_plot_psth.
mean_rate = [];
psth_centers = [];
n_trials = 0;
if ~isfield(w, 'raster') || isempty(w.raster) || ~isfield(w, 'cond_mask') ...
        || isempty(w.cond_mask) || ci > size(w.cond_mask, 2)
    return
end
keep_mask = w.cond_mask(:, ci);
if ~any(keep_mask)
    return
end

core_cols = (w.psth_pad_bins + 1):(size(w.raster, 2) - w.psth_pad_bins);
psth_edges = w.t_start_s:cfg.psth.bin_size_s:w.t_end_s;
if numel(psth_edges) < 2 || psth_edges(end) < w.t_end_s
    psth_edges(end+1) = w.t_end_s; %#ok<AGROW>
end
psth_centers = psth_edges(1:end-1).' + diff(psth_edges(1:2))/2;
n_psth_bins = numel(psth_centers);
raster_to_psth_bin = discretize(w.bin_centers(core_cols), psth_edges);

cond_rate = double(w.raster(keep_mask, :)) ./ w.bin_size_s;
cond_rate = conv2(cond_rate, conv_kernel, 'same');
cond_rate = cond_rate(:, core_cols);
cond_rate_psth = NaN(size(cond_rate, 1), n_psth_bins);
for bi = 1:n_psth_bins
    src_cols = raster_to_psth_bin == bi;
    if any(src_cols)
        cond_rate_psth(:, bi) = mean(cond_rate(:, src_cols), 2, 'omitnan');
    end
end
mean_rate = mean(cond_rate_psth, 1, 'omitnan');
n_trials = nnz(keep_mask);
end

function conv_kernel = psth_conv_kernel(cfg)
bin_size_s = cfg.raster_bin_size_s;
switch lower(char(cfg.psth.smoothing_kernel))
    case 'gaussian'
        sigma_bins = max(cfg.psth.smoothing_width_s / bin_size_s, eps);
        radius_bins = max(1, ceil(4 * sigma_bins));
        x = -radius_bins:radius_bins;
        conv_kernel = exp(-0.5 * (x / sigma_bins).^2);
        conv_kernel = conv_kernel ./ sum(conv_kernel);
    case 'box'
        box_bins = max(1, round(cfg.psth.smoothing_width_s / bin_size_s));
        conv_kernel = ones(1, box_bins) / box_bins;
    otherwise
        conv_kernel = 1;
end
end

function color = normalize_plot_color(color_in)
color = color_in(:).';
if numel(color) ~= 3
    color = [0 0 0];
end
if max(color) > 1
    color = color / 255;
end
end
