function epp_build_rasters(cfg, unit)
%EPP_BUILD_RASTERS  Stage 2: one unit's event-aligned raster.
%
%   epp_build_rasters(cfg, unit)
%
%   Inputs:
%     cfg  - struct from epp_load_cfg. Uses cfg.WINDOWS, cfg.EPOCHS,
%            cfg.CONDITIONS (baked into cond_mask), cfg.psth (pad width),
%            cfg.raster_bin_size_s, cfg.roots.raster.
%     unit - one unit from epp_load_unit_combined.
%
%   Outputs:
%     none. Writes cfg.roots.raster/<Neuron_ID>_raster.mat with top-level
%     variables so callers can load selectively:
%       Neuron_ID, windows, behavior, epochs
%
%   windows  - cfg.WINDOWS plus, per window:
%       bin_edges, bin_centers, bin_size_s, psth_pad_bins, event_times,
%       trial_n, cond_mask (n_events x n_conditions logical), raster
%       (logical, spike present in bin). Multi-spike bins collapse to true.
%       Changing cfg.CONDITIONS requires a rebuild.
%   behavior - parallel to windows: .channels, .data
%       (single, n_events x n_bin_edges x n_channels). Plot-only load.
%   epochs   - cfg.EPOCHS plus event_time, rate_hz, trial_n, cond_mask.
%       trial_n is only for paired tests across epochs.

behavior_to_align = {'x_hnd', 'y_hnd', 'x_eye', 'y_eye'};
behavior_to_downsample = {'x_eye', 'y_eye'};
n_beh_ch = numel(behavior_to_align);
conditions = cfg.CONDITIONS;

neuron_id = unit.Neuron_ID;
spike_times = unit.spike_times;
samples = unit.sample;
trials = unit.trial(:);
if isempty(trials)
    trial_timestamps = zeros(1, 0);
else
    trial_timestamps = [trials.timestamp];
end

window_rasters = cfg.WINDOWS(:);
n_windows = numel(window_rasters);
behavior = struct('channels', cell(n_windows, 1), 'data', cell(n_windows, 1));
pad_bins = psth_smooth_pad_bins(cfg);

for wi = 1:n_windows
    w = window_rasters(wi);
    bin_edges = w.t_start_s:cfg.raster_bin_size_s:w.t_end_s;
    bin_size_s = diff(bin_edges(1:2));
    left_edges = bin_edges(1) - (pad_bins:-1:1) * bin_size_s;
    right_edges = bin_edges(end) + (1:pad_bins) * bin_size_s;
    spike_edges = [left_edges, bin_edges, right_edges];
    event_times = unit.state.timestamps(unit.state.Value == w.align_state);
    event_times = event_times(isfinite(event_times));
    event_times = event_times(:);
    n_events = numel(event_times);
    n_bins_pad = numel(spike_edges) - 1;
    n_bins_beh = numel(bin_edges);

    raster_counts = zeros(n_events, n_bins_pad);
    trial_n = zeros(n_events, 1);
    beh_data = zeros(n_events, n_bins_beh, n_beh_ch, 'single');

    for k = 1:n_events
        rel_spike_times = spike_times - event_times(k);
        raster_counts(k, :) = histcounts(rel_spike_times, spike_edges);
        if isempty(trial_timestamps)
            trial_n(k) = 0;
        else
            idx = find(trial_timestamps <= event_times(k), 1, 'last');
            if isempty(idx), idx = 1; end
            trial_n(k) = idx;
        end

        % One sample of pad so interp1 has anchors on each side of the window.
        t_edges = event_times(k) + [bin_edges(1) bin_edges(end)];
        sample_idx_s = find(samples.timestamps > t_edges(1), 1, 'first') - 1;
        sample_idx_e = find(samples.timestamps < t_edges(2), 1, 'last') + 1;
        sample_idx = sample_idx_s:sample_idx_e;
        ts = samples.timestamps(sample_idx) - event_times(k);
        for bv = 1:n_beh_ch
            FN = behavior_to_align{bv};
            beh_data(k, :, bv) = resample_trace( ...
                ts, samples.(FN)(sample_idx), bin_edges, ...
                ismember(FN, behavior_to_downsample), cfg.raster_smooth_samples);
        end
    end

    window_rasters(wi).bin_edges = bin_edges(:);
    window_rasters(wi).bin_centers = spike_edges(1:end-1).' + bin_size_s/2;
    window_rasters(wi).bin_size_s = bin_size_s;
    window_rasters(wi).psth_pad_bins = pad_bins;
    window_rasters(wi).event_times = event_times;
    window_rasters(wi).trial_n = trial_n;
    window_rasters(wi).cond_mask = event_condition_mask(trials, trial_n, conditions);
    window_rasters(wi).raster = raster_counts > 0;
    behavior(wi).channels = behavior_to_align;
    behavior(wi).data = beh_data;
end

% One rate per align event. trial_n pairs epochs; cond_mask selects conditions.
epochs = cfg.EPOCHS(:);
for ei = 1:numel(epochs)
    e = epochs(ei);
    event_times = unit.state.timestamps(unit.state.Value == e.align_state);
    event_times = event_times(isfinite(event_times));
    event_times = event_times(:);
    n_events = numel(event_times);
    duration_s = e.t_end_s - e.t_start_s;

    rate_hz = nan(n_events, 1);
    trial_n = zeros(n_events, 1);

    for k = 1:n_events
        win_start = event_times(k) + e.t_start_s;
        win_end = event_times(k) + e.t_end_s;
        rate_hz(k) = sum(spike_times >= win_start & spike_times < win_end) / duration_s;
        if isempty(trial_timestamps)
            trial_n(k) = 0;
        else
            idx = find(trial_timestamps <= event_times(k), 1, 'last');
            if isempty(idx), idx = 1; end
            trial_n(k) = idx;
        end
    end

    epochs(ei).event_time = event_times;
    epochs(ei).rate_hz = rate_hz;
    epochs(ei).trial_n = trial_n;
    epochs(ei).cond_mask = event_condition_mask(trials, trial_n, conditions);
end

Neuron_ID = neuron_id; %#ok<NASGU>
windows = window_rasters; %#ok<NASGU>
raster_file = fullfile(cfg.roots.raster, sprintf('%s_raster.mat', neuron_id));
save(raster_file, 'Neuron_ID', 'windows', 'behavior', 'epochs', '-v7');
end

function cond_mask = event_condition_mask(trials, trial_n, conditions)
% n_events x n_conditions. An event can match several conditions.
n_events = numel(trial_n);
n_cond = numel(conditions);
cond_mask = false(n_events, n_cond);
valid = trial_n(:) >= 1 & trial_n(:) <= numel(trials);
if ~any(valid) || isempty(trials) || n_cond == 0
    return
end
ev_trials = trials(trial_n(valid));
for ci = 1:n_cond
    keep = true(nnz(valid), 1);
    param_names = fieldnames(conditions(ci).parameters);
    for pi = 1:numel(param_names)
        pn = param_names{pi};
        target_vals = conditions(ci).parameters.(pn);
        if ~isfield(ev_trials, pn)
            keep(:) = false;
            break
        end
        trial_vals = [ev_trials.(pn)].';
        keep = keep & ismember(trial_vals, target_vals(:));
    end
    cond_mask(valid, ci) = keep;
end
end

function pad_bins = psth_smooth_pad_bins(cfg)
% Bins of real spikes on each side. Same half-width as psth_conv_kernel
% in epp_plot_psth, so conv(...,'same') stays inside the requested window.
bin_size_s = cfg.raster_bin_size_s;
switch lower(char(cfg.psth.smoothing_kernel))
    case 'gaussian'
        sigma_bins = max(cfg.psth.smoothing_width_s / bin_size_s, eps);
        kernel_bins = 2 * max(1, ceil(4 * sigma_bins)) + 1;
    case 'box'
        kernel_bins = max(1, round(cfg.psth.smoothing_width_s / bin_size_s));
    otherwise
        kernel_bins = 1;
end
pad_bins = ceil((kernel_bins - 1) / 2);
end

function vq = resample_trace(ts, v, grid, downsample, smooth_samples)
% Drop held samples, interp1 onto grid, then boxcar-smooth.
ts = ts(:); v = v(:);
if downsample
    keep = unique([find([true; diff(v) ~= 0]); numel(v)]);
    ts = ts(keep); v = v(keep);
end
vq = single(interp1(ts, v, grid(:), 'linear'));
if smooth_samples > 1
    vq = movmean(vq, smooth_samples, 'omitnan');
end
end
