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
%       Neuron_ID, windows, epochs
%
%   windows - cfg.WINDOWS plus, per window:
%       bin_edges, bin_centers, bin_size_s, psth_pad_bins,
%       event_times (block-local ephys time), trial_block, trial_n
%       (within-block trial number), cond_mask, raster (logical).
%       Monkey/Session come from the unit row at plot time.
%       No behavior traces — plot resamples from accumulated_blocks.
%   epochs  - cfg.EPOCHS plus event_time (block-local), rate_hz,
%       trial_n (unit-local index into unit.trial), cond_mask.
%       Saved under cfg.roots.epochs. No behavior backtrack.

conditions = cfg.CONDITIONS;

neuron_id = unit.Neuron_ID;
spike_times = unit.spike_times;
trials = unit.trial(:);
if isempty(trials)
    trial_timestamps = zeros(1, 0);
else
    trial_timestamps = [trials.timestamp];
end

window_rasters = cfg.WINDOWS(:);
pad_bins = psth_smooth_pad_bins(cfg);

for wi = 1:numel(window_rasters)
    w = window_rasters(wi);
    bin_edges = w.t_start_s:cfg.raster_bin_size_s:w.t_end_s;
    bin_size_s = diff(bin_edges(1:2));
    left_edges = bin_edges(1) - (pad_bins:-1:1) * bin_size_s;
    right_edges = bin_edges(end) + (1:pad_bins) * bin_size_s;
    spike_edges = [left_edges, bin_edges, right_edges];
    event_times = unit.state.timestamps(unit.state.Value == w.align_state);
    %event_times = event_times(isfinite(event_times));
    %event_times = event_times(:);
    n_events = numel(event_times);
    n_bins_pad = numel(spike_edges) - 1;

    raster_counts = zeros(n_events, n_bins_pad);
    trial_n = zeros(n_events, 1);
    trial_block = zeros(n_events, 1);
    unit_trial_idx = zeros(n_events, 1);
    event_times_local = nan(n_events, 1);

    for k = 1:n_events
        rel_spike_times = spike_times - event_times(k);
        raster_counts(k, :) = histcounts(rel_spike_times, spike_edges);
        if isempty(trial_timestamps)
            continue
        end
        idx = find(trial_timestamps <= event_times(k), 1, 'last');
        if isempty(idx), idx = 1; end
        tr = trials(idx);
        unit_trial_idx(k) = idx;
        trial_n(k) = tr.n;
        trial_block(k) = tr.block;
        event_times_local(k) = event_time_block_local(event_times(k), tr.block, unit);
    end

    window_rasters(wi).bin_edges = bin_edges(:);
    window_rasters(wi).bin_centers = spike_edges(1:end-1).' + bin_size_s/2;
    window_rasters(wi).bin_size_s = bin_size_s;
    window_rasters(wi).psth_pad_bins = pad_bins;
    window_rasters(wi).event_times = event_times_local;
    window_rasters(wi).trial_n = trial_n;
    window_rasters(wi).trial_block = trial_block;
    window_rasters(wi).cond_mask = event_condition_mask(trials, unit_trial_idx, conditions);
    window_rasters(wi).raster = raster_counts > 0;
end

% One rate per align event. trial_n = unit-local index; cond_mask selects.
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
    event_times_local = nan(n_events, 1);

    for k = 1:n_events
        win_start = event_times(k) + e.t_start_s;
        win_end = event_times(k) + e.t_end_s;
        rate_hz(k) = sum(spike_times >= win_start & spike_times < win_end) / duration_s;
        if isempty(trial_timestamps)
            continue
        end
        idx = find(trial_timestamps <= event_times(k), 1, 'last');
        if isempty(idx), idx = 1; end
        tr = trials(idx);
        trial_n(k) = idx;
        event_times_local(k) = event_time_block_local(event_times(k), tr.block, unit);
    end

    epochs(ei).event_time = event_times_local;
    epochs(ei).rate_hz = rate_hz;
    epochs(ei).trial_n = trial_n;
    epochs(ei).cond_mask = event_condition_mask(trials, trial_n, conditions);
end

Neuron_ID = neuron_id; %#ok<NASGU>
windows = window_rasters; %#ok<NASGU>
raster_file = fullfile(cfg.roots.raster, sprintf('%s_raster.mat', neuron_id));
epoch_file = fullfile(cfg.roots.epochs, sprintf('%s_epochs.mat', neuron_id));
save(raster_file, 'Neuron_ID', 'windows');
save(epoch_file, 'Neuron_ID', 'epochs');
end

function t_local = event_time_block_local(event_time_concat, block, unit)
% Undo the segment_start shift for this block.
bi = find(unit.Blocks == block, 1, 'first');
if isempty(bi)
    t_local = NaN;
    return
end
t_local = event_time_concat - unit.segment_starts(bi);
end

function cond_mask = event_condition_mask(trials, unit_trial_idx, conditions)
% n_events x n_conditions. unit_trial_idx indexes into trials (not trial.n).
n_events = numel(unit_trial_idx);
n_cond = numel(conditions);
cond_mask = false(n_events, n_cond);
valid = unit_trial_idx(:) >= 1 & unit_trial_idx(:) <= numel(trials);
if ~any(valid) || isempty(trials) || n_cond == 0
    return
end
ev_trials = trials(unit_trial_idx(valid));
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
