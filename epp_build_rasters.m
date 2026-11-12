function epp_build_rasters(cfg, unit_row, processed_blocks)
%EPP_BUILD_RASTERS  Stage 2: one unit's event-aligned raster.
%
%   epp_build_rasters(cfg, unit_row, processed_blocks)
%
%   Inputs:
%     cfg - struct from epp_load_cfg. Uses cfg.WINDOWS, cfg.EPOCHS,
%           cfg.raster_bin_size_s, cfg.roots.raster. Spike alignment is
%           epp_load_unit_combined.
%     unit_row - one element of unit_info.
%     processed_blocks - stage-1 payloads from epp_load_saved_processed_blocks.
%           epp_initiation loads these once and passes them in.
%
%   Outputs:
%     none. Per unit, writes cfg.roots.raster/<Neuron_ID>_raster.mat
%     (variable raster_data).
%
%   raster_data.Neuron_ID - char.
%   raster_data.windows   - cfg.WINDOWS plus, per window:
%       bin_edges, bin_centers, bin_size_s, event_times, raster (counts),
%       rate_hz, trial (behavior struct at each align event), behavior
%       (x_hnd, y_hnd, x_eye, y_eye on bin_edges).
%   raster_data.epochs    - cfg.EPOCHS plus, per epoch:
%       event_time, rate_hz, trial, trial_n. rate_hz is spikes/s for that
%       align event alone. trial_n is the trial the event fell in (1 when
%       the event is before every trial). epp_compute_statistics matches
%       paired tests on trial_n.

behavior_to_align = {'x_hnd','y_hnd','x_eye','y_eye'};
behavior_to_downsample = {'x_eye','y_eye'};

unit = epp_load_unit_combined(unit_row, processed_blocks, cfg);
neuron_id = unit_row.Neuron_ID;

spike_times = unit.spike_times;
samples = unit.sample;
window_rasters = cfg.WINDOWS(:);

for wi = 1:numel(window_rasters)
    clear behavior w_trials raster_counts
    w = window_rasters(wi);
    bin_edges = w.t_start_s:cfg.raster_bin_size_s:w.t_end_s;
    bin_size_s = diff(bin_edges(1:2));
    event_times = unit.state.timestamps(unit.state.Value == w.align_state);
    event_times = event_times(isfinite(event_times));
    event_times = event_times(:);
    n_events = numel(event_times);

    raster_counts = zeros(n_events, numel(bin_edges)-1);
    for k = 1:n_events
        rel_spike_times = spike_times - event_times(k);
        raster_counts(k,:) = histcounts(rel_spike_times, bin_edges);
        idx = find([unit.trial.timestamp] <= event_times(k), 1, 'last');
        if isempty(idx), idx = 1; end
        w_trials(k) = unit.trial(idx); %#ok<AGROW>

        % One bin of pad so interp1 has an anchor on each side of the window.
        t_edges = event_times(k) + [bin_edges(1) bin_edges(end)];
        sample_idx_s = find(samples.timestamps > t_edges(1),1,'first')-1;
        sample_idx_e = find(samples.timestamps < t_edges(2),1,'last')+1;
        sample_idx = sample_idx_s : sample_idx_e;
        ts = samples.timestamps(sample_idx) - event_times(k);
        for bv = 1:numel(behavior_to_align)
            FN = behavior_to_align{bv};
            behavior(k).(FN) = resample_trace(ts, samples.(FN)(sample_idx), bin_edges, ismember(FN, behavior_to_downsample), cfg.raster_smooth_samples); %#ok<AGROW>
        end
    end

    window_rasters(wi).bin_edges = bin_edges(:);
    window_rasters(wi).bin_centers = bin_edges(1:end-1).' + bin_size_s/2;
    window_rasters(wi).bin_size_s = bin_size_s;
    window_rasters(wi).event_times = event_times;
    window_rasters(wi).raster = raster_counts;
    window_rasters(wi).rate_hz = raster_counts ./ bin_size_s;
    window_rasters(wi).trial = w_trials;
    window_rasters(wi).behavior = behavior;
end

% One rate per align event. trial_n is the trial index used for paired tests.
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
    trial_n = nan(n_events, 1);

    for k = 1:n_events
        win_start = event_times(k) + e.t_start_s;
        win_end = event_times(k) + e.t_end_s;
        rate_hz(k) = sum(spike_times >= win_start & spike_times < win_end) / duration_s;
        idx = find([unit.trial.timestamp] <= event_times(k), 1, 'last');
        if isempty(idx), idx = 1; end
        epoch_trials(k) = unit.trial(idx); %#ok<AGROW>
        trial_n(k) = idx;
    end
    
    epochs(ei).event_time = event_times;
    epochs(ei).rate_hz = rate_hz;
    epochs(ei).trial = epoch_trials;
    epochs(ei).trial_n = trial_n;
end

raster_data = struct();
raster_data.Neuron_ID = neuron_id;
raster_data.windows = window_rasters;
raster_data.epochs = epochs;
raster_file = sprintf('%s_raster.mat', neuron_id);
save(fullfile(cfg.roots.raster, raster_file), 'raster_data');
end

function vq = resample_trace(ts, v, grid, downsample, smooth_samples)
% Drop held samples, interp1 onto grid, then boxcar-smooth.
ts = ts(:); v = v(:);
if downsample
    keep = unique([find([true; diff(v)~=0]); numel(v)]);
    ts = ts(keep); v = v(keep);
end
vq = interp1(ts, v, grid(:), 'linear');
if smooth_samples > 1
    vq = movmean(vq, smooth_samples, 'omitnan');
end
end


