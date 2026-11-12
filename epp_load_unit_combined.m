function unit_combined = epp_load_unit_combined(unit_row, processed_blocks, cfg)
%EPP_LOAD_UNIT_COMBINED  Align one unit's spikes with its behavioral blocks.
%
%   unit_combined = epp_load_unit_combined(unit_row, processed_blocks, cfg)
%
%   Called by epp_build_rasters. epp_compute_statistics reads the epoch
%   rates that stage saves and does not call this again.
%
%   Spike file (UltraSort): dataspikes_rb<FFF>_ch<CCC>_<cfg.spike_threshold>.mat with
%   cluster_class (sort code, time in ms) and par.sr. segmentends come from
%   datafused_rb<FFF>_ch<CCC>.mat (par.segmentends), same as UltraLoadData.
%   Block i of the unit is time-shifted by segment start i. Trials then pass
%   through epp_enrich_trials_for_unit. Edge trials below cfg.min_edge_rate_hz
%   are removed (0 disables that cut).
%
%   unit_combined.trial, .state, .sample, .spike_times

session = unit_row.Session;
blocks = unit_row.Blocks;
perturbations = unit_row.Perturbation;
monkey = unit_row.monkey;
channel = unit_row.Channel;
filenumber = unit_row.Filenumber;
sortcode = unit_row.Unit;
threshold = cfg.spike_threshold;

day_token = sprintf('%08d', session);
session_folder = fullfile(cfg.roots.ephys_tanks, [monkey '_phys'], day_token);
spike_file = fullfile(session_folder, sprintf('dataspikes_rb%03d_ch%03d_%s.mat', filenumber, channel, threshold));
spike_data = load(spike_file, 'cluster_class', 'par');
cc = spike_data.cluster_class;
par = spike_data.par;

fused_file = fullfile(session_folder, sprintf('datafused_rb%03d_ch%03d.mat', filenumber, channel));
fused = load(fused_file, 'par');
segment_ends = fused.par.segmentends(:) ./ par.sr;
segment_starts = [0; segment_ends(1:end-1)];
spike_times = cc(cc(:, 1) == sortcode, 2) / 1000;

unit_combined = struct( ...
    'trial', struct([]), ...
    'state', struct('Value', [], 'timestamps', []), ...
    'sample', struct('x_hnd', [], 'y_hnd', [], 'x_eye', [], 'y_eye', [], 'timestamps', []), ...
    'spike_times', spike_times);

for bi = 1:numel(blocks)
    b = blocks(bi);
    idx_block = find([processed_blocks.session] == session & [processed_blocks.block] == b & ...
        strcmp({processed_blocks.monkey}, monkey), 1, 'first');
    if isempty(idx_block), continue; end
    bp = processed_blocks(idx_block);

    segment_start_s = segment_starts(bi);
    trial_aligned = bp.trial;
    [trial_aligned.block] = deal(b);
    [trial_aligned.perturbation] = deal(perturbations(bi));
    ts_cell = num2cell([trial_aligned.timestamp] + segment_start_s);
    [trial_aligned.timestamp] = ts_cell{:};
    state_aligned = bp.state;
    state_aligned.timestamps = state_aligned.timestamps + segment_start_s;
    sample_aligned = bp.sample;
    sample_aligned.timestamps = sample_aligned.timestamps + segment_start_s;

    aligned_chunk = struct('trial', trial_aligned, 'state', state_aligned, 'sample', sample_aligned);
    unit_combined = epp_concatenate_aligned_data(unit_combined, aligned_chunk);
end

unit_combined.trial = epp_enrich_trials_for_unit(unit_combined.trial, unit_row);
unit_combined = trim_low_rate_edge_trials(unit_combined, cfg);
end

function unit_combined = trim_low_rate_edge_trials(unit_combined, cfg)
% Drop leading and trailing trials below min_rate_hz (spikes/s).
% A trial spans [timestamp, next trial timestamp). The last trial ends at
% the last finite state or sample time. Trials in the middle are kept.
% min_rate_hz <= 0 does nothing. State and sample rows outside the kept
% span are removed. Spike times stay, so a kept trial can still see spikes
% in its pre-event window.

min_rate_hz = cfg.min_edge_rate_hz;
if isempty(min_rate_hz) || ~isscalar(min_rate_hz) || ~isfinite(min_rate_hz) || min_rate_hz <= 0
    return
end
if isempty(unit_combined.trial)
    return
end

[trial_starts, trial_ends] = unit_trial_bounds(unit_combined);
rates = trial_span_rates(unit_combined.spike_times, trial_starts, trial_ends);
below = ~isfinite(rates) | rates < min_rate_hz;
first_keep = find(~below, 1, 'first');
if isempty(first_keep)
    unit_combined.trial = unit_combined.trial([]);
    unit_combined = crop_state_and_sample(unit_combined, inf, -inf);
    return
end
last_keep = find(~below, 1, 'last');
if first_keep == 1 && last_keep == numel(below)
    return
end

if first_keep == 1
    t0 = -inf;
else
    t0 = trial_starts(first_keep);
end
if last_keep == numel(below)
    t1 = inf;
else
    t1 = trial_starts(last_keep + 1);
end
unit_combined.trial = unit_combined.trial(first_keep:last_keep);
unit_combined = crop_state_and_sample(unit_combined, t0, t1);
end

function [trial_starts, trial_ends] = unit_trial_bounds(unit_combined)
trial_starts = [unit_combined.trial.timestamp].';
n_trials = numel(trial_starts);
data_end = last_data_time(unit_combined);
if n_trials > 1
    trial_ends = [trial_starts(2:end); data_end];
else
    trial_ends = data_end;
end
end

function t_end = last_data_time(unit_combined)
candidates = [unit_combined.state.timestamps(:); unit_combined.sample.timestamps(:)];
candidates = candidates(isfinite(candidates));
if isempty(candidates)
    t_end = NaN;
else
    t_end = max(candidates);
end
end

function rates = trial_span_rates(spike_times, trial_starts, trial_ends)
n_trials = numel(trial_starts);
duration = trial_ends(:) - trial_starts(:);
rates = NaN(n_trials, 1);
edges = [trial_starts(:); trial_ends(end)];
monotonic = numel(edges) == n_trials + 1 && all(isfinite(edges)) && all(diff(edges) > 0);
spikes = spike_times(isfinite(spike_times));
if monotonic
    rates = histcounts(spikes, edges).' ./ duration;
    return
end
for i = 1:n_trials
    if ~isfinite(duration(i)) || duration(i) <= 0
        continue
    end
    rates(i) = sum(spikes >= trial_starts(i) & spikes < trial_ends(i)) / duration(i);
end
end

function unit_combined = crop_state_and_sample(unit_combined, t0, t1)
keep = unit_combined.state.timestamps >= t0 & unit_combined.state.timestamps < t1;
unit_combined.state.Value = unit_combined.state.Value(keep);
unit_combined.state.timestamps = unit_combined.state.timestamps(keep);

sample_keep = unit_combined.sample.timestamps >= t0 & unit_combined.sample.timestamps < t1;
sample_fields = {'x_hnd', 'y_hnd', 'x_eye', 'y_eye', 'timestamps'};
for i = 1:numel(sample_fields)
    fn = sample_fields{i};
    if numel(unit_combined.sample.(fn)) == numel(sample_keep)
        unit_combined.sample.(fn) = unit_combined.sample.(fn)(sample_keep);
    end
end
end

function combined = epp_concatenate_aligned_data(combined, chunk)
combined.trial = [combined.trial; chunk.trial(:)];                                          % #ok<AGROW>
combined.state.Value = [combined.state.Value; chunk.state.Value(:)];                        % #ok<AGROW>
combined.state.timestamps = [combined.state.timestamps; chunk.state.timestamps(:)];         % #ok<AGROW>
combined.sample.x_hnd = [combined.sample.x_hnd; chunk.sample.x_hnd(:)];                     % #ok<AGROW>
combined.sample.y_hnd = [combined.sample.y_hnd; chunk.sample.y_hnd(:)];                     % #ok<AGROW>
combined.sample.x_eye = [combined.sample.x_eye; chunk.sample.x_eye(:)];                     % #ok<AGROW>
combined.sample.y_eye = [combined.sample.y_eye; chunk.sample.y_eye(:)];                     % #ok<AGROW>
combined.sample.timestamps = [combined.sample.timestamps; chunk.sample.timestamps(:)];      % #ok<AGROW>
end
