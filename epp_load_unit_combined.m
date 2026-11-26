function unit = epp_load_unit_combined(UI, Blocks_in, cfg)
%EPP_LOAD_UNIT_COMBINED  Spikes + concatenated trials/states for one unit.
%
%   unit = epp_load_unit_combined(unit_row, accumulated_blocks, cfg)
%
%   Inputs:
%     UI        - one unit_info row (Monkey, Session, Filenumber, Channel,
%                 Unit, Blocks, Perturbation, Hemisphere, Neuron_ID, ...).
%     Blocks_in - block_payload array; filtered to this Monkey/Session.
%     cfg       - spike_threshold, min_edge_rate_hz, roots.ephys_tanks.
%
%   Spikes: dataspikes_rb<FFF>_ch<CCC>_<threshold>.mat (ms -> s).
%   Segments: datafused_rb<FFF>_ch<CCC>.mat par.segmentends -> segment_starts.
%   Trials/states shifted by segment_starts(bi) into concatenated time;
%   trial.block / perturbation stamped; hemisphere flip on trial fields.
%   Edge trials with FR <= cfg.min_edge_rate_hz trimmed (states cropped too).
%   Samples are NOT kept on unit (plot uses accumulated_blocks).
%
%   unit: UI fields + segment_starts, spike_times, trial, state.

Blocks_in =Blocks_in([Blocks_in.session] == UI.Session & strcmp({Blocks_in.monkey}, UI.Monkey));

day_token = sprintf('%08d', UI.Session);
session_folder = fullfile(cfg.roots.ephys_tanks, [UI.Monkey '_phys'], day_token);
filename_spikes= sprintf('dataspikes_rb%03d_ch%03d_%s.mat', UI.Filenumber, UI.Channel, cfg.spike_threshold);
filename_fused = sprintf('datafused_rb%03d_ch%03d.mat', UI.Filenumber, UI.Channel);

spike_data = load(fullfile(session_folder, filename_spikes), 'cluster_class', 'par');
cc = spike_data.cluster_class;
par = spike_data.par;

fused = load(fullfile(session_folder, filename_fused), 'par');
segment_ends = fused.par.segmentends(:) ./ par.sr;

unit=UI;
unit.segment_starts = [0; segment_ends(1:end-1)];
unit.spike_times = cc(cc(:, 1) == UI.Unit, 2) / 1000;
unit.trial=struct([]);
unit.state=struct('Value', [], 'timestamps', []);

%% concatenate_across_blocks
unit=concatenate_blocks(unit,Blocks_in);

%% add hemisphere and adjust space/hand to contra/ipsi
unit = hemisphere_transformation(unit);

%% trim trials (beginning and end across ALL blocks !)
unit = trim_low_rate_edge_trials(unit, cfg);
end

function unit = concatenate_blocks(unit,Blocks_in)

for bi = 1:numel(unit.Blocks)
    b = unit.Blocks(bi);
    idx_block = [Blocks_in.block] == b;
    %if isempty(idx_block), continue; end
    bp = Blocks_in(idx_block);

    segment_start_s = unit.segment_starts(bi);
    trial_aligned = bp.trial;
    [trial_aligned.block] = deal(b);
    [trial_aligned.perturbation] = deal(unit.Perturbation(bi));
    ts_cell = num2cell([trial_aligned.timestamp] + segment_start_s);
    [trial_aligned.timestamp] = ts_cell{:};

    state_timestamps = bp.state.timestamps + segment_start_s;
    bp.sample.timestamps = bp.sample.timestamps + segment_start_s;

    unit.trial = [unit.trial; trial_aligned(:)];                                          % #ok<AGROW>

    unit.state.Value = [unit.state.Value; bp.state.Value(:)];                        % #ok<AGROW>
    unit.state.timestamps = [unit.state.timestamps; state_timestamps(:)];         % #ok<AGROW>
end

end

function unit = hemisphere_transformation(unit)
trial=unit.trial;
rec_hemi_sign=NaN;
if ~isempty(unit.Hemisphere)
    hemisphere = lower(unit.Hemisphere(1));
    if hemisphere == 'l'; rec_hemi_sign = 1; elseif hemisphere == 'r';  rec_hemi_sign = -1; end
end

if isnan(rec_hemi_sign)
    error('epp_enrich_trials_for_unit:InvalidHemisphere', ...
        'Unit %s has missing or invalid hemisphere. Expected left (L/l) or right (R/r).', unit.Neuron_ID);
end

target_pos_fields = {'tar_pos', 'nct_pos', 'fix_pos', 'reach_tar_pos_closest', 'saccade_tar_pos_closest', 'hemifield'};
hand_fields = {'demanded_hand', 'used_hand', 'reach_hand'};

for f = target_pos_fields
    if ~isfield(trial, f{1}), continue; end
    IN = [trial.(f{:})];
    if isreal(IN)
        OUT = num2cell(real(IN) * rec_hemi_sign);
    else
        OUT = num2cell(real(IN) * rec_hemi_sign + 1i * imag(IN));
    end
    [trial.(f{:})] = deal(OUT{:});
end

for f = hand_fields
    if ~isfield(trial, f{1}), continue; end
    IN = [trial.(f{:})];
    OUT = IN;
    OUT(IN == 1) = -1; % left hand
    OUT(IN == 2) = 1;  % right hand
    OUT = num2cell(OUT * rec_hemi_sign);
    [trial.(f{:})] = deal(OUT{:});
end
unit.trial=trial;
end

function unit = trim_low_rate_edge_trials(unit, cfg)
% Drop leading and trailing trials at or below min_rate_hz (spikes/s).
% A trial spans [timestamp, next trial timestamp). The last trial ends at
% the last finite state or sample time. Trials in the middle are kept.
% State and sample rows outside the kept span are removed. Spike times
% stay, so a kept trial can still see spikes in its pre-event window.

min_rate_hz = cfg.min_edge_rate_hz;
if isempty(min_rate_hz) || ~isscalar(min_rate_hz) || ~isfinite(min_rate_hz) || isempty(unit.trial)
    return
end

trial_starts = [unit.trial.timestamp].';
data_end = max(unit.state.timestamps);

if numel(trial_starts) > 1
    trial_ends = [trial_starts(2:end); data_end];
else
    trial_ends = data_end;
end

duration = trial_ends(:) - trial_starts(:);
edges = [trial_starts(:); trial_ends(end)];
rates = histcounts(unit.spike_times, edges).' ./ duration;

below = rates <= min_rate_hz;
first_keep = find(~below, 1, 'first');
if isempty(first_keep)
    return
end
last_keep = find(~below, 1, 'last');

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

%% crop now

% trials
unit.trial = unit.trial(first_keep:last_keep);
%states
keep = unit.state.timestamps >= t0 & unit.state.timestamps < t1;
unit.state.Value = unit.state.Value(keep);
unit.state.timestamps = unit.state.timestamps(keep);
end
