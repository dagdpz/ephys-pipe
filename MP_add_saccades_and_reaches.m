function trial = MP_add_saccades_and_reaches(trial, cfg)
% MP_add_saccades_and_reaches
% Standalone extraction of movement onsets/offsets from a trial struct.
%
% This function takes a MonkeyPsych trial struct array, applies the core
% correction/detection steps used in monkeypsych_analyze_working, and
% returns a structure with added per-trial fields:
%   - saccade_ini
%   - saccade_end
%   - reach_onset
%   - reach_offset
%
% Usage:
%   trial = MP_add_saccades_and_reaches(trial)
%   trial = MP_add_saccades_and_reaches(trial, cfg)
%
% Key cfg options (string modes):
%   cfg.saccade_definition:
%       'closest_big_enough'    (old 1, default)
%       'biggest_close_enough'  (old 2)
%       'last_in_state'         (old 3)
%       'first_in_state'        (old 4)
%       'first_big_enough'      (old 5)
%       'first_end_in_target'   (old 10)
%       'first_start_in_target' (old 11)
%
%   cfg.reach_definition:
%       'first_touch'
%       'first_touch_in_target' (default)
%       'state_change_position'

if nargin < 2; cfg = struct(); end
cfg = apply_default_cfg(cfg);
if cfg.enable_cleaning; trial = clean_trials(trial); end

% Constants/parameters aligned to monkeypsych defaults.
MS = get_ma_states();

trial = initialize_trial_nan_fields(trial, cfg);
for n = 1:numel(trial)
    tr = trial(n);
    
    effector_sr = get_expected_effectors(tr);
    [observed_state, state_1ao, can_infer_target] = get_state_definition(tr,MS);
    [trial(n), all_eye_target_pos, all_hnd_target_pos] = fill_target_fields(tr, effector_sr, observed_state, can_infer_target);
    
    % Reach-hand bookkeeping/debug (monkeypsych-style semantics):
    % demanded = trial.task.reach_hand, used = trial.reach_hand
    trial(n).demanded_hand = NaN;
    trial(n).used_hand = NaN;
    if isfield(tr, 'task') && isfield(tr.task, 'reach_hand') && ~isempty(tr.task.reach_hand)
        trial(n).demanded_hand = tr.task.reach_hand;
    end
    if isfield(tr, 'reach_hand')
        if ~isempty(tr.reach_hand)
            trial(n).used_hand = tr.reach_hand;
        end
    end
        
    % --- Correction step: optional eye offset correction at fixation hold ---
    if cfg.correct_offset && isfield(tr, 'eye') && isfield(tr.eye, 'fix') && isfield(tr.eye.fix, 'pos') && numel(tr.eye.fix.pos) >= 2 && max(tr.state) > 2
        [tr.x_eye, tr.y_eye] = offset_corrected(tr.x_eye, tr.y_eye, tr.state, MS.FIX_HOL, tr.eye.fix.pos(1), tr.eye.fix.pos(2));
    end
    
    % Trial time axis as in monkeypsych (relative to FIX_ACQ onset).
    idx_fix_acq = find(tr.state == MS.FIX_ACQ, 1, 'first');
    
    t_state2 = tr.tSample_from_time_start(idx_fix_acq);
    tr.time_axis = tr.tSample_from_time_start - t_state2;
    
    if numel(tr.time_axis) < 3
        continue;
    end
    
    % Observed state from task-type definition (monkeypsych-style).
    smp_state_target = find(tr.state == observed_state);
    obs_state_onset = NaN;
    if ~isempty(smp_state_target)
        obs_state_onset = tr.time_axis(smp_state_target(1));
    end
    
    %% Reach onset/offset (same logic core as monkeypsych)
        
    % Type-1 fallback from sensor release in FIX_ACQ (same correction path).
    if effector_sr(2) && isfield(tr, 'type') && tr.type == 1 && all(isfield(tr, {'sen_L','sen_R'}))
        idx_release_fixacq = find(logical([diff(tr.sen_L) ~= 0 | diff(tr.sen_R) ~= 0; 0]) & tr.state == MS.FIX_ACQ, 1, 'first');
        idx_fix = find(tr.state == MS.FIX_HOL, 1, 'first');
        if ~isempty(idx_release_fixacq)
            trial(n).reach_onset = tr.time_axis(idx_release_fixacq);
            if ~isempty(idx_fix)
                trial(n).reach_offset = tr.time_axis(idx_fix);
                trial(n).reach_dur = trial(n).reach_offset - trial(n).reach_onset;
            end
            trial(n).reach_lat = trial(n).reach_onset;
        end
        
    elseif effector_sr(2)
        smp_reaching_nan = find(isnan(tr.x_hnd) & tr.state == observed_state)';
        if ~isempty(smp_state_target) && ~isempty(smp_reaching_nan) && smp_reaching_nan(end) ~= numel(tr.state)
            smp_reaching = [smp_reaching_nan(1)-1, smp_reaching_nan, smp_reaching_nan(end)+1];
            times_reach = tr.time_axis(smp_reaching);
            reach_ini = times_reach(2);
            
            smp_state_1ao = find(tr.state == state_1ao);
            if isempty(smp_state_1ao)
                smp_state_1ao = numel(tr.state);
            end
            
            
            switch cfg.reach_definition
                case 'first_touch'
                    reach_endpos = tr.x_hnd(smp_reaching(end)) + 1i*tr.y_hnd(smp_reaching(end));
                case 'first_touch_in_target'
                    reach_endpos = tr.x_hnd(smp_state_target(end)) + 1i*tr.y_hnd(smp_state_target(end));
                case 'state_change_position'
                    reach_endpos = tr.x_hnd(smp_state_1ao(end)) + 1i*tr.y_hnd(smp_state_1ao(end));
            end
            
            reach_offset = tr.time_axis(smp_state_target(end));
            if isfinite(real(reach_endpos))
                trial(n).reach_endpos = reach_endpos;
                trial(n).reach_onset = reach_ini;
                trial(n).reach_offset = reach_offset;
                trial(n).reach_lat = reach_ini - obs_state_onset;
                trial(n).reach_dur = reach_offset - reach_ini;
                trial(n).reach_tar_pos_closest = get_closest_target_pos(reach_endpos, all_hnd_target_pos);
            end
        end
    end
    
    
    %% Saccade onset/offset (same thresholding core as monkeypsych)
    if effector_sr(1)
        
        
        interp_idx = 1:numel(tr.x_eye);
        interp_end = tr.time_axis(end);
        
        % Optional downsampling keeps only "real" eye-samples (non-repeated x/y pairs),
        % then interpolates from those anchors onto the uniform axis.
        if cfg.downsampling
            repeats_idx = ([NaN; diff(tr.x_eye(:))] == 0) & ([NaN; diff(tr.y_eye(:))] == 0);
            interp_idx = unique([find(~repeats_idx); numel(repeats_idx)])';
            if numel(interp_idx) < 2
                interp_idx = [1 numel(tr.x_eye)];
            end
            interp_end = tr.time_axis(interp_idx(end));
        end
        
        % Interpolation + smoothing.
        time_axis = ceil(tr.time_axis(1)*cfg.interpolated_SR)/cfg.interpolated_SR : 1/cfg.interpolated_SR : interp_end;
        x = tr.x_eye(interp_idx);
        y = tr.y_eye(interp_idx);
        t = tr.time_axis(interp_idx);
        if numel(x) < 2 || all(diff(x) == 0) || all(diff(y) == 0)
            x_eye = repmat(x(1), 1, numel(time_axis));
            y_eye = repmat(y(1), 1, numel(time_axis));
        else
            x_eye = interp1(t, x, time_axis, 'linear');
            y_eye = interp1(t, y, time_axis, 'linear');
        end
        x_eye = filter_et(x_eye, cfg.smoothing_samples);
        y_eye = filter_et(y_eye, cfg.smoothing_samples);
        eye_vel = [0 sqrt((diff(x_eye).*cfg.interpolated_SR).^2 + (diff(y_eye).*cfg.interpolated_SR).^2)];
        eye_vel = filter_et(eye_vel, cfg.smoothing_samples);
        
        % Extend state labels onto interpolated time base.
        idx_state_changes = find([true; diff(tr.state) ~= 0]);
        states_present = tr.state(idx_state_changes);
        times_state_changed = [tr.time_axis(idx_state_changes); time_axis(end)];
        state = zeros(size(time_axis));
        for s = 1:numel(states_present)
            state_start = floor(times_state_changed(s)*cfg.interpolated_SR)/cfg.interpolated_SR;
            state_end = floor(times_state_changed(s+1)*cfg.interpolated_SR)/cfg.interpolated_SR;
            state(time_axis >= state_start & time_axis < state_end) = states_present(s);
        end
        
        
        not_iti = ~(state == MS.INI_TRI | state == MS.ITI);
        sac_above = (eye_vel >= cfg.sac_ini_t) & not_iti;
        sac_under = (eye_vel <= cfg.sac_end_t) & not_iti;
        between_start = find([diff(sac_above) == -1 | diff(sac_under) == -1 false]);
        between_end = find([false diff(sac_above) == 1 | diff(sac_under) == 1]);
        
        if ~isempty(between_end) && ~isempty(between_start) && between_end(1) <= between_start(1)
            between_end(1) = [];
        end
        if isempty(between_start) || isempty(between_end)
            continue;
        end
        
        startfromlow = eye_vel(between_start) <= cfg.sac_end_t;
        startfromhigh = eye_vel(between_start) >= cfg.sac_ini_t;
        endinlow = eye_vel(between_end) <= cfg.sac_end_t;
        endinhigh = eye_vel(between_end) >= cfg.sac_ini_t;
        
        if numel(between_end) + 1 == numel(between_start)
            between_end(end+1) = numel(time_axis);
            endinlow(end+1) = true;
            endinhigh(end+1) = false;
        end
        
        sac_start = between_end(endinhigh & startfromlow);
        sac_end = between_end(endinlow & startfromhigh);
        if isempty(sac_start) && ~isempty(sac_end) && ~sac_above(1)
            sac_start = 1;
        elseif ~isempty(sac_end) && ~isempty(sac_start) && sac_end(1) < sac_start(1)
            sac_end(1) = [];
        end
        if numel(sac_end) + 1 == numel(sac_start)
            sac_end(end+1) = numel(time_axis);
        end
        if isempty(sac_start) || isempty(sac_end)
            continue;
        end
        
        durations = time_axis(sac_end) - time_axis(sac_start);
        sacidx_dur = durations .* 1000 >= cfg.sac_min_dur; % intentionally same unit handling as original code
        sac_start = sac_start(sacidx_dur);
        sac_end = sac_end(sacidx_dur);
        if isempty(sac_start)
            continue;
        end
        
        nsacc = min(cfg.nsacc_max, numel(sac_start));
        sac_start = sac_start(1:nsacc);
        sac_end = sac_end(1:nsacc);
        saccade_vel = nan(1, nsacc);
        for t = 1:nsacc
            saccade_vel(t) = max(eye_vel(sac_start(t):sac_end(t)));
        end
        
        states_of_saccades = state(sac_start);
        % Selection criterion: saccade initiation must be in observed state.
        % End point is taken from paired detection and can fall outside observed state.
        idx_candidates = (states_of_saccades == observed_state);
        if ~any(idx_candidates)
            continue;
        end
        ini_candidates = time_axis(sac_start(idx_candidates));
        end_candidates = time_axis(sac_end(idx_candidates));
        startpos_candidates = x_eye(sac_start(idx_candidates)) + 1i*y_eye(sac_start(idx_candidates));
        endpos_candidates = x_eye(sac_end(idx_candidates)) + 1i*y_eye(sac_end(idx_candidates));
        vel_candidates = saccade_vel(idx_candidates);
        amp_candidates = abs(endpos_candidates - startpos_candidates);
        
        sel_idx = select_saccade_index( ...
            cfg, amp_candidates, endpos_candidates, startpos_candidates, trial(n).tar_pos, ...
            trial(n).fix_pos);
        
        trial(n).saccade_ini = ini_candidates(sel_idx);
        trial(n).saccade_end = end_candidates(sel_idx);
        trial(n).saccade_startpos = startpos_candidates(sel_idx);
        trial(n).saccade_endpos = endpos_candidates(sel_idx);
        trial(n).saccade_vel = vel_candidates(sel_idx);
        trial(n).saccade_amplitude = amp_candidates(sel_idx);
        trial(n).saccade_lat = trial(n).saccade_ini - obs_state_onset;
        trial(n).saccade_dur = trial(n).saccade_end - trial(n).saccade_ini;
        trial(n).saccade_tar_pos_closest = get_closest_target_pos(trial(n).saccade_endpos, all_eye_target_pos);
    end

    % Insert movement event states into existing trial state timeline and
    % keep it ordered by state onset.
    trial(n) = insert_movement_states(trial(n), MS);
end

% Drop trial.task from output while keeping computations above intact.
if ~isempty(trial) && isfield(trial, 'task')
    trial = rmfield(trial, 'task');
end

end

function trial = initialize_trial_nan_fields(trial, cfg)
nan_templates = struct( ...
    'scalar', NaN, ...
    'complex_scalar', NaN + 1i*NaN, ...
    'rgb', [NaN NaN NaN]);

field_groups = struct( ...
    'scalar', {{ ...
    'saccade_ini','saccade_end','saccade_vel','saccade_amplitude','saccade_lat','saccade_dur', ...
    'reach_onset','reach_offset','reach_lat','reach_dur','demanded_hand','used_hand', ...
    'choice','target_selected_idx','target_nonselected_idx','hemifield'}}, ...
    'complex_scalar', {{ ...
    'saccade_startpos','saccade_endpos', ...
    'reach_endpos','saccade_tar_pos_closest','reach_tar_pos_closest','tar_pos','nct_pos','fix_pos'}}, ...
    'rgb', {{ ...
    'col_dim','col_bri','nct_col_dim','nct_col_bri'}});

group_names = fieldnames(field_groups);
for g = 1:numel(group_names)
    group_name = group_names{g};
    template = nan_templates.(group_name);
    fields = field_groups.(group_name);
    for f = 1:numel(fields)
        if ~isfield(trial, fields{f})
            [trial.(fields{f})] = deal(template);
        end
    end
end
end

function tr_out = insert_movement_states(tr_in, MS)
tr_out = tr_in;
state_values = tr_in.states;
state_onsets = tr_in.states_onset;

movement_state_values = [MS.SAC_INI MS.SAC_END MS.REA_INI MS.REA_END];
movement_onsets = [tr_in.saccade_ini tr_in.saccade_end tr_in.reach_onset tr_in.reach_offset];

state_values_merged = [state_values movement_state_values];
state_onsets_merged = [state_onsets movement_onsets + state_onsets(state_values==2)]; %movement_onsets were realtive to state 2 of this trial so far
[state_onsets_sorted, sort_idx] = sort(state_onsets_merged, 'ascend');
state_values_sorted = state_values_merged(sort_idx);

tr_out.states = state_values_sorted;
tr_out.states_onset = state_onsets_sorted;
end

%% Task/state helper functions
function effector_sr = get_expected_effectors(trial)
switch trial.effector
    case {0,3}
        effector_sr = [1 0];
    case {1,4,6}
        effector_sr = [0 1];
    case 2
        effector_sr = [1 1];
    otherwise
    effector_sr = [NaN NaN];
end
end

function [observed_state, state_1ao, can_infer_target] = get_state_definition(trial,MS)
% State sequence depends on trial.type.

switch trial.type
    case 1 % fixation only
        all_states = [0 MS.FIX_ACQ MS.FIX_HOL MS.SUCCESS_ABORT MS.REWARD MS.ITI];
        observed_state = MS.FIX_HOL;
        information_state = MS.FIX_ACQ;
    case 2 % direct movement
        all_states = [MS.FIX_ACQ MS.FIX_HOL MS.TAR_ACQ MS.TAR_HOL MS.SUCCESS_ABORT MS.REWARD MS.ITI];
        observed_state = MS.TAR_ACQ;
        information_state = MS.TAR_ACQ;
    case 2.5 % direct movement with dimmed targets
        all_states = [MS.FIX_ACQ MS.FIX_HOL MS.CUE_ON MS.MEM_PER MS.TAR_ACQ MS.TAR_HOL MS.SUCCESS_ABORT MS.REWARD MS.ITI];
        observed_state = MS.TAR_ACQ;
        information_state = MS.CUE_ON;
    case 3 % memory tasks
        all_states = [MS.FIX_ACQ MS.FIX_HOL MS.CUE_ON MS.MEM_PER MS.TAR_ACQ_INV MS.TAR_HOL_INV MS.TAR_ACQ MS.TAR_HOL MS.SUCCESS_ABORT MS.REWARD MS.ITI];
        observed_state = MS.TAR_ACQ_INV;
        information_state = MS.CUE_ON;
    case 4 % delay response
        all_states = [MS.INI_TRI MS.FIX_ACQ MS.FIX_HOL MS.CUE_ON MS.DEL_PER MS.TAR_ACQ MS.TAR_HOL MS.SUCCESS MS.REWARD MS.ITI];
        observed_state = MS.TAR_ACQ;
        information_state = MS.CUE_ON;
    case 5 % match-to-sample
        all_states = [MS.INI_TRI MS.FIX_ACQ MS.FIX_HOL MS.CUE_ON MS.MEM_PER MS.MAT_ACQ MS.MAT_HOL MS.SUCCESS MS.REWARD MS.ITI];
        observed_state = MS.MAT_ACQ;
        information_state = MS.CUE_ON;
    case 6 % match-to-sample masked
        all_states = [MS.INI_TRI MS.FIX_ACQ MS.FIX_HOL MS.CUE_ON MS.MEM_PER MS.MAT_ACQ_MSK MS.MAT_HOL_MSK MS.SUCCESS MS.REWARD MS.ITI];
        observed_state = MS.MAT_ACQ_MSK;
        information_state = MS.CUE_ON;
    case 9 % rotation difference, masked sample
        all_states = [MS.INI_TRI MS.FIX_ACQ MS.FIX_HOL MS.CUE_ON MS.MSK_HOL MS.MEM_PER MS.TAR_ACQ MS.TAR_HOL MS.SUCCESS MS.REWARD MS.ITI];
        observed_state = MS.TAR_ACQ;
        information_state = MS.TAR_ACQ;
    otherwise
        observed_state = NaN;
        state_1ao = NaN;
        can_infer_target = false;
        return;
end

state_1ao = all_states(find(all_states == observed_state, 1, 'first')+1);

if trial.aborted_state == -1
    can_infer_target = true;
    return;
end

if isfield(trial, 'choice') && ~isempty(trial.choice) && ~isnan(trial.choice)
    if trial.choice
        threshold_state = observed_state;
    else
        threshold_state = information_state;
    end
end

can_infer_target = find(all_states == trial.aborted_state, 1, 'first') >= find(all_states == threshold_state, 1, 'first');
end

function MS = get_ma_states()
global MA_STATES
MPA_get_states;
MS = MA_STATES;
end

%% Trial target helpers
function [T, all_eye_target_pos, all_hnd_target_pos] = fill_target_fields(T, effector_sr, observed_state, can_infer_target)
% Build eye/hand target catalogs and initial selected indices.
target_names = {'eye','hnd'};
sample_xy = {'x_eye','y_eye'; 'x_hnd','y_hnd'};
n_targets = [0 0];
all_targets = {NaN(1,0)+1i*NaN(1,0), NaN(1,0)+1i*NaN(1,0)};
selected_idx = [NaN NaN];
nonselected_idx = [NaN NaN];
for e = 1:2
    eff = target_names{e};
    if isfield(T, eff) && isfield(T.(eff), 'tar') && ~isempty(T.(eff).tar)
        n_targets(e) = numel(T.(eff).tar);
        pos_e = nan(1,n_targets(e)) + 1i*nan(1,n_targets(e));
        for k = 1:n_targets(e)
            if isfield(T.(eff).tar(k), 'pos') && numel(T.(eff).tar(k).pos) >= 2
                pos_e(k) = T.(eff).tar(k).pos(1) + 1i*T.(eff).tar(k).pos(2);
            end
        end
        all_targets{e} = unique(pos_e);
    end
    if n_targets(e) >= 1 && isfield(T, 'target_selected') && numel(T.target_selected) >= e && ~isnan(T.target_selected(e))
        candidate = round(T.target_selected(e));
        if candidate >= 1 && candidate <= n_targets(e)
            selected_idx(e) = candidate;
        end
    end
end
all_eye_target_pos = all_targets{1};
all_hnd_target_pos = all_targets{2};

% Fallback target inference from observed position.
if can_infer_target && isfield(T,'state') && ~isempty(T.state)
    obs_samples = find(T.state == observed_state);
    if ~isempty(obs_samples)
        for e = 1:2
            if ~isnan(selected_idx(e)) || n_targets(e) < 1
                continue;
            end
            x_field = sample_xy{e,1};
            y_field = sample_xy{e,2};
            if ~isfield(T, x_field) || ~isfield(T, y_field)
                continue;
            end
            observed_pos = median(T.(x_field)(obs_samples)) + 1i*median(T.(y_field)(obs_samples));
            if isfinite(real(observed_pos)) && ~any(isnan(all_targets{e}))
                [~, selected_idx(e)] = min(abs(all_targets{e} - observed_pos));
                selected_idx(e) = min(max(selected_idx(e),1),n_targets(e));
            end
        end
    end
end

% Non-selected target indices.
for e = 1:2
    if ~isnan(selected_idx(e)) && n_targets(e) >= 2
        remaining = 1:n_targets(e);
        remaining = remaining(remaining ~= selected_idx(e));
        if ~isempty(remaining)
            nonselected_idx(e) = remaining(1);
        end
    end
end

% Unified target definition:
% - saccade-only trials: use eye target
% - reach-only trials: use hand target
% - dual-effector trials: use hand/reach target for both effectors
is_saccade_expected = logical(effector_sr(1));
is_reach_expected = logical(effector_sr(2));
source_e = 1;
if is_reach_expected
    source_e = 2;
end
idx_selected_unified = selected_idx(source_e);
idx_nonselected_unified = nonselected_idx(source_e);
source_targets = all_targets{source_e};
source_count = n_targets(source_e);
tar_pos_unified = NaN + 1i*NaN;
nct_pos_unified = NaN + 1i*NaN;
if ~isnan(idx_selected_unified) && idx_selected_unified >= 1 && idx_selected_unified <= source_count
    tar_pos_unified = source_targets(idx_selected_unified);
end
if ~isnan(idx_nonselected_unified) && idx_nonselected_unified >= 1 && idx_nonselected_unified <= source_count
    nct_pos_unified = source_targets(idx_nonselected_unified);
end

if is_saccade_expected && is_reach_expected
    all_eye_target_pos = all_hnd_target_pos;
end

color_effector = target_names{source_e};
T.fix_pos = NaN + 1i*NaN;
if isfield(T, color_effector) && isfield(T.(color_effector), 'fix') && ...
        isfield(T.(color_effector).fix, 'pos') && numel(T.(color_effector).fix.pos) >= 2
    T.fix_pos = T.(color_effector).fix.pos(1) + 1i*T.(color_effector).fix.pos(2);
end

T.target_selected_idx = idx_selected_unified;
T.target_nonselected_idx = idx_nonselected_unified;
T.tar_pos = tar_pos_unified;
T.nct_pos = nct_pos_unified;
if isfield(T, 'type') && T.type == 1
    T.tar_pos = T.fix_pos;
    T.nct_pos = T.fix_pos;
end
T.hemifield = sign(real(T.tar_pos));

target_for_color = 'tar';
if     T.type == 1
    target_for_color = 'fix';
elseif T.type == 3
    target_for_color = 'cue';
end

T.col_dim = [NaN NaN NaN];
T.col_bri = [NaN NaN NaN];
T.nct_col_dim = [NaN NaN NaN];
T.nct_col_bri = [NaN NaN NaN];
if isfield(T, 'task') && isfield(T.task, color_effector) && isfield(T.task.(color_effector), target_for_color)
    target_struct = T.task.(color_effector).(target_for_color);
    if ~isempty(target_struct)
        if any(strcmp(target_for_color, {'fix','cue'}))
            source_sel = target_struct(1);
            source_nct = target_struct(1);
        else
            if ~isnan(idx_selected_unified) && idx_selected_unified >= 1 && idx_selected_unified <= numel(target_struct)
                source_sel = target_struct(idx_selected_unified);
            else
                source_sel = target_struct(1);
            end
            if ~isnan(idx_nonselected_unified) && idx_nonselected_unified >= 1 && idx_nonselected_unified <= numel(target_struct)
                source_nct = target_struct(idx_nonselected_unified);
            else
                source_nct = target_struct(1);
            end
        end
        if isfield(source_sel, 'color_dim') && numel(source_sel.color_dim) >= 3
            T.col_dim = source_sel.color_dim(1:3);
        end
        if isfield(source_sel, 'color_bright') && numel(source_sel.color_bright) >= 3
            T.col_bri = source_sel.color_bright(1:3);
        end
        if isfield(source_nct, 'color_dim') && numel(source_nct.color_dim) >= 3
            T.nct_col_dim = source_nct.color_dim(1:3);
        end
        if isfield(source_nct, 'color_bright') && numel(source_nct.color_bright) >= 3
            T.nct_col_bri = source_nct.color_bright(1:3);
        end
    end
end
end

function closest_pos = get_closest_target_pos(end_pos, all_target_pos)
closest_pos = NaN + 1i*NaN;
if ~isfinite(real(end_pos)) || isempty(all_target_pos)
    return;
end
valid_targets = all_target_pos(isfinite(real(all_target_pos)) & isfinite(imag(all_target_pos)));
if isempty(valid_targets)
    return;
end
[~, idx] = min(abs(valid_targets - end_pos));
closest_pos = valid_targets(idx);
end

%% Definition/selection helpers
function sel_idx = select_saccade_index(cfg, amplitudes, endpos_candidates, startpos_candidates, tar_pos, fix_pos)
sel_idx = 1;
n_candidates = numel(amplitudes);
if n_candidates == 0
    return;
end
all_idx = 1:n_candidates;

switch cfg.saccade_definition
    case 'closest_big_enough'
        big_idx = amplitudes >= cfg.sac_min_amp;
        if any(big_idx) && isfinite(real(tar_pos))
            distances = abs(endpos_candidates - tar_pos);
            candidate = all_idx(big_idx);
            [~, min_local] = min(distances(big_idx));
            sel_idx = candidate(min_local);
        end
    case 'biggest_close_enough'
        if isfinite(real(tar_pos))
            close_idx = abs(endpos_candidates - tar_pos) <= cfg.sac_max_off;
            if any(close_idx)
                candidate = all_idx(close_idx);
                [~, max_local] = max(amplitudes(close_idx));
                sel_idx = candidate(max_local);
            end
        end
    case 'last_in_state'
        sel_idx = n_candidates;
    case 'first_in_state'
        sel_idx = 1;
    case 'first_big_enough'
        big_idx = find(amplitudes >= cfg.sac_min_amp, 1, 'first');
        if ~isempty(big_idx)
            sel_idx = big_idx;
        end
    case 'first_end_in_target'
        if isfinite(real(tar_pos))
            idx = find(abs(endpos_candidates - tar_pos) <= cfg.sac_max_off, 1, 'first');
            if ~isempty(idx)
                sel_idx = idx;
            end
        end
    case 'first_start_in_target'
        if isfinite(real(tar_pos)) && isfinite(real(fix_pos))
            radius = abs(tar_pos - fix_pos);
            idx = find(abs(startpos_candidates - tar_pos) <= radius, 1, 'first');
            if ~isempty(idx)
                sel_idx = idx;
            end
        end
end
end

function cfg = apply_default_cfg(cfg)
defaults = struct( ...
    'interpolated_SR', 1000, ...
    'downsampling', 1, ...
    'smoothing_samples', 12, ...
    'sac_ini_t', 200, ...
    'sac_end_t', 50, ...
    'sac_min_dur', 0.03, ...
    'sac_min_amp', 2, ...
    'sac_max_off', 50, ...
    'nsacc_max', 20, ...
    'saccade_definition', 'closest_big_enough', ...
    'reach_definition', 'first_touch_in_target', ...
    'correct_offset', 1, ...
    'enable_cleaning', 1);

fields = fieldnames(defaults);
for i = 1:numel(fields)
    fn = fields{i};
    if ~isfield(cfg, fn) || isempty(cfg.(fn))
        cfg.(fn) = defaults.(fn);
    end
end
end

function trial = clean_trials(trial)
% File-level cleanup used in MPA_clean_data:
% remove last trial if empty state or duplicate time samples.
if isfield(trial(end), 'state') && isempty(trial(end).state) && numel(trial) > 1
    trial = trial(1:end-1);
end
if ~isempty(trial)
    if isfield(trial(end), 'tSample_from_time_start') && numel(trial(end).tSample_from_time_start) > 1
        if any(diff(trial(end).tSample_from_time_start) == 0) && numel(trial) > 1
            trial = trial(1:end-1);
        end
    end
end

for t = 1:numel(trial)
    % Use task target definitions if eye/hnd targets are absent.
    if (~isfield(trial(t), 'eye') || ~isfield(trial(t).eye, 'tar') || isempty(trial(t).eye.tar)) && ...
            isfield(trial(t), 'task') && isfield(trial(t).task, 'eye') && isfield(trial(t).task.eye, 'tar')
        trial(t).eye.tar = trial(t).task.eye.tar;
    end
    if (~isfield(trial(t), 'hnd') || ~isfield(trial(t).hnd, 'tar') || isempty(trial(t).hnd.tar)) && ...
            isfield(trial(t), 'task') && isfield(trial(t).task, 'hnd') && isfield(trial(t).task.hnd, 'tar')
        trial(t).hnd.tar = trial(t).task.hnd.tar;
    end

    % Remove duplicate time samples (monkeypsych branch).
    if isfield(trial(t), 'tSample_from_time_start') && numel(trial(t).tSample_from_time_start) > 1
        idx_to_remove = [false; diff(trial(t).tSample_from_time_start(:)) == 0];
        if any(idx_to_remove)
            sample_fields = {'x_eye','y_eye','x_hnd','y_hnd','state','tSample_from_time_start','sen_L','sen_R'};
            for sf = 1:numel(sample_fields)
                fn = sample_fields{sf};
                if isfield(trial(t), fn) && numel(trial(t).(fn)) == numel(idx_to_remove)
                    trial(t).(fn)(idx_to_remove) = [];
                end
            end
        end
    end
    
    
    % Repair missing used reach hand in known valid cases.
    if ~isfield(trial(t), 'reach_hand') || isempty(trial(t).reach_hand) || numel(trial(t).reach_hand)==2 
        demanded_hand = NaN;
        if isfield(trial(t), 'task') && isfield(trial(t).task, 'reach_hand') && ~isempty(trial(t).task.reach_hand)
            demanded_hand = trial(t).task.reach_hand;
        end
        if isfield(trial(t), 'effector') && isfield(trial(t), 'success') && ...
                trial(t).effector == 6 && trial(t).success == 1 && any(demanded_hand == [1 2])
            trial(t).reach_hand = demanded_hand;
        elseif numel(trial(t).reach_hand)==2 && trial(t).success==0 % this is new and bad! -> setup 1 specific? sensor input is broken there...
            
            trial(t).reach_hand = NaN;
        else
            trial(t).reach_hand = NaN;
        end        
    end
    
%     % State patch - this might actually be critical
%     if any(trial(t).states == 13) && any(trial(t).states == 12)
%         trial(t).states(trial(t).states == 12) = 14;
%     end
end
end

%% Signal processing and correction helpers
function [x_filt] = filter_et(x, flen)
if numel(x) < flen
    flen = numel(x);
end
rav = ones(1,flen)./flen;
x_filt = conv(x, rav);
if flen > 1
    start_edge_factor = flen./(ceil(flen/2):flen-1);
    start_values = x_filt(ceil(flen/2):flen-1).*start_edge_factor;
    end_edge_factor = flen./(flen-1:-1:floor(flen/2));
    end_values = x_filt(length(x_filt)+2-flen:length(x_filt)-floor(flen/2)+1).*end_edge_factor;
    x_filt = [start_values x_filt(flen:length(x_filt)-flen) end_values];
end
end

function [new_pos_x, new_pos_y] = offset_corrected(old_pos_x, old_pos_y, state_current, fixation_state_temp, fixation_x, fixation_y)
temp_pos_x = median(old_pos_x(state_current == fixation_state_temp));
temp_pos_y = median(old_pos_y(state_current == fixation_state_temp));
offset_x = abs(temp_pos_x - fixation_x);
offset_y = abs(temp_pos_y - fixation_y);

if temp_pos_x <= fixation_x
    new_pos_x = old_pos_x + offset_x;
elseif temp_pos_x > fixation_x
    new_pos_x = old_pos_x - offset_x;
else
    new_pos_x = old_pos_x;
end

if temp_pos_y <= fixation_y
    new_pos_y = old_pos_y + offset_y;
elseif temp_pos_y > fixation_y
    new_pos_y = old_pos_y - offset_y;
else
    new_pos_y = old_pos_y;
end
end
