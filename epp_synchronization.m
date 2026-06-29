function [sync_data, report]=epp_synchronization(ephys_data,behavioral_data,debug_on)
%epp_synchronization Align behavioral timing to ephys block time.
%   This function maps timing from behavioral trial structure to timestamps
%   relative to the start of a given TDT ephys block.
%
%   Stream mapping used in this function:
%   - Behavioral stream: behavioral_data.trial(t).state and
%     behavioral_data.trial(t).tSample_from_time_start
%   - Ephys trial/run stream: ephys_data.epocs.Tnum and ephys_data.epocs.RunN
%   - Ephys state stream: ephys_data.epocs.SVal

% INPUTS:
% behavioral_data - struct containing:
%                   behavioral_data.trial (trial structure array)
%                   behavioral_data.run   (run number)
% ephys_data      - struct created by reading a TDT block (TDTbin2mat or
%                   TDTbin2mat_working)

% OUTPUTS:
% sync_data             - struct with synchronization outputs:
%                         .continuous_timestamps
%                         .Trial_timestamps
%                         .retained_trial_numbers
%                         .retained_behavior_trial_numbers
%                         .state_onsets
%                         .state_timestamps
%                         .state_values
% report                - string array with all displayed messages in order

if nargin < 3 || isempty(debug_on)
    debug_on = true;
end

% Initialize outputs for safe early returns.
sync_data = struct( ...
    'continuous_timestamps', [], ...
    'Trial_timestamps', [], ...
    'retained_trial_numbers', [], ...
    'retained_behavior_trial_numbers', [], ...
    'state_onsets', [], ...
    'state_timestamps', [], ...
    'state_values', []);

% Collect all displayed messages for optional downstream logging/QA.
report = '';

% Early input sanity checks: behavioral/ephys trial and block identifiers
% must be available and non-empty.
if ~isfield(behavioral_data,'trial') || isempty(behavioral_data.trial)
    log_message('Behavioral trial data is missing or empty.');
    return
end
if ~isfield(behavioral_data,'run') || isempty(behavioral_data.run)
    log_message('Behavioral run/block identifier is missing or empty.');
    return
end
if ~isfield(ephys_data,'epocs') || ~isfield(ephys_data.epocs,'Tnum') || ~isfield(ephys_data.epocs.Tnum,'data') || isempty(ephys_data.epocs.Tnum.data)
    log_message('Ephys trial numbers (Tnum) are missing or empty.');
    return
end
if ~isfield(ephys_data,'epocs') || ~isfield(ephys_data.epocs,'RunN') || ~isfield(ephys_data.epocs.RunN,'data') || isempty(ephys_data.epocs.RunN.data)
    log_message('Ephys run/block identifiers (RunN) are missing or empty.');
    return
end
if ~isfield(ephys_data,'epocs') || ~isfield(ephys_data.epocs,'SVal') || ...
        ~isfield(ephys_data.epocs.SVal,'onset') || isempty(ephys_data.epocs.SVal.onset) || ...
        ~isfield(ephys_data.epocs.SVal,'data') || isempty(ephys_data.epocs.SVal.data)
    log_message('Ephys state stream (SVal) is missing or empty.');
    return
end

% Extract behavioral and ephys metadata used for alignment
behavior_trials = behavioral_data.trial;
behavior_trial_numbers = [behavior_trials.n];
behavior_run_number = behavioral_data.run;
ephys_trial_numbers = ephys_data.epocs.Tnum.data;
ephys_run_numbers = ephys_data.epocs.RunN.data;
ephys_trial_onsets = [ephys_data.epocs.Tnum.onset];


% Compute sample-wise and trial-wise timestamps in one pass
% (same state-2 alignment anchor used for both outputs).
ephys_state_onsets = ephys_data.epocs.SVal.onset;
ephys_state_values = ephys_data.epocs.SVal.data;

if debug_on 
    % Optional handling for known historical acquisition anomalies:
    
    Session     =ephys_data.epocs.Sess.data;
    % Correct trial/run counters if the first trial was initialized incorrectly.
    if numel(ephys_trial_numbers)>1 && ephys_trial_numbers(1)~=1
        log_message('First incorrectly initialized trial corrected');
        ephys_trial_numbers(1) = ephys_trial_numbers(2)-1;
        ephys_run_numbers(1) = ephys_run_numbers(2);
        ephys_trial_onsets(1) = 0;
        Session(1)  = Session(2);
    end
    
    if ephys_state_values(end)==1 || (ephys_trial_onsets(end) >= ephys_state_onsets(find(ephys_state_values==2,1,'last')))
        log_message('Last incorrectly initialized trial corrected');
        ephys_trial_numbers(end) = [];
        ephys_run_numbers(end) = [];
        ephys_trial_onsets(end) = [];
        Session(end)  = [];
    end
    
    % Remove one spurious initial trial in a known legacy recording issue (Linus_20150703, Block-5).
    if numel(ephys_trial_numbers)>1 && ephys_trial_numbers(2)==1
        N_to_remove=find(diff(find(ephys_trial_numbers==1))==1,1,'last');
        
        ephys_trial_numbers(1:N_to_remove) = [];
        ephys_run_numbers(1:N_to_remove) = [];
        Session(1:N_to_remove)      =[];
        ephys_trial_onsets(1:N_to_remove) = [];
        log_message(['Additional ' mat2str(N_to_remove) ' trial(s) in the beginning removed']);
    end
    
    % If multiple initial trial counters are invalid, reject this block.
    % Example known case: Bac_20210826.
    if numel(ephys_trial_numbers)>1 && (any(ephys_trial_numbers<1) || any( Session<100000 | Session>800000) || any(Session~=Session(end)))
        log_message('Synchronization impossible due to corrupted ephys state information - entire run invalid');
        return;
    end
        
    % Remove ephys trial onsets that do not correspond to the behavioral run.
    % This can occur when an ephys block spans multiple behavioral runs.
    if any(ephys_run_numbers~=behavior_run_number)
        log_message(['Warning: multiple runs in one block! Run onsets at TDT trials: ' mat2str(find(ephys_trial_numbers==1))]);
        matching_trials=ephys_run_numbers==behavior_run_number;
        ephys_trial_numbers = ephys_trial_numbers(matching_trials);
        ephys_trial_onsets = ephys_trial_onsets(matching_trials);
        log_message(['Retained TDT trials corresponding to the requested behavioral run: ' mat2str(find(matching_trials))]);
    end
    
    % Keep only behavioral trials that are in the ephys data as well - rare case
    % of a last behavioral trial was initiated, but not streamed to ephys
    if numel(behavior_trial_numbers)>numel(ephys_trial_numbers)
        overlapping_trials=arrayfun(@(x) any(ephys_trial_numbers==x),behavior_trial_numbers);
        behavior_trial_numbers = behavior_trial_numbers(overlapping_trials);
        log_message(['Too many behavioral trials: ' mat2str(numel(overlapping_trials) - sum(overlapping_trials)) ' behavioral trial(s) removed']);
    end
    
end

continuous_timestamps=[];
Trial_timestamps=[];
behavior_state_onsets_all = [];
behavior_state_values_all = [];
for t=behavior_trial_numbers
    % Align to state 2 (not state 1), because state 1 is used for trial
    % initiation/signaling and has different onset timing properties.
    behavior_state2_time = behavior_trials(t).tSample_from_time_start(find(behavior_trials(t).state==2,1));
    ephys_trial_onset = ephys_trial_onsets(find(ephys_trial_numbers==t, 1, 'first'));
    ephys_state2_onset = ephys_state_onsets(find(ephys_state_values==2 & ephys_state_onsets>ephys_trial_onset, 1, 'first'));
    behavior_timestamps_in_ephys_time = behavior_trials(t).tSample_from_time_start - behavior_state2_time + ephys_state2_onset;
    continuous_timestamps=[continuous_timestamps; behavior_timestamps_in_ephys_time];
    Trial_timestamps=[Trial_timestamps; ephys_state2_onset];

    behavior_state_values_all = [behavior_state_values_all; behavior_trials(t).states(:)]; %#ok<AGROW>
    behavior_state_onsets_all = [behavior_state_onsets_all; behavior_trials(t).states_onset(:) - behavior_state2_time + ephys_state2_onset]; %#ok<AGROW>
end

sync_data.continuous_timestamps = continuous_timestamps;
sync_data.Trial_timestamps = Trial_timestamps;
sync_data.retained_trial_numbers = ephys_trial_numbers(:).';
sync_data.retained_behavior_trial_numbers = behavior_trial_numbers(:).';
sync_data.state_onsets = behavior_state_onsets_all;
sync_data.state_timestamps = behavior_state_onsets_all;
sync_data.state_values = behavior_state_values_all;

    function log_message(msg)
        % Keep legacy console output behavior while storing a report copy.
        disp(msg);
        report = [report ' | ' msg];
    end
end
