function epp_prepare_blocks(cfg)
%EPP_PREPARE_BLOCKS  Stage 1: sync behavior to TDT blocks and save payloads.
%
%   epp_prepare_blocks(cfg)
%
%   Inputs:
%     cfg - struct from epp_load_cfg. Uses cfg.monkeys, cfg.datasets,
%           cfg.movements (passed to MP_add_saccades_and_reaches),
%           cfg.roots.sorting_tables, ephys_tanks, behavior, Trials.
%
%   Outputs:
%     none. Writes cfg.unit_info_file (sorting-table units, plus .monkey)
%     and, into cfg.roots.Trials:
%       <monkey>_<YYYYMMDD>_Block-<NNN>.mat
%                                  - variable block_payload with fields
%                                    monkey, session, block, trial, state, sample
%       prepare_blocks_report.txt  - one line per run: ephys path, behavior file,
%                                    run number, sync messages
%
%   block_payload.trial   - behavior trials kept by sync, with .timestamp set to
%                           the ephys time of state 2. Sample-wise fields
%                           (state, x_hnd, y_hnd, x_eye, y_eye, tSample_...) removed.
%   block_payload.state   - .Value and .timestamps (ephys time) from behavior states.
%   block_payload.sample  - x_hnd, y_hnd, x_eye, y_eye, timestamps (ephys time).
%   Per-trial sample vectors must be columns; timestamps are zipped by position,
%   not by trial number.
%
% Stage 1: build unit_info and save processed block payloads.

epp_ensure_paths();
unit_info = struct([]);
report = repmat(struct('ephys_block', '', 'behavior_file', '', 'run', NaN, 'messages', ''), 0, 1);

fields_to_remove = {'tSample_from_time_start','trial_number','state','sen_L','sen_R','jaw','body','states','states_onset','eye','hnd'};
per_sample_fields = {'x_hnd','y_hnd','x_eye','y_eye'};
per_sample_struct=[per_sample_fields' repmat({[]},numel(per_sample_fields),1)]';


for m = 1:numel(cfg.monkeys)
    monkey = cfg.monkeys{m};
    monkey_units = epp_sorting_table_to_units(monkey, cfg);
    [monkey_units.monkey] = deal(monkey);
    unit_info = [unit_info; monkey_units]; %#ok<AGROW>

    % One row per (session, block) that any unit of this monkey actually is found in   
    sessions_blocks=zeros(0,2);
    for i = 1:numel(monkey_units)
        pair = [repmat(monkey_units(i).Session, numel(monkey_units(i).Blocks), 1), monkey_units(i).Blocks(:)];
        sessions_blocks = union(sessions_blocks, pair, 'rows');
    end

    for b = 1:size(sessions_blocks,1)
        session = sessions_blocks(b,1);
        block = sessions_blocks(b,2);
        day_token = sprintf('%08d', session);
        ephys_block_path = fullfile(cfg.roots.ephys_tanks, [monkey '_phys'], day_token, sprintf('Block-%d', block));

        ephys_data = TDTbin2mat_working(ephys_block_path, 'EXCLUSIVELYREAD', {'SVal','Tnum','RunN','Sess'});       

        if isfield(ephys_data, 'epocs') && isfield(ephys_data.epocs, 'RunN') && isfield(ephys_data.epocs.RunN, 'data') && numel(ephys_data.epocs.RunN.data)>1
            runs = ephys_data.epocs.RunN.data;
            runs = unique(runs(2:end)); % skipping first potentially corrupted value
            runs = runs(:)';
        else
            report(end+1) = make_report_entry(ephys_block_path, '', NaN, 'No run candidates in ephys block'); %#ok<AGROW>
            continue;
        end

        block_payload = struct('monkey', monkey, 'session', session, 'block', block, ...
            'trial', struct([]), ...
            'state', struct('Value', [], 'timestamps', []), ...
            'sample', struct(per_sample_struct{:}, 'timestamps', []));

        for r = 1:numel(runs)
            run = runs(r);
            fname = sprintf('%s%s-%s-%s_%02d.mat', monkey(1:3), day_token(1:4), day_token(5:6), day_token(7:8), run);
            behavior_file = fullfile(cfg.roots.behavior, monkey, day_token, fname);

            if ~isfile(behavior_file)
                report(end+1) = make_report_entry(ephys_block_path, behavior_file, run, 'Behavior file missing'); %#ok<AGROW>
                continue;
            end

            behavioral_data = load(behavior_file, 'trial');
            Trials = MP_add_saccades_and_reaches(behavioral_data.trial, cfg.movements);
            [sync_data, sync_report] = epp_synchronization(ephys_data, Trials, run);

            run_messages = sync_report;
            if isempty(sync_data.Trial_timestamps)                
                run_messages = [run_messages ' | Synchronization returned no trial timestamps']; %#ok<AGROW>
            end
            report(end+1) = make_report_entry(ephys_block_path, behavior_file, run, run_messages); %#ok<AGROW>
            
            keep_trials = ismember([Trials.n], sync_data.retained_behavior_trial_numbers);
            Trials = Trials(keep_trials);
            timestamps = num2cell(sync_data.Trial_timestamps);
            [Trials.timestamp]=timestamps{:};

            existing_fields = fieldnames(Trials);
            remove_fields = intersect([fields_to_remove per_sample_fields], existing_fields);
            trial_export = rmfield(Trials, remove_fields);

            block_payload.trial              = [block_payload.trial; trial_export(:)];                                  % #ok<AGROW>
            block_payload.state.Value        = [block_payload.state.Value; sync_data.state_values(:)];                  % #ok<AGROW>
            block_payload.state.timestamps   = [block_payload.state.timestamps; sync_data.state_onsets(:)];             % #ok<AGROW>
            block_payload.sample.timestamps  = [block_payload.sample.timestamps; sync_data.continuous_timestamps(:)];   % #ok<AGROW>
            
            for f = 1:numel(per_sample_fields)
                FN = per_sample_fields{f};
                block_payload.sample.(FN)       = [block_payload.sample.(FN); vertcat(Trials.(FN))];                    % #ok<AGROW>
            end
        end

        if isempty(block_payload.trial)
            report(end+1) = make_report_entry(ephys_block_path, '', NaN, 'Block saved with no trial data');             % #ok<AGROW>
        end

        output_name = sprintf('%s_%s_Block-%03d.mat', monkey, day_token, block);
        save(fullfile(cfg.roots.processed_trials, output_name), 'block_payload');
    end
end

save(cfg.unit_info_file, 'unit_info');
write_prepare_blocks_report(fullfile(cfg.roots.processed_trials, 'prepare_blocks_report.txt'), report);
end

function entry = make_report_entry(ephys_block, behavior_file, run_number, messages)
entry = struct('ephys_block', ephys_block, 'behavior_file', behavior_file, 'run', run_number, 'messages', messages);
end

function write_prepare_blocks_report(report_file, report)
fid = fopen(report_file, 'wt');
if fid == -1
    warning('Could not create report file: %s', report_file);
    return;
end
cleanup_obj = onCleanup(@() fclose(fid)); % #ok<NASGU>

for i = 1:numel(report)
    entry = report(i);
    [~, behavior_name, behavior_ext] = fileparts(entry.behavior_file);
    if isempty(behavior_name)
        behavior_file_name = 'N/A';
    else
        behavior_file_name = [behavior_name behavior_ext];
    end

    if isempty(entry.messages)
        message_text = 'OK';
    else
        message_text = entry.messages;
    end
    message_text = strrep(message_text, sprintf('\r'), ' ');
    message_text = strrep(message_text, sprintf('\n'), ' ');

    if isnan(entry.run)
        run_text = 'run=N/A';
    else
        run_text = sprintf('run=%d', entry.run);
    end

    fprintf(fid, '%s - %s - %s - %s\n', entry.ephys_block, behavior_file_name, run_text, message_text);
end
end
