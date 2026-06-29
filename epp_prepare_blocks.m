function epp_prepare_blocks(cfg)
% Stage 1: build unit_info and save processed block payloads.

unit_info = struct([]);
report = repmat(struct('ephys_block', '', 'behavior_file', '', 'run', NaN, 'messages', ''), 0, 1);

for m = 1:numel(cfg.monkeys)
    monkey = cfg.monkeys{m};
    monkey_units = epp_sorting_table_to_units(monkey, cfg);
    [monkey_units.monkey] = deal(monkey);
    unit_info = [unit_info; monkey_units]; %#ok<AGROW>

    relevant_blocks = epp_get_relevant_blocks(monkey_units);
    for b = 1:size(relevant_blocks,1)
        session = relevant_blocks(b,1);
        block = relevant_blocks(b,2);
        day_token = sprintf('%08d', session);
        ephys_block_path = fullfile(cfg.roots.ephys_tanks, [monkey '_phys'], day_token, sprintf('Block-%d', block));

        if ~isfolder(ephys_block_path)
            report(end+1) = make_report_entry(ephys_block_path, '', NaN, 'Ephys block folder missing'); %#ok<AGROW>
            continue;
        end
        ephys_data = TDTbin2mat_working(ephys_block_path, 'EXCLUSIVELYREAD', {'SVal','Tnum','RunN','Sess'});
        run_candidates = epp_get_run_candidates(ephys_data);
        if isempty(run_candidates)
            report(end+1) = make_report_entry(ephys_block_path, '', NaN, 'No run candidates in ephys block'); %#ok<AGROW>
            continue;
        end

        block_payload = struct('monkey', monkey, 'session', session, 'block', block, ...
            'trial', struct([]), ...
            'state', struct('Value', [], 'timestamps', []), ...
            'sample', struct('x_hnd', [], 'y_hnd', [], 'x_eye', [], 'y_eye', [], 'timestamps', []));

        for r = 1:numel(run_candidates)
            run_number = run_candidates(r);
            fname = sprintf('%s%s-%s-%s_%02d.mat', monkey(1:3), day_token(1:4), day_token(5:6), day_token(7:8), run_number);
            behavior_file = fullfile(cfg.roots.behavior, monkey, day_token, fname);

            if ~isfile(behavior_file)
                report(end+1) = make_report_entry(ephys_block_path, behavior_file, run_number, 'Behavior file missing'); %#ok<AGROW>
                continue;
            end

            behavioral_data = load(behavior_file, 'trial');
            processed_run.trial = MP_add_saccades_and_reaches(behavioral_data.trial);
            processed_run.run = run_number;
            [sync_data, sync_report] = epp_synchronization(ephys_data, processed_run);

            run_messages = sync_report;
            if isempty(sync_data.Trial_timestamps)
                run_messages = append_report_message(run_messages, 'Synchronization returned no trial timestamps');
            end
            report(end+1) = make_report_entry(ephys_block_path, behavior_file, run_number, run_messages); %#ok<AGROW>

            run_chunk = epp_build_run_chunk(processed_run, sync_data);
            block_payload = epp_concatenate_aligned_data(block_payload, run_chunk);
        end

        if isempty(block_payload.trial)
            report(end+1) = make_report_entry(ephys_block_path, '', NaN, 'Block saved with no trial data'); %#ok<AGROW>
        end

        output_name = sprintf('%s_%s_Block-%03d.mat', monkey, day_token, block);
        save(fullfile(cfg.roots.processed_trials, output_name), 'block_payload');
    end
end

save(fullfile(cfg.roots.processed_trials, 'unit_info.mat'), 'unit_info');
write_prepare_blocks_report(fullfile(cfg.roots.processed_trials, 'prepare_blocks_report.txt'), report);
end

function relevant_blocks = epp_get_relevant_blocks(units)
relevant_blocks = zeros(0,2);
for i = 1:numel(units)
    session = units(i).Session;
    blocks = units(i).Blocks;
    for b = blocks
        relevant_blocks(end+1,:) = [session b]; %#ok<AGROW>
    end
end
relevant_blocks = unique(relevant_blocks, 'rows');
relevant_blocks = sortrows(relevant_blocks, [1 2]);
end

function runs = epp_get_run_candidates(ephys_data)
runs = [];
if ~isfield(ephys_data, 'epocs') || ~isfield(ephys_data.epocs, 'RunN') || ...
        ~isfield(ephys_data.epocs.RunN, 'data') || isempty(ephys_data.epocs.RunN.data)
    return;
end
if numel(ephys_data.epocs.RunN.data)>1 && ephys_data.epocs.RunN.data(1) ~= ephys_data.epocs.RunN.data(2)
    ephys_data.epocs.RunN.data(1) = ephys_data.epocs.RunN.data(2);
end

runs = unique(ephys_data.epocs.RunN.data(:));
runs = runs(isfinite(runs));
runs = runs(:)';
end

function run_chunk = epp_build_run_chunk(processed_run, sync_data)

keep_trials = ismember([processed_run.trial.n], sync_data.retained_behavior_trial_numbers);
trial_processed = processed_run.trial(keep_trials);
n_assign = min(numel(trial_processed), numel(sync_data.Trial_timestamps));
for tt = 1:numel(trial_processed)
    trial_processed(tt).timestamp = NaN;
end
for tt = 1:n_assign
    trial_processed(tt).timestamp = sync_data.Trial_timestamps(tt);
end

per_sample_fields = {'x_hnd','y_hnd','x_eye','y_eye'};
sample_run = struct();
for f = 1:numel(per_sample_fields)
    FN = per_sample_fields{f};
    sample_run.(FN) = vertcat(trial_processed.(FN));
end
sample_run.timestamps = sync_data.continuous_timestamps;

fields_to_remove = {'tSample_from_time_start','trial_number','state','sen_L','sen_R','jaw','body','states','states_onset','eye','hnd'};
existing_fields = fieldnames(trial_processed);
remove_fields = intersect(fields_to_remove, existing_fields);
trial_export = rmfield(trial_processed, remove_fields);

run_chunk = struct( ...
    'trial', trial_export, ...
    'state', struct('Value', sync_data.state_values(:), 'timestamps', sync_data.state_onsets(:)), ...
    'sample', sample_run);
end

function entry = make_report_entry(ephys_block, behavior_file, run_number, messages)
entry = struct('ephys_block', ephys_block, 'behavior_file', behavior_file, ...
    'run', run_number, 'messages', messages);
end

function messages = append_report_message(messages, msg)
if isempty(messages)
    messages = msg;
else
    messages = [messages ' | ' msg];
end
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
