function [synchronization_results, synchronization_report] = epp_general_synchronization_check(ephys_main_folder, behavior_main_folder)
% DAG_SYNCHRONIZATION_ALL_BLOCKS_SPLIT_SOURCES_EXAMPLE
% Example call:
%   [results, report] = DAG_synchronization_all_blocks_split_sources_example('D:\...\ephys\Bac','D:\...\behavior\Bac');
%
% This function:
% 1) Loops through day folders (YYYYMMDD) in ephys_main_folder
% 2) Loops through Block-* folders in each ephys day folder
% 3) Reads each ephys block and extracts unique RunN values
% 4) Finds behavioral file candidates in behavior_main_folder by day + run
% 5) Runs epp_synchronization for each matched run
%
% Notes:
% - Ephys and behavior files are not required to be in the same folder.
% - Matching uses subject/day prefix + parsed run number from MAT filename.

%ephys_main_folder='Y:\Data\TDTtanks\Linus_phys';
%behavior_main_folder='Y:\Data\Linus';


if nargin < 2 || isempty(ephys_main_folder) || isempty(behavior_main_folder)
    warning('Provide both ephys_main_folder and behavior_main_folder. Returning empty results.');
    synchronization_results = empty_results_struct();
    synchronization_report = empty_report_struct();
    return;
end

[~, subject_name] = fileparts(ephys_main_folder);

%% Discover and sort day folders from ephys root (YYYYMMDD)
day_dirs = dir(ephys_main_folder);
day_dirs = day_dirs([day_dirs.isdir]);
day_names = {day_dirs.name};
is_day = ~ismember(day_names, {'.', '..'}) & ~cellfun(@isempty, regexp(day_names, '^\d{8}$', 'once'));
day_dirs = day_dirs(is_day);

if isempty(day_dirs)
    warning('No day folders (YYYYMMDD) found in: %s. Returning empty results.', ephys_main_folder);
    synchronization_results = empty_results_struct();
    synchronization_report = empty_report_struct();
    return;
end

[~, day_sort_idx] = sort({day_dirs.name});
day_dirs = day_dirs(day_sort_idx);

%% Preallocate results containers
synchronization_results = empty_results_struct();
synchronization_report = empty_report_struct();

%% Process all ephys days and blocks
for d = 1:numel(day_dirs)
    day_token = day_dirs(d).name;
    ephys_day_folder = fullfile(ephys_main_folder, day_token);
    behavior_day_folder = fullfile(behavior_main_folder, day_token);
    day_token_dashed = sprintf('%s-%s-%s', day_token(1:4), day_token(5:6), day_token(7:8));
    %behavior_prefix = [subject_name day_token_dashed '_'];
    behavior_prefix = [day_token_dashed ];

    % Collect all behavior candidates for this day (day root + Block-* children).
    [behavior_files, behavior_run_numbers] = collect_behavior_day_candidates(behavior_day_folder, behavior_prefix);

    block_dirs = dir(fullfile(ephys_day_folder, 'Block-*'));
    block_dirs = block_dirs([block_dirs.isdir]);
    if isempty(block_dirs)
        day_warning = sprintf('[%s] No Block-* folders found in ephys day folder.', day_token);
        warning('%s', day_warning);

        result_entry = initialize_result_entry(ephys_day_folder, '', '', NaN);
        report_entry = initialize_report_entry(ephys_day_folder, '', '', NaN);
        report_entry.messages = append_report(report_entry.messages, day_warning);
        synchronization_results(end+1) = result_entry; 
        synchronization_report(end+1) = report_entry; 
        continue;
    end

    block_numbers = nan(numel(block_dirs), 1);
    for i = 1:numel(block_dirs)
        block_numbers(i) = parse_block_number(block_dirs(i).name);
    end
    [~, block_sort_idx] = sort(block_numbers);
    block_dirs = block_dirs(block_sort_idx);

    for i = 1:numel(block_dirs)
        block_name = block_dirs(i).name;
        block_path = fullfile(ephys_day_folder, block_name);

        % Read ephys first, then use unique RunN values as requested.
        try
            ephys_data = TDTbin2mat_working(block_path, 'EXCLUSIVELYREAD', {'SVal','Tnum','RunN','Sess'});
        catch eeee
            msg = sprintf('[%s/%s] Could not read ephys. probably no TSQ file?', day_token, block_name);
            warning('%s', msg);
            result_entry = initialize_result_entry(ephys_day_folder, block_path, '', NaN);
            report_entry = initialize_report_entry(ephys_day_folder, block_path, '', NaN);
            report_entry.messages = append_report(report_entry.messages, msg);
            synchronization_results(end+1) = result_entry;
            synchronization_report(end+1) = report_entry;
            continue;
        end

        if ~isfield(ephys_data,'epocs') || ~isfield(ephys_data.epocs,'RunN') || ~isfield(ephys_data.epocs.RunN,'data') || isempty(ephys_data.epocs.RunN.data)
            msg = sprintf('[%s/%s] Missing RunN in ephys data.', day_token, block_name);
            warning('%s', msg);
            result_entry = initialize_result_entry(ephys_day_folder, block_path, '', NaN);
            report_entry = initialize_report_entry(ephys_day_folder, block_path, '', NaN);
            report_entry.messages = append_report(report_entry.messages, msg);
            synchronization_results(end+1) = result_entry; 
            synchronization_report(end+1) = report_entry; 
            continue;
        end

        run_candidates = unique(ephys_data.epocs.RunN.data(2:end));
        run_candidates = run_candidates(~isnan(run_candidates));
        run_candidates = run_candidates(:)';
        if isempty(run_candidates)
            msg = sprintf('[%s/%s] RunN present but no valid run numbers.', day_token, block_name);
            warning('%s', msg);
            result_entry = initialize_result_entry(ephys_day_folder, block_path, '', NaN);
            report_entry = initialize_report_entry(ephys_day_folder, block_path, '', NaN);
            report_entry.messages = append_report(report_entry.messages, msg);
            synchronization_results(end+1) = result_entry; 
            synchronization_report(end+1) = report_entry; 
            continue;
        end

        for run_number = run_candidates
            result_entry = initialize_result_entry(ephys_day_folder, block_path, '', run_number);
            report_entry = initialize_report_entry(ephys_day_folder, block_path, '', run_number);

            matched_idx = find(behavior_run_numbers == run_number);
            if isempty(matched_idx)
                msg = sprintf('[%s/%s] No behavior MAT candidate found for run %d (prefix %s*.mat).', day_token, block_name, run_number, behavior_prefix);
                report_entry.messages = append_report(report_entry.messages, msg);
                warning('%s', msg);
                synchronization_results(end+1) = result_entry; 
                synchronization_report(end+1) = report_entry; 
                continue;
            end

            matched_files = sort(behavior_files(matched_idx));
            behavior_file = matched_files{1};
            result_entry.behavior_file = behavior_file;
            report_entry.behavior_file = behavior_file;

            if numel(matched_files) > 1
                [~, chosen_name, chosen_ext] = fileparts(behavior_file);
                msg = sprintf('[%s/%s] Multiple behavior candidates for run %d. Using: %s%s', day_token, block_name, run_number, chosen_name, chosen_ext);
                report_entry.messages = append_report(report_entry.messages, msg);
                warning('%s', msg);
            end

            behavioral_data = load(behavior_file, 'trial');
            behavioral_data.run = run_number;
            [sync_data, report] = epp_synchronization(ephys_data, behavioral_data);

            result_entry.continuous_timestamps = sync_data.continuous_timestamps;
            result_entry.Trial_timestamps = sync_data.Trial_timestamps;
            report_entry.messages = [report_entry.messages; report];

            synchronization_results(end+1) = result_entry; 
            synchronization_report(end+1) = report_entry; 
        end
    end
end

%% Optional quick summary
num_fail = sum(arrayfun(@has_failure, synchronization_report));
num_ok = numel(synchronization_results) - num_fail;
fprintf('Done. %d ephys days scanned, %d block/run entries processed: %d successful, %d failed.\n', numel(day_dirs), numel(synchronization_results), num_ok, num_fail);

%% Save one-line-per-entry report text file in ephys_main_folder
report_file = fullfile(ephys_main_folder, 'synchronization_report_split_sources.txt');
write_text_report(report_file, synchronization_report);
result_file = fullfile(ephys_main_folder, 'synchronization_info');
save(result_file,'synchronization_results');
fprintf('Text report saved: %s\n', report_file);
end

% -------------------------- Local helper functions --------------------------
function n = parse_block_number(block_name)
token = regexp(block_name, '^Block-(\d+)$', 'tokens', 'once');
if isempty(token)
    n = inf;
else
    n = str2double(token{1});
end
end

function run_number = parse_run_number(file_name)
token = regexp(file_name, '_(\d+)\.mat$', 'tokens', 'once');
if isempty(token)
    run_number = NaN;
else
    run_number = str2double(token{1});
end
end

function [files_out, runs_out] = collect_behavior_day_candidates(day_folder, prefix)
files_out = {};
runs_out = [];
if ~isfolder(day_folder)
    return;
end

% Candidates directly in day folder.
root_candidates = dir(fullfile(day_folder, ['*' prefix '*.mat']));
for j = 1:numel(root_candidates)
    current_path = fullfile(day_folder, root_candidates(j).name);
    files_out{end+1,1} = current_path; %#ok<AGROW>
    runs_out(end+1,1) = parse_run_number(root_candidates(j).name); %#ok<AGROW>
end

% Candidates inside behavior Block-* folders for same day.
behavior_block_dirs = dir(fullfile(day_folder, 'Block-*'));
behavior_block_dirs = behavior_block_dirs([behavior_block_dirs.isdir]);
for b = 1:numel(behavior_block_dirs)
    current_block = fullfile(day_folder, behavior_block_dirs(b).name);
    block_candidates = dir(fullfile(current_block, [prefix '*.mat']));
    for j = 1:numel(block_candidates)
        current_path = fullfile(current_block, block_candidates(j).name);
        files_out{end+1,1} = current_path; %#ok<AGROW>
        runs_out(end+1,1) = parse_run_number(block_candidates(j).name); %#ok<AGROW>
    end
end

% Keep only candidates with valid run numbers.
valid_idx = ~isnan(runs_out);
files_out = files_out(valid_idx);
runs_out = runs_out(valid_idx);
end

function report = append_report(report, msg)
%report(end+1,1) = string(msg);
report = [report '|' msg];
end

function tf = has_failure(report_entry)
if isempty(report_entry.messages)
    tf = false;
    return;
end
report_text = lower(strjoin(cellstr(report_entry.messages), ' '));
failure_markers = {'corrupted', 'failed', 'no behavioral mat file', 'no block-* folders found', 'could not parse run number'};
for f=1:numel(failure_markers)
    tf = any(strfind(report_text, failure_markers{f}));
    if tf
        break;
    end
end
end

function entry = initialize_result_entry(day_folder, block_folder, behavior_file, run_number)
entry.day_folder = day_folder;
entry.block_folder = block_folder;
entry.behavior_file = behavior_file;
entry.run_number = run_number;
entry.continuous_timestamps = [];
entry.Trial_timestamps = [];
end

function entry = initialize_report_entry(day_folder, block_folder, behavior_file, run_number)
entry.day_folder = day_folder;
entry.block_folder = block_folder;
entry.behavior_file = behavior_file;
entry.run_number = run_number;
%entry.messages = strings(0,1);
entry.messages = '';
end

function s = empty_results_struct()
s = struct( ...
    'day_folder', {}, ...
    'block_folder', {}, ...
    'behavior_file', {}, ...
    'run_number', {}, ...
    'continuous_timestamps', {}, ...
    'Trial_timestamps', {} );
end

function s = empty_report_struct()
s = struct( ...
    'day_folder', {}, ...
    'block_folder', {}, ...
    'behavior_file', {}, ...
    'run_number', {}, ...
    'messages', {} );
end

function write_text_report(report_file, synchronization_report)
fid = fopen(report_file, 'wt');
if fid == -1
    warning('Could not create report file: %s', report_file);
    return;
end
cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

for i = 1:numel(synchronization_report)
    entry = synchronization_report(i);
    output_folder = entry.block_folder;
    if isempty(output_folder)
        output_folder = entry.day_folder;
    end

    [~, behavior_name, behavior_ext] = fileparts(entry.behavior_file);
    if isempty(behavior_name)
        behavior_file_name = 'N/A';
    else
        behavior_file_name = [behavior_name behavior_ext];
    end

    if isempty(entry.messages)
        message_text = 'OK';
    else
        message_text = strjoin(cellstr(entry.messages), ' | ');
    end

    message_text = strrep(message_text, sprintf('\r'), ' ');
    message_text = strrep(message_text, sprintf('\n'), ' ');
    if isnan(entry.run_number)
        run_text = 'run=N/A';
    else
        run_text = sprintf('run=%d', entry.run_number);
    end

    fprintf(fid, '%s - %s - %s - %s\n', output_folder, behavior_file_name, run_text, message_text);
end
end
