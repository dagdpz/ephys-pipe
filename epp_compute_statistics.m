function epp_compute_statistics(cfg)
% Stage 3: per-unit epoch firing rates and statistical comparisons.

load(fullfile(cfg.roots.processed_trials, 'unit_info.mat'), 'unit_info');
processed_blocks = epp_load_saved_processed_blocks(cfg.roots.processed_trials);

unit_statistics = repmat(struct('Neuron_ID', '', 'comparisons', struct()), numel(unit_info), 1);

for u = 1:numel(unit_info)
    unit_statistics(u).Neuron_ID = unit_info(u).Neuron_ID;
    unit_combined = load_unit_combined(unit_info(u), processed_blocks, cfg);
    if isempty(unit_combined.trial)
        continue;
    end

    trial_rates = compute_trial_epoch_rates(unit_combined, cfg);
    unit_statistics(u).comparisons = build_unit_comparisons(trial_rates, unit_combined, cfg);
end

statistics_table = unit_statistics_to_table(unit_statistics);
save(fullfile(cfg.roots.statistics, 'unit_statistics.mat'), 'unit_statistics', 'statistics_table');
writetable(statistics_table, fullfile(cfg.roots.statistics, 'unit_statistics.xlsx'));
end

function unit_combined = load_unit_combined(unit_row, processed_blocks, cfg)
session = unit_row.Session;
blocks = unit_row.Blocks;
monkey = unit_row.monkey;
channel = unit_row.Channel;
filenumber = unit_row.Filenumber;
sortcode = unit_row.Unit;
threshold = 'negthr';

day_token = sprintf('%08d', session);
session_folder = fullfile(cfg.roots.ephys_tanks, [monkey '_phys'], day_token);
spike_file = fullfile(session_folder, sprintf('dataspikes_rb%03d_ch%03d_%s.mat', filenumber, channel, threshold));
spike_data = load(spike_file, 'cluster_class', 'par');
cc = spike_data.cluster_class;
par = spike_data.par;

segment_ends = par.segmentends(:) ./ par.sr;
segment_starts = [0; segment_ends(1:end-1)];
spike_times = cc(cc(:,1) == sortcode, 2) / 1000;

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
end

function trial_rates = compute_trial_epoch_rates(unit_combined, cfg)
n_trials = numel(unit_combined.trial);
n_epochs = numel(cfg.EPOCHS);
trial_rates = NaN(n_trials, n_epochs);
spike_times = unit_combined.spike_times;

trial_starts = [unit_combined.trial.timestamp].';
if n_trials > 1
    trial_ends = [trial_starts(2:end); inf];
else
    trial_ends = inf;
end

for ti = 1:n_trials
    for ei = 1:n_epochs
        trial_rates(ti, ei) = trial_epoch_rate( ...
            unit_combined, spike_times, trial_starts(ti), trial_ends(ti), cfg.EPOCHS(ei));
    end
end
end

function rate_hz = trial_epoch_rate(unit_combined, spike_times, trial_start, trial_end, epoch)
state_vals = unit_combined.state.Value;
state_times = unit_combined.state.timestamps;
align_mask = (state_vals == epoch.align_state) & ...
    state_times >= trial_start & state_times < trial_end;
align_times = state_times(align_mask);

if isempty(align_times)
    rate_hz = NaN;
    return;
end

event_rates = arrayfun(@(align_time) mean_rate_in_window( ...
    spike_times, align_time, epoch.t_start_s, epoch.t_end_s), align_times);
rate_hz = mean(event_rates, 'omitnan');
end

function rate_hz = mean_rate_in_window(spike_times, align_time, t_start_s, t_end_s)
win_start = align_time + t_start_s;
win_end = align_time + t_end_s;
duration_s = win_end - win_start;

if ~isfinite(align_time) || duration_s <= 0
    rate_hz = NaN;
    return;
end

n_spikes = sum(spike_times >= win_start & spike_times < win_end);
rate_hz = n_spikes / duration_s;
end

function comparisons = build_unit_comparisons(trial_rates, unit_combined, cfg)
comparisons = struct();
n_trials = size(trial_rates, 1);
n_epochs = numel(cfg.EPOCHS);
n_conditions = numel(cfg.CONDITIONS);

condition_masks = false(n_trials, n_conditions);
for ci = 1:n_conditions
    for ti = 1:n_trials
        condition_masks(ti, ci) = trial_matches_condition(unit_combined.trial(ti), cfg.CONDITIONS(ci));
    end
end

epoch_idx_by_name = containers.Map({cfg.EPOCHS.name}, num2cell(1:n_epochs));

for ci = 1:n_conditions
    cond_name = cfg.CONDITIONS(ci).name;
    cond_rates = trial_rates(condition_masks(:, ci), :);

    for ei = 1:n_epochs
        ep = cfg.EPOCHS(ei);
        if ~isKey(epoch_idx_by_name, ep.baseline)
            continue;
        end
        baseline_idx = epoch_idx_by_name(ep.baseline);
        if baseline_idx == ei
            continue;
        end

        field_name = comparison_field_name('baseline', cond_name, ep.name, ep.baseline);
        comparisons.(field_name) = run_statistical_test( ...
            cfg.statistics.baseline_test, cond_rates(:, ei), cond_rates(:, baseline_idx), ...
            struct('comparison_type', 'baseline', 'condition', cond_name, ...
            'epoch', ep.name, 'baseline_epoch', ep.baseline));
    end
end

for ei = 1:n_epochs
    ep = cfg.EPOCHS(ei);
    for ci_a = 1:n_conditions
        for ci_b = (ci_a + 1):n_conditions
            rates_a = trial_rates(condition_masks(:, ci_a), ei);
            rates_b = trial_rates(condition_masks(:, ci_b), ei);

            field_name = comparison_field_name('condition', ep.name, ...
                cfg.CONDITIONS(ci_a).name, cfg.CONDITIONS(ci_b).name);
            comparisons.(field_name) = run_statistical_test( ...
                cfg.statistics.condition_test, rates_a, rates_b, ...
                struct('comparison_type', 'condition', 'epoch', ep.name, ...
                'condition_a', cfg.CONDITIONS(ci_a).name, ...
                'condition_b', cfg.CONDITIONS(ci_b).name));
        end
    end
end
end

function tf = trial_matches_condition(trial_row, cond_def)
tf = true;
param_names = fieldnames(cond_def.parameters);
for pi = 1:numel(param_names)
    pn = param_names{pi};
    if ~isfield(trial_row, pn)
        tf = false;
        return;
    end
    target_vals = cond_def.parameters.(pn);
    tf = tf && ismember(trial_row.(pn), target_vals(:));
end
end

function result = run_statistical_test(test_name, sample_a, sample_b, meta)
result = meta;
result.test = test_name;
result.n_a = 0;
result.n_b = 0;
result.p = NaN;
result.tstat = NaN;
result.df = NaN;
result.mean_a = NaN;
result.mean_b = NaN;

sample_a = sample_a(:);
sample_b = sample_b(:);

switch lower(char(test_name))
    case 'paired_ttest'
        valid = isfinite(sample_a) & isfinite(sample_b);
        sample_a = sample_a(valid);
        sample_b = sample_b(valid);
        if numel(sample_a) < 2
            result.n_a = numel(sample_a);
            result.n_b = numel(sample_b);
            return;
        end
        [~, p, ~, stats] = ttest(sample_a, sample_b);
        result.n_a = numel(sample_a);
        result.n_b = numel(sample_b);
        result.p = p;
        result.tstat = stats.tstat;
        result.df = stats.df;
        result.mean_a = mean(sample_a);
        result.mean_b = mean(sample_b);

    case 'unpaired_ttest'
        sample_a = sample_a(isfinite(sample_a));
        sample_b = sample_b(isfinite(sample_b));
        if numel(sample_a) < 2 || numel(sample_b) < 2
            result.n_a = numel(sample_a);
            result.n_b = numel(sample_b);
            return;
        end
        [~, p, ~, stats] = ttest2(sample_a, sample_b);
        result.n_a = numel(sample_a);
        result.n_b = numel(sample_b);
        result.p = p;
        result.tstat = stats.tstat;
        result.df = stats.df;
        result.mean_a = mean(sample_a);
        result.mean_b = mean(sample_b);

    otherwise
        error('epp_compute_statistics:UnknownTest', 'Unknown statistical test: %s', test_name);
end
end

function field_name = comparison_field_name(comparison_type, varargin)
switch comparison_type
    case 'baseline'
        field_name = sprintf('%s_%s_vs_baseline_%s', varargin{:});
    case 'condition'
        field_name = sprintf('%s_%s_vs_%s', varargin{:});
    otherwise
        error('epp_compute_statistics:UnknownComparisonType', ...
            'Unknown comparison type: %s', comparison_type);
end
field_name = matlab.lang.makeValidName(field_name);
end

function statistics_table = unit_statistics_to_table(unit_statistics)
row_cells = {};
for u = 1:numel(unit_statistics)
    comparison_names = fieldnames(unit_statistics(u).comparisons);
    for c = 1:numel(comparison_names)
        cmp_name = comparison_names{c};
        cmp = unit_statistics(u).comparisons.(cmp_name);
        row_cells{end+1} = comparison_to_table_row(unit_statistics(u).Neuron_ID, cmp_name, cmp); %#ok<AGROW>
    end
end

if isempty(row_cells)
    statistics_table = empty_statistics_table();
else
    statistics_table = struct2table(vertcat(row_cells{:}));
end
end

function row = comparison_to_table_row(neuron_id, comparison_name, cmp)
row = struct( ...
    'Neuron_ID', neuron_id, ...
    'comparison', comparison_name, ...
    'comparison_type', '', ...
    'test', '', ...
    'epoch', '', ...
    'condition', '', ...
    'condition_a', '', ...
    'condition_b', '', ...
    'baseline_epoch', '', ...
    'n_a', 0, ...
    'n_b', 0, ...
    'p', NaN, ...
    'tstat', NaN, ...
    'df', NaN, ...
    'mean_a', NaN, ...
    'mean_b', NaN);

row.comparison_type = cmp.comparison_type;
row.test = cmp.test;
row.epoch = cmp.epoch;
row.n_a = cmp.n_a;
row.n_b = cmp.n_b;
row.p = cmp.p;
row.tstat = cmp.tstat;
row.df = cmp.df;
row.mean_a = cmp.mean_a;
row.mean_b = cmp.mean_b;

switch row.comparison_type
    case 'baseline'
        row.condition = cmp.condition;
        row.condition_a = row.condition;
        row.baseline_epoch = cmp.baseline_epoch;
    case 'condition'
        row.condition_a = cmp.condition_a;
        row.condition_b = cmp.condition_b;
end
end

function statistics_table = empty_statistics_table()
statistics_table = table( ...
    strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    zeros(0, 1), zeros(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), ...
    'VariableNames', { ...
    'Neuron_ID', 'comparison', 'comparison_type', 'test', 'epoch', ...
    'condition', 'condition_a', 'condition_b', 'baseline_epoch', ...
    'n_a', 'n_b', 'p', 'tstat', 'df', 'mean_a', 'mean_b'});
end

function blocks_out = epp_load_saved_processed_blocks(output_root)
blocks_out = [];
files = dir(fullfile(output_root, '*_Block-*.mat'));
for i = 1:numel(files)
    file_path = fullfile(output_root, files(i).name);
    payload = load(file_path, 'block_payload');
    blocks_out = [blocks_out payload.block_payload]; %#ok<AGROW>
end

if ~isempty(blocks_out)
    sort_key = [[blocks_out.session].' [blocks_out.block].'];
    [~, order_idx] = sortrows(sort_key, [1 2]);
    blocks_out = blocks_out(order_idx);
end
end
