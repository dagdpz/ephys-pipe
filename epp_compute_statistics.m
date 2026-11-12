function unit_statistic = epp_compute_statistics(cfg, unit_row)
%EPP_COMPUTE_STATISTICS  Epoch firing-rate comparisons for one unit.
%
%   unit_statistic = epp_compute_statistics(cfg, unit_row)
%
%   Inputs:
%     cfg - struct from epp_load_cfg. Reads
%           cfg.roots.raster/<Neuron_ID>_raster.mat (field epochs, written
%           by epp_build_rasters). Uses cfg.CONDITIONS and
%           cfg.statistics.comparisons. Does not read spike files.
%     unit_row - one element of unit_info.
%
%   Outputs:
%     unit_statistic - struct with .Neuron_ID and .comparisons.<name>
%           (n_a, n_b, p, tstat, df, mean_a, mean_b, comparison metadata).
%           epp_save_unit_statistics writes the table for every unit.
%
%   Each epochs row is one align event. Unpaired tests use those events
%   after the condition mask. Paired tests match events on trial_n and keep
%   a trial only when both epochs have the same number of events from it.

unit_statistic = struct('Neuron_ID', unit_row.Neuron_ID, 'comparisons', struct());
raster_file = fullfile(cfg.roots.raster, sprintf('%s_raster.mat', unit_row.Neuron_ID));
if ~isfile(raster_file)
    error('epp_compute_statistics:MissingRaster', ...
        'No raster for %s. Run epp_build_rasters first.', unit_row.Neuron_ID);
end
loaded = load(raster_file, 'raster_data');
if ~isfield(loaded.raster_data, 'epochs')
    error('epp_compute_statistics:NoEpochRates', ...
        'Raster %s has no epochs. Re-run epp_build_rasters.', unit_row.Neuron_ID);
end

unit_statistic.comparisons = build_unit_comparisons(loaded.raster_data.epochs, cfg);
end

function comparisons = build_unit_comparisons(epochs, cfg)
comparisons = struct();
if isempty(epochs)
    epoch_idx_by_name = containers.Map('KeyType', 'char', 'ValueType', 'double');
else
    epoch_idx_by_name = containers.Map({epochs.name}, num2cell(1:numel(epochs)));
end

for i = 1:numel(cfg.statistics.comparisons)
    cmp_def = cfg.statistics.comparisons(i);
    fr_a = epochs(resolve_epoch_index(epoch_idx_by_name, cmp_def.epoch));
    fr_b = epochs(resolve_epoch_index(epoch_idx_by_name, cmp_def.baseline_epoch));
    [rates_a, rates_b] = comparison_rates(fr_a, fr_b, cmp_def, cfg.CONDITIONS);

    meta = struct( ...
        'comparison_scope', cmp_def.comparison_scope, ...
        'name', cmp_def.name, ...
        'epoch', cmp_def.epoch, ...
        'conditions', cmp_def.conditions, ...
        'baseline_epoch', cmp_def.baseline_epoch, ...
        'baseline_conditions', cmp_def.baseline_conditions);

    field_name = comparison_field_name(cmp_def.name);
    comparisons.(field_name) = run_statistical_test(cmp_def.test, rates_a, rates_b, meta);
end
end

function [rates_a, rates_b] = comparison_rates(fr_a, fr_b, cmp_def, conditions)
idx_a = condition_index_for_events(fr_a, conditions);
idx_b = condition_index_for_events(fr_b, conditions);
mask_a = fr_a.trial_n(:) > 0 & ismember(idx_a, cmp_def.conditions(:));
mask_b = fr_b.trial_n(:) > 0 & ismember(idx_b, cmp_def.baseline_conditions(:));

if strcmpi(cmp_def.test, 'paired_ttest')
    [rates_a, rates_b] = pair_epoch_rates_by_trial(fr_a, fr_b, mask_a, mask_b);
else
    rates_a = fr_a.rate_hz(mask_a);
    rates_b = fr_b.rate_hz(mask_b);
end
end

function idx = condition_index_for_events(fr, conditions)
idx = zeros(numel(fr.rate_hz), 1);
valid = fr.trial_n(:) > 0;
if ~any(valid) || isempty(fr.trial)
    return
end
idx(valid) = assign_trial_condition_indices(fr.trial(valid), conditions);
end

function [rates_a, rates_b] = pair_epoch_rates_by_trial(fr_a, fr_b, mask_a, mask_b)
% A trial contributes only when both epochs have the same event count for it.
% Those events are paired in event_time order.
idx_a = find(mask_a);
idx_b = find(mask_b);
shared = intersect(fr_a.trial_n(idx_a), fr_b.trial_n(idx_b), 'stable');
rates_a = zeros(0, 1);
rates_b = zeros(0, 1);
for s = 1:numel(shared)
    trial_id = shared(s);
    sel_a = idx_a(fr_a.trial_n(idx_a) == trial_id);
    sel_b = idx_b(fr_b.trial_n(idx_b) == trial_id);
    if numel(sel_a) ~= numel(sel_b)
        continue
    end
    [~, order_a] = sort(fr_a.event_time(sel_a));
    [~, order_b] = sort(fr_b.event_time(sel_b));
    rates_a = [rates_a; fr_a.rate_hz(sel_a(order_a))]; %#ok<AGROW>
    rates_b = [rates_b; fr_b.rate_hz(sel_b(order_b))]; %#ok<AGROW>
end
end

function trial_condition_idx = assign_trial_condition_indices(trials, conditions)
trial_condition_idx = zeros(numel(trials), 1);
for ti = 1:numel(trials)
    for ci = 1:numel(conditions)
        if trial_matches_condition(trials(ti), conditions(ci))
            trial_condition_idx(ti) = ci;
            break;
        end
    end
end
end

function epoch_idx = resolve_epoch_index(epoch_idx_by_name, epoch_name)
if ~isKey(epoch_idx_by_name, epoch_name)
    error('epp_compute_statistics:UnknownEpoch', 'Unknown epoch: %s', epoch_name);
end
epoch_idx = epoch_idx_by_name(epoch_name);
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

switch lower(test_name)
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

function field_name = comparison_field_name(comparison_name)
field_name = matlab.lang.makeValidName(comparison_name);
end
