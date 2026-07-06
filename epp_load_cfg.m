function cfg = epp_load_cfg(project, version)

if nargin < 1 || isempty(project)
    error('Provide project name.');
end
if nargin < 2 || isempty(version)
    error('Provide version name.');
end

cfg = struct;
cfg.project_version = version;
cfg = epp_general_settings(project, cfg);

project_settings = [cfg.roots.settings 'epp_project_settings.m'];
if exist(project_settings, 'file')
    run(project_settings);
else
    warning('Missing project settings file: %s', project_settings);
end

if ~isempty(cfg.project_version)
    version_settings = [cfg.roots.settings sprintf('epp_%s_settings.m', cfg.project_version)];
    if exist(version_settings, 'file')
        run(version_settings);
    else
        warning('Missing version settings file: %s', version_settings);
    end
end

cfg.EPOCHS = epochs_table_to_struct(cfg.EPOCHS);
cfg.statistics.comparisons = build_statistics_comparisons(cfg.statistics);

root_fields = fieldnames(cfg.roots);
for i = 1:numel(root_fields)
    fn = root_fields{i};
    if ~isempty(cfg.roots.(fn)) && ~isfolder(cfg.roots.(fn))
        mkdir(cfg.roots.(fn));
    end
end

end

function epochs = epochs_table_to_struct(epochs_table)
% Convert cfg.EPOCHS cell table to struct array (same field style as cfg.WINDOWS).
%
% Settings rows: name | state | start | end
% Empty rows (blank name) are skipped.

if isstruct(epochs_table)
    epochs = epochs_table(:);
    return;
end

epochs = repmat(struct('name', '', 'align_state', [], 't_start_s', [], 't_end_s', []), 0, 1);

if isempty(epochs_table) || ~iscell(epochs_table)
    return;
end

if size(epochs_table, 2) ~= 4
    error('epp_load_cfg:InvalidEpochsColumns', ...
        'EPOCHS table must have 4 columns: name, state, start, end.');
end

for r = 1:size(epochs_table, 1)
    name_val = epochs_table{r, 1};
    if isempty(name_val)
        continue;
    end

    ep = struct();
    ep.name = name_val;
    ep.align_state = epochs_table{r, 2};
    ep.t_start_s = epochs_table{r, 3};
    ep.t_end_s = epochs_table{r, 4};
    epochs = [epochs; ep]; %#ok<AGROW>
end
end

function comparisons = build_statistics_comparisons(statistics)
% Merge within_epoch and across_epochs tables into one comparison struct array.

comparisons = repmat(empty_comparison_struct(), 0, 1);

for r = 1:size(statistics.within_epoch, 1)
    if isempty(statistics.within_epoch{r, 1})
        continue;
    end
    cmp = comparison_row_to_struct(statistics.within_epoch(r, :), ...
        'within_epoch', statistics.within_epoch_test);
    comparisons = [comparisons; cmp]; %#ok<AGROW>
end

for r = 1:size(statistics.across_epochs, 1)
    if isempty(statistics.across_epochs{r, 1})
        continue;
    end
    cmp = comparison_row_to_struct(statistics.across_epochs(r, :), ...
        'across_epochs', statistics.across_epochs_test);
    comparisons = [comparisons; cmp]; %#ok<AGROW>
end
end

function cmp = comparison_row_to_struct(row, comparison_scope, test_name)
if numel(row) ~= 4
    error('epp_load_cfg:InvalidComparisonColumns', ...
        'Statistics comparison rows need 4 columns (see within_epoch / across_epochs format).');
end

cmp = empty_comparison_struct();
cmp.comparison_scope = comparison_scope;
cmp.name = row{1};
cmp.epoch = row{2};
cmp.conditions = condition_index_vector(row{3});
cmp.test = test_name;

switch comparison_scope
    case 'within_epoch'
        cmp.baseline_epoch = row{2};
        cmp.baseline_conditions = condition_index_vector(row{4});
    case 'across_epochs'
        cmp.baseline_epoch = row{4};
        cmp.baseline_conditions = cmp.conditions;
    otherwise
        error('epp_load_cfg:UnknownComparisonScope', ...
            'Unknown comparison scope: %s', comparison_scope);
end
end

function cmp = empty_comparison_struct()
cmp = struct( ...
    'comparison_scope', '', ...
    'name', '', ...
    'epoch', '', ...
    'conditions', [], ...
    'baseline_epoch', '', ...
    'baseline_conditions', [], ...
    'test', '');
end

function idx = condition_index_vector(value)
idx = double(value(:).');
end

