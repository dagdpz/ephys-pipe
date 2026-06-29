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

if isfield(cfg, 'EPOCHS')
    cfg.EPOCHS = epochs_table_to_struct(cfg.EPOCHS);
end

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
% Settings files define rows as:
%   name | state | start | end | baseline
% Empty rows (blank name) are skipped so lines can be added/removed easily.


epochs = repmat(struct( ...
    'name', '', 'align_state', [], 't_start_s', [], 't_end_s', [], 'baseline', ''), 0, 1);

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
    ep.baseline = epochs_table{r, 5};
    epochs = [epochs; ep]; %#ok<AGROW>
end
end

