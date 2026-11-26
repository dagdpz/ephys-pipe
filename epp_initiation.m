function epp_initiation(project, version)
%EPP_INITIATION  Run the ephys-pipe stages in order.
%
%   epp_initiation(project, version)
%
%   Inputs:
%     project - char/string. Settings folder under settings/<project>/.
%     version - char/string. Loads settings/<project>/<version>.m and writes
%               under Y:\Projects\<project>\<version>\.
%
%   Pipeline (each gated by cfg.todo.*):
%     1. unit_info           sorting tables -> cfg.unit_info_file
%     2. session_range filter  keep units whose Session is in
%                            cfg.session_range.<Monkey> = [start end]
%                            (always applied when ranges are set)
%     3. prepare_blocks      epp_prepare_blocks -> processed trial mats
%     4. build_rasters       epp_load_unit_combined + epp_build_rasters
%     5. plot_psth           epp_plot_psth (needs accumulated_blocks)
%     6. compute_statistics  collect this run's units, then merge into
%                            unit_statistics.mat by Neuron_ID in one write
%     7. plot_population_psth epp_plot_population_psth from rasters
%
%   Outputs:
%     none. Side effects are the stage outputs under cfg.roots.*.

repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);
addpath(genpath(fullfile(repo_root, 'external')));

cfg = epp_load_cfg(project, version);

% --- unit_info ---
if cfg.todo.unit_info
    unit_info = struct([]);
    for m = 1:numel(cfg.monkeys)
        monkey = cfg.monkeys{m};
        monkey_units = epp_sorting_table_to_units(monkey, cfg);
        unit_info = [unit_info; monkey_units]; %#ok<AGROW>
    end
    save(cfg.unit_info_file, 'unit_info');
else
    loaded_units = load(cfg.unit_info_file, 'unit_info');
    unit_info = loaded_units.unit_info;
end

% --- reduce to session range ---
monkeys=intersect(cfg.monkeys,fieldnames(cfg.session_range));
for m = 1:numel(monkeys)
    rng = cfg.session_range.(monkeys{m});
    sess = [unit_info.Session]';
    to_skip = (sess < rng(1) | sess > rng(2)) & ismember({unit_info.Monkey}',monkeys{m});
    unit_info(to_skip) = [];
end

% --- prepare_blocks ---
if cfg.todo.prepare_blocks
    epp_prepare_blocks(cfg, unit_info);
end

if cfg.todo.build_rasters || cfg.todo.plot_psth
    accumulated_blocks = [];
    files = dir(fullfile(cfg.roots.processed_trials, '*_Block-*.mat'));
    for i = 1:numel(files)
        file_path = fullfile(cfg.roots.processed_trials, files(i).name);
        payload = load(file_path, 'block_payload');
        accumulated_blocks = [accumulated_blocks payload.block_payload]; %#ok<AGROW>
    end
end

% Stats for units processed this run; merged into the combined file after the loop.

stats_file = fullfile(cfg.roots.statistics, 'unit_statistics.mat');
if isfile(stats_file)
    unit_statistics = load(stats_file).unit_statistics;
else
    unit_statistics = repmat(struct('Neuron_ID', '', 'comparisons', struct()), 0, 1);
end

% ---  per unit processing ---
for u = 1:numel(unit_info)
    % --- build_rasters ---
    if cfg.todo.build_rasters
        unit = epp_load_unit_combined(unit_info(u), accumulated_blocks, cfg);
        epp_build_rasters(cfg, unit);
    end

    % --- plot_psth ---
    if cfg.todo.plot_psth
        epp_plot_psth(cfg, unit_info(u), accumulated_blocks);
    end

    % --- compute_statistics ---
    if cfg.todo.compute_statistics
        stats = epp_compute_statistics(cfg, unit_info(u)); 
        hit = ismember({unit_statistics.Neuron_ID},unit_info(u).Neuron_ID);
        idx = min([find(hit) numel(unit_statistics)+1]);
        unit_statistics(idx)=stats;
    end
end
% --- save statistics ---
epp_save_unit_statistics(cfg, unit_statistics);


% --- population PSTH ---
if cfg.todo.plot_population_psth
    epp_plot_population_psth(cfg);
end
end
