function cfg = epp_general_settings(project, cfg)
%EPP_GENERAL_SETTINGS  Default paths and analysis parameters.
%
%   cfg = epp_general_settings(project, cfg)
%
%   Inputs:
%     project - char/string. Project folder name under settings/<project>/.
%     cfg     - struct, optional. Must already contain cfg.project_version
%               (set by epp_load_cfg). Missing cfg starts from struct().
%
%   Outputs:
%     cfg - same struct with cfg.roots.* (settings, sorting_tables,
%           ephys_tanks, behavior, project_version, processed_trials, raster,
%           psth, population_psth, statistics), cfg.sorting_sheet,
%           cfg.spike_threshold, cfg.unit_info_file, cfg.monkeys, cfg.datasets,
%           cfg.raster_bin_size_s, cfg.raster_smooth_samples,
%           cfg.WINDOWS, cfg.EPOCHS (cell table),
%           cfg.CONDITIONS, cfg.psth, cfg.movements, cfg.min_edge_rate_hz,
%           and cfg.statistics filled in.
%           Existing fields are overwritten. Project/version scripts may
%           override them afterward.
%
% Minimal general settings for epp_initiation.
%
% This only defines paths required by the unit initializer.

if nargin < 2 || isempty(cfg)
    cfg = struct;
end

if nargin < 1 || isempty(project)
    error('Provide project name.');
end

% Path roots. Project settings live in this repo under settings/<project>/.
repo_root = fileparts(mfilename('fullpath'));
cfg.roots.settings = [fullfile(repo_root, 'settings', project) filesep];
cfg.roots.sorting_tables = 'Y:\Data\Sorting_tables';
cfg.roots.ephys_tanks = 'Y:\Data\TDTtanks';
cfg.roots.behavior = 'Y:\Data';
cfg.roots.project_version = fullfile('Y:\Projects', project, cfg.project_version);
cfg.roots.processed_trials = fullfile(cfg.roots.project_version, 'behavior');
cfg.roots.raster = fullfile(cfg.roots.project_version, 'unit_rasters');
cfg.roots.psth = fullfile(cfg.roots.project_version, 'unit_psth');
cfg.roots.population_psth = fullfile(cfg.roots.project_version, 'population_psth');
cfg.roots.statistics = fullfile(cfg.roots.project_version, 'statistics');

% Sorting workbook sheet, spike-file threshold token, and unit_info.mat path.
% Spike files are dataspikes_rb<FFF>_ch<CCC>_<spike_threshold>.mat.
cfg.sorting_sheet = 'automatic_sorting';
cfg.spike_threshold = 'negthr';
cfg.unit_info_file = fullfile(cfg.roots.project_version, 'unit_info.mat');

% Default pipeline behavior (project/version settings may override).
cfg.monkeys = {};
cfg.datasets = [];

% Drop leading and trailing trials below this rate (spikes/s). 0 disables.
% Rate is that unit's spikes from one trial timestamp to the next.
cfg.min_edge_rate_hz = 0;

% Default raster settings used by epp_initiation.
cfg.raster_bin_size_s = 0.001;
% Boxcar width in raster bins after interpolation. 1 leaves the trace unchanged.
% Which traces are downsampled first is behavior_to_downsample in epp_build_rasters.
cfg.raster_smooth_samples = 1;
cfg.WINDOWS = struct( ...
    'name', {'Default'}, ...
    'align_state', {2}, ...
    't_start_s', {-0.5}, ...
    't_end_s', {1.0});


% EPOCHS table: name | state | start | end
% Converted to struct array in epp_load_cfg.
cfg.EPOCHS = { ...
    'INI',   2,  -0.4,  -0.1; ...
    'Facq',  3,  -0.4,  -0.1; ...
    'Fhol',  6,  -0.3,   0; ...
    'Cue',   6,   0.06,  0.12; ...
    'Del',   4,  -0.3,   0; ...
    'PreS',  60, -0.22, -0.02; ...
    'PeriS', 60, -0.02,  0.08; ...
    'PostS', 61,  0.05,  0.2; ...
    'PreR',  62, -0.4,  -0.1; ...
    'PeriR', 62, -0.05,  0.25; ...
    'PostR', 63,  0.1,   0.4; ...
    'Thol',  20, -0.3,   0; ...
    };

cfg.CONDITIONS = struct( ...
    'name', {'Ipsi', 'Contra'}, ...
    'color', {[200 55 12], [100 55 22]}, ...
    'parameters', num2cell(struct('choice', {0, 0}, 'hemifield', {-1, 1})));


% Saccade/reach detection. Passed to MP_add_saccades_and_reaches as cfg.movements.
% saccade_definition: 'closest_big_enough' (default), 'biggest_close_enough',
%   'last_in_state', 'first_in_state', 'first_big_enough',
%   'first_end_in_target', 'first_start_in_target'
% reach_definition: 'first_touch_in_target' (default), 'first_touch',
%   'state_change_position'
cfg.movements.interpolated_SR = 1000;
cfg.movements.downsampling = 1;
cfg.movements.smoothing_samples = 12;
cfg.movements.sac_ini_t = 200;       % saccade onset velocity threshold
cfg.movements.sac_end_t = 50;        % saccade offset velocity threshold
cfg.movements.sac_min_dur = 0.03;
cfg.movements.sac_min_amp = 2;
cfg.movements.sac_max_off = 50;
cfg.movements.nsacc_max = 20;
cfg.movements.saccade_definition = 'closest_big_enough';
cfg.movements.reach_definition = 'first_touch_in_target';
cfg.movements.correct_offset = 1;
cfg.movements.enable_cleaning = 1;

% Default PSTH plotting settings.
cfg.psth.bin_size_s = 0.01;
cfg.psth.smoothing_kernel = 'gaussian'; % 'gaussian' or 'box'
cfg.psth.smoothing_width_s = 0.03;

% Statistics (cell tables; merged to cfg.statistics.comparisons in epp_load_cfg).
%
% within_epoch:  name | epoch | conditions | baseline_conditions
% across_epochs: name | epoch | conditions | baseline_epoch
cfg.statistics.within_epoch = { ...
    'PreSpace',  'PreR', [1 2], [3 4]; ...
    'PeriSpace', 'PeriR', [1 2], [3 4]; ...
    };

cfg.statistics.across_epochs = { ...
    'Cue_vs_INI',  'Cue',   [1 2], 'INI'; ...
    'PeriR_vs_INI','PeriR', [1 2], 'INI'; ...
    };

cfg.statistics.within_epoch_test = 'unpaired_ttest';
cfg.statistics.across_epochs_test = 'paired_ttest';

end
