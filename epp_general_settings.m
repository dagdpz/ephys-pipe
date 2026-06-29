function cfg = epp_general_settings(project, cfg)
% Minimal general settings for epp_initiation.
%
% This only defines paths required by the unit initializer.

if nargin < 2 || isempty(cfg)
    cfg = struct;
end

if nargin < 1 || isempty(project)
    error('Provide project name.');
end

settings_root = fileparts(mfilename('fullpath'));
fileseps=strfind(settings_root,filesep);
settings_root =[settings_root(1:fileseps(end)) 'Settings'];

cfg.roots.settings = [settings_root filesep project filesep 'epp' filesep ];

% Path roots.
cfg.roots.sorting_tables = 'Y:\Data\Sorting_tables';
cfg.roots.ephys_tanks = 'Y:\Data\TDTtanks';
cfg.roots.behavior = 'Y:\Data';
cfg.roots.project_version = fullfile('Y:\Projects', project, cfg.project_version);
cfg.roots.processed_trials = fullfile(cfg.roots.project_version, 'behavior');
cfg.roots.raster = fullfile(cfg.roots.project_version, 'unit_rasters');
cfg.roots.psth = fullfile(cfg.roots.project_version, 'unit_psth');
cfg.roots.population_psth = fullfile(cfg.roots.project_version, 'population_psth');
cfg.roots.statistics = fullfile(cfg.roots.project_version, 'statistics');


% Default pipeline behavior (project/version settings may override).
cfg.monkeys = {};
cfg.datasets = [];

% Default raster settings used by epp_initiation.
cfg.raster_bin_size_s = 0.001;
cfg.WINDOWS = struct( ...
    'name', {'Default'}, ...
    'align_state', {2}, ...
    't_start_s', {-0.5}, ...
    't_end_s', {1.0});



% EPOCHS table: one row per epoch — name | state | start | end | baseline
% Converted to struct array in epp_load_cfg (fields align with cfg.WINDOWS + baseline).
cfg.EPOCHS = { ...
    'INI',   2,  -0.4,  -0.1,  'INI'; ...
    'Facq',  3,  -0.4,  -0.1,  'INI'; ...
    'Fhol',  6,  -0.3,   0,    'INI'; ...
    'Cue',   6,   0.06,  0.12,  'INI'; ...
    'Del',   4,  -0.3,   0,    'INI'; ...
    'PreS',  60, -0.22, -0.02,  'INI'; ...
    'PeriS', 60, -0.02,  0.08,  'INI'; ...
    'PostS', 61,  0.05,  0.2,   'INI'; ...
    'PreR',  62, -0.4,  -0.1,   'INI'; ...
    'PeriR', 62, -0.05,  0.25,  'INI'; ...
    'PostR', 63,  0.1,   0.4,   'INI'; ...
    'Thol',  20, -0.3,   0,    'INI'; ...
    };


% Default PSTH plotting settings.
cfg.psth.bin_size_s = 0.01;
cfg.psth.smoothing_kernel = 'gaussian'; % 'gaussian' or 'box'
cfg.psth.smoothing_width_s = 0.03;

% Default epoch statistics settings.
cfg.statistics.baseline_test = 'paired_ttest';
cfg.statistics.condition_test = 'unpaired_ttest';

end
