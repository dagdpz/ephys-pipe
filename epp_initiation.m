function epp_initiation(project, version)
%EPP_INITIATION  Run the ephys-pipe stages in order.
%
%   epp_initiation(project, version)
%
%   Inputs:
%     project - char/string. Settings folder name under settings/<project>/.
%     version - char/string. Selects settings/<project>/<version>.m and the output
%               folder Y:\Projects\<project>\<version>\.
%
%   Outputs:
%     none. Calls epp_load_cfg and epp_prepare_blocks, loads unit_info once,
%     then for each unit calls epp_build_rasters, epp_compute_statistics, and
%     epp_plot_psth. Writes the statistics table and runs
%     epp_plot_population_psth after that loop.
%
% Orchestrate the ephys-alignment pipeline stages.

repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);
addpath(genpath(fullfile(repo_root, 'external')));

cfg = epp_load_cfg(project, version);
epp_prepare_blocks(cfg);

loaded_units = load(cfg.unit_info_file, 'unit_info');
unit_info = loaded_units.unit_info;

processed_blocks = [];
files = dir(fullfile(output_root, '*_Block-*.mat'));
for i = 1:numel(files)
    file_path = fullfile(output_root, files(i).name);
    payload = load(file_path, 'block_payload');
    processed_blocks = [processed_blocks payload.block_payload]; %#ok<AGROW>
end

unit_statistics = repmat(struct('Neuron_ID', '', 'comparisons', struct()), numel(unit_info), 1);

for u = 1:numel(unit_info)
    epp_build_rasters(cfg, unit_info(u), processed_blocks);
    unit_statistics(u) = epp_compute_statistics(cfg, unit_info(u));
    epp_plot_psth(cfg, unit_info(u));
end

epp_save_unit_statistics(cfg, unit_statistics);
epp_plot_population_psth(cfg);
end



