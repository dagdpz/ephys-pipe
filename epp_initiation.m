function epp_initiation(project, version)
%EPP_INITIATION  Run the four ephys-pipe stages in order.
%
%   epp_initiation(project, version)
%
%   Inputs:
%     project - char/string. Settings folder name under ../Settings/<project>/epp/.
%     version - char/string. Selects epp_<version>_settings.m and the output
%               folder Y:\Projects\<project>\<version>\.
%
%   Outputs:
%     none. Calls epp_load_cfg, then epp_prepare_blocks, epp_build_rasters,
%     epp_plot_psth, and epp_plot_population_psth. Files written by those
%     stages are the pipeline outputs.
%
% Orchestrate the ephys-alignment pipeline stages.

cfg = epp_load_cfg(project, version);
epp_prepare_blocks(cfg);
epp_build_rasters(cfg);
epp_plot_psth(cfg);
epp_plot_population_psth(cfg);
end



