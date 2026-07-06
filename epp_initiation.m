function epp_initiation(project, version)
% Orchestrate the ephys-alignment pipeline stages.

cfg = epp_load_cfg(project, version);
epp_prepare_blocks(cfg);
epp_build_rasters(cfg);
epp_plot_psth(cfg);
epp_plot_population_psth(cfg);
end



