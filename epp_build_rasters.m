function epp_build_rasters(cfg)
% Stage 2: build per-unit rasters from processed blocks and spikes.

load(fullfile(cfg.roots.processed_trials, 'unit_info.mat'), 'unit_info');
processed_blocks = epp_load_saved_processed_blocks(cfg.roots.processed_trials);

for u = 1:numel(unit_info)
    session = unit_info(u).Session;
    Neuron_ID = unit_info(u).Neuron_ID;
    blocks = unit_info(u).Blocks;
    monkey = unit_info(u).monkey;
    channel = unit_info(u).Channel;
    filenumber = unit_info(u).Filenumber;
    sortcode = unit_info(u).Unit;
    threshold = 'negthr';

    day_token = sprintf('%08d', session);
    session_folder = fullfile(cfg.roots.ephys_tanks, [monkey '_phys'], day_token);
    spike_file = [session_folder filesep 'dataspikes_rb', sprintf('%03d',filenumber), '_ch' sprintf('%03d',channel) '_' threshold '.mat'];
    spike_data = load(spike_file, 'cluster_class','par');
    cc = spike_data.cluster_class;
    par = spike_data.par;

    segment_ends = par.segmentends(:) ./ par.sr;
    segment_starts = [0; segment_ends(1:end-1)];
    spike_times = cc(cc(:,1) == sortcode, 2)/1000;

    unit_combined = struct( ...
        'trial', struct([]), ...
        'state', struct('Value', [], 'timestamps', []), ...
        'sample', struct('x_hnd', [], 'y_hnd', [], 'x_eye', [], 'y_eye', [], 'timestamps', []));

    for bi = 1:numel(blocks)
        b = blocks(bi);
        idx_block = find([processed_blocks.session] == session & [processed_blocks.block] == b & strcmp({processed_blocks.monkey}, monkey), 1, 'first');
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

    unit_combined.trial = epp_enrich_trials_for_unit(unit_combined.trial, unit_info(u));

    window_rasters = cfg.WINDOWS(:);
    for wi = 1:numel(window_rasters)
        w = window_rasters(wi);
        bin_edges = w.t_start_s:cfg.raster_bin_size_s:w.t_end_s;
        bin_size_s = diff(bin_edges(1:2));
        event_times = unit_combined.state.timestamps(unit_combined.state.Value == w.align_state);
        event_times = event_times(isfinite(event_times));

        raster_counts = zeros(numel(event_times), numel(bin_edges)-1);
        for k = 1:numel(event_times)
            rel_spike_times = spike_times - event_times(k);
            raster_counts(k,:) = histcounts(rel_spike_times, bin_edges);
        end

        event_trials = struct([]);
        if ~isempty(event_times) && ~isempty(unit_combined.trial)
            trial_timestamps = [unit_combined.trial.timestamp].';
            event_trials = repmat(unit_combined.trial(1), numel(event_times), 1);
            for k = 1:numel(event_times)
                idx = find(trial_timestamps <= event_times(k), 1, 'last');
                if isempty(idx), idx = 1; end
                event_trials(k,1) = unit_combined.trial(idx);
            end
        end

        window_rasters(wi).bin_edges = bin_edges(:);
        window_rasters(wi).bin_centers = bin_edges(1:end-1).' + bin_size_s/2;
        window_rasters(wi).bin_size_s = bin_size_s;
        window_rasters(wi).event_times = event_times(:);
        window_rasters(wi).raster = raster_counts;
        window_rasters(wi).rate_hz = raster_counts ./ bin_size_s;
        window_rasters(wi).trial = event_trials;
    end

    raster_data = struct();
    raster_data.Neuron_ID = Neuron_ID;
    raster_data.windows = window_rasters;
    raster_file = sprintf('%s_raster.mat', Neuron_ID);
    save(fullfile(cfg.roots.raster, raster_file), 'raster_data');
end
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