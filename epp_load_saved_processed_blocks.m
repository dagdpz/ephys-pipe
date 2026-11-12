function blocks_out = epp_load_saved_processed_blocks(output_root)
%EPP_LOAD_SAVED_PROCESSED_BLOCKS  Load stage-1 block payloads.
%
%   blocks_out = epp_load_saved_processed_blocks(output_root)
%
%   Reads every *_Block-*.mat in output_root (variable block_payload) and
%   returns them sorted by session, then block. epp_initiation loads them
%   once and passes them into epp_build_rasters.

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
