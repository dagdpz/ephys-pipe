function epp_ensure_paths()
%EPP_ENSURE_PATHS  Put this repo and external/ on the MATLAB path.
%
%   epp_ensure_paths()
%
%   external/ holds vendored dependencies and their own dependencies:
%     DAGtools/         DAG_get_server_IP, DAG_find_column_index,
%                       DAG_update_cell_table, DAG_update_mastertable_cell
%     shadedErrorBar/   shadedErrorBar, shadedErrorBar_2014, shadedErrorBar_2011
%     TDT/              SEV2mat_working (called by TDTbin2mat_working for streams)

repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);
addpath(genpath(fullfile(repo_root, 'external')));
end
