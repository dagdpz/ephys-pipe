function out_file = epp_convert_electrode_depths(depth_file, neurons_file, out_file)
%EPP_CONVERT_ELECTRODE_DEPTHS  Old Electrode_depths_*.m to the current format.
%
%   out_file = epp_convert_electrode_depths(depth_file, neurons_file)
%   out_file = epp_convert_electrode_depths(depth_file, neurons_file, out_file)
%
%   Adds perturbation{k} from neurons_file sheet final_sorting, matching
%   Date and Block (the Session column there is a counter, not the date).
%   Within each session, a channel depth within 50 um of that channel's
%   previous depth is rewritten to that previous depth.
%   Comment lines are kept. depth_file is not modified.
%   Default output is <depth_file>_new.m in the same folder.

if nargin < 3 || isempty(out_file)
    [folder, name] = fileparts(depth_file);
    out_file = fullfile(folder, [name '_new.m']);
end
if strcmpi(char(depth_file), char(out_file))
    error('epp_convert_electrode_depths:SameFile', ...
        'Refusing to overwrite %s.', depth_file);
end

same_depth_um = 50;
depth_src = fileread(depth_file);
run(depth_file);
n_rows = numel(Session);
perturbation = perturbation_by_date_block(neurons_file, Session, block);

last_z = containers.Map('KeyType', 'double', 'ValueType', 'double');
current_session = NaN;
for k = 1:n_rows
    if Session{k} ~= current_session
        last_z = containers.Map('KeyType', 'double', 'ValueType', 'double');
        current_session = Session{k};
    end
    ch = channels{k}(:).';
    zz = z{k}(:).';
    if numel(ch) ~= numel(zz)
        error('epp_convert_electrode_depths:ChannelDepthMismatch', ...
            'Session %.0f block %.0f: %d channels and %d depths.', ...
            Session{k}, block{k}, numel(ch), numel(zz));
    end
    for c = 1:numel(ch)
        if isKey(last_z, ch(c)) && abs(zz(c) - last_z(ch(c))) < same_depth_um
            zz(c) = last_z(ch(c));
        end
        last_z(ch(c)) = zz(c);
    end
    z{k} = zz;
end

lines = regexp(depth_src, '\r\n|\n|\r', 'split');
fid = fopen(out_file, 'w');
if fid < 0
    error('epp_convert_electrode_depths:CannotWrite', 'Cannot write %s.', out_file);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
k = 0;
for i = 1:numel(lines)
    line = lines{i};
    code = regexprep(line, '%.*$', '');
    if isempty(regexp(code, 'k\s*=\s*k\s*\+\s*1', 'once'))
        fprintf(fid, '%s\n', line);
        continue
    end
    k = k + 1;
    comment = regexp(line, '%.*', 'match', 'once');
    if isempty(comment)
        comment = '';
    else
        comment = [' ' comment];
    end
    fprintf(fid, ['k=k+1; Session{k}=%.0f; block{k}=%.0f; perturbation{k}=%.0f; ' ...
        'channels{k}=[%s]; z{k}=[%s];%s\n'], ...
        Session{k}, block{k}, perturbation(k), format_vec(channels{k}), format_vec(z{k}), comment);
end
if k ~= n_rows
    error('epp_convert_electrode_depths:RowCount', ...
        'File has %d assignment lines, run() produced %d rows.', k, n_rows);
end
end

function perturbation = perturbation_by_date_block(neurons_file, session, block)
[~, ~, raw] = xlsread(neurons_file, 'final_sorting');
date_col = find(strcmp(raw(1,:), 'Date'), 1);
block_col = find(strcmp(raw(1,:), 'Block'), 1);
pert_col = find(strcmp(raw(1,:), 'Perturbation'), 1);
if isempty(date_col) || isempty(block_col) || isempty(pert_col)
    error('epp_convert_electrode_depths:MissingColumn', ...
        'final_sorting needs Date, Block, and Perturbation.');
end

pert_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
for r = 2:size(raw, 1)
    date_val = cell_scalar(raw{r, date_col});
    block_val = cell_scalar(raw{r, block_col});
    pert_val = cell_scalar(raw{r, pert_col});
    if ~isfinite(date_val) || ~isfinite(block_val) || ~isfinite(pert_val)
        continue
    end
    key = sprintf('%.0f|%.0f', date_val, block_val);
    if isKey(pert_map, key) && pert_map(key) ~= pert_val
        warning('epp_convert_electrode_depths:PerturbationConflict', ...
            'Date %.0f block %.0f has perturbation %.0f and %.0f. Keeping %.0f.', ...
            date_val, block_val, pert_map(key), pert_val, pert_map(key));
        continue
    end
    pert_map(key) = pert_val;
end

perturbation = zeros(1, numel(session));
for k = 1:numel(session)
    key = sprintf('%.0f|%.0f', session{k}, block{k});
    if isKey(pert_map, key)
        perturbation(k) = pert_map(key);
    else
        warning('epp_convert_electrode_depths:NoPerturbation', ...
            'No Perturbation for date %.0f block %.0f. Using 0.', session{k}, block{k});
    end
end
end

function value = cell_scalar(x)
if isnumeric(x) && isscalar(x)
    value = double(x);
elseif ischar(x) || isstring(x)
    value = str2double(x);
else
    value = NaN;
end
end

function text = format_vec(values)
text = strtrim(sprintf('%g ', values(:).'));
end
