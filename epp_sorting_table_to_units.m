function units = epp_sorting_table_to_units(monkey, cfg)
% epp_sorting_table_to_units
% Convert sorting-table rows to unit struct array with dataset/unit filtering.

units = struct([]);
filename = [monkey(1:3) '_UltraSort'];
sheetname = 'final_sorting';
foldername = fullfile(cfg.roots.sorting_tables, monkey);

xlsx_matches = dir(fullfile(foldername, [filename '.xls*']));
if isempty(xlsx_matches)
    warning('No sorting table found for %s in %s', filename, foldername);
    return;
end

excel_sheet = fullfile(foldername, xlsx_matches(1).name);
[~, ~, sorting_table] = xlsread(excel_sheet, sheetname);
if isempty(sorting_table) || numel(sorting_table) < 2
    warning('Sorting table empty in %s (sheet %s).', excel_sheet, sheetname);
    return;
end
header_row = sorting_table(1,:);

data_rows = sorting_table(2:end,:);
keep_mask = true(size(data_rows,1),1);

if ~isempty(cfg.datasets)
    idx_set = strcmpi(header_row, 'Dataset');
    set_numeric = cell2mat(data_rows(:, idx_set));
    keep_mask = ismember(set_numeric, cfg.datasets);
end

idx_unit = strcmpi(header_row, 'unit');
unit_values = data_rows(:, idx_unit);
is_zero_unit = cellfun(@(x) (isnumeric(x) || islogical(x)) && isscalar(x) && double(x) == 0, unit_values);
keep_mask = keep_mask & ~is_zero_unit;

filtered_rows = data_rows(keep_mask, :);
if isempty(filtered_rows)
    return;
end

field_names = cell(size(header_row));
for c = 1:numel(header_row)
    col_name = header_row{c};
    if ~(ischar(col_name) || isstring(col_name)) || isempty(col_name)
        col_name = sprintf('column_%d', c);
    end
    field_names{c} = matlab.lang.makeValidName(strtrim(char(col_name)));
end

typed_by_col = cell(1, size(filtered_rows,2));
type_by_col = cell(1, size(filtered_rows,2));
for c = 1:size(filtered_rows,2)
    [type_by_col{c}, typed_by_col{c}] = infer_and_convert_column_type(filtered_rows(:,c));
end

units = repmat(struct(), size(filtered_rows,1), 1);
for r = 1:size(filtered_rows,1)
    for c = 1:size(filtered_rows,2)
        if strcmp(type_by_col(c), 'logical') || strcmp(type_by_col(c), 'numeric')
            value_out = typed_by_col{c}(r);
        else
            value_out = filtered_rows{r,c};
        end
        if strcmpi(field_names{c}, 'Blocks') || strcmpi(field_names{c}, 'Perturbation')
            value_out = parse_blocks_cell(value_out);
        end
        units(r).(field_names{c}) = value_out;
    end
    units(r) = normalize_unit_blocks_and_perturbation(units(r));
end

end

function unit_row = normalize_unit_blocks_and_perturbation(unit_row)
% Perturbation is a numeric vector parallel to Blocks (same pipe-list format).
blocks = unit_row.Blocks(:).';
n_blocks = numel(blocks);

if ~isfield(unit_row, 'Perturbation') || isempty(unit_row.Perturbation)
    unit_row.Perturbation = zeros(1, n_blocks);
    return;
end

perturbations = unit_row.Perturbation(:).';
if isscalar(perturbations)
    perturbations = repmat(perturbations, 1, n_blocks);
elseif numel(perturbations) ~= n_blocks
    error('epp_sorting_table_to_units:PerturbationBlockMismatch', ...
        'Unit %s: %d Perturbation value(s) for %d Blocks.', ...
        unit_row.Neuron_ID, numel(perturbations), n_blocks);
end

unit_row.Blocks = blocks;
unit_row.Perturbation = perturbations;
end

function blocks = parse_blocks_cell(value)
if isnumeric(value) || islogical(value)
    blocks = value;
    return;
end
if ~(ischar(value) || isstring(value))
    blocks = value;
    return;
end

s = strtrim(char(value));
if isempty(s)
    blocks = [];
    return;
end

tokens = regexp(s, '[|\s;]+', 'split');
tokens = tokens(~cellfun(@isempty, tokens));
nums = str2double(tokens);
if any(isnan(nums))
    blocks = value;
else
    blocks = nums(:).';
end
end

function [column_type, converted_col] = infer_and_convert_column_type(col_cells)
n_rows = numel(col_cells);
logical_col = false(n_rows, 1);
numeric_col = nan(n_rows, 1);
logical_ok = true;
numeric_ok = true;

for i = 1:numel(col_cells)
    value = col_cells{i};
    if islogical(value) && isscalar(value)
        logical_col(i) = value;
    elseif isnumeric(value) && isscalar(value) && any(value == [0, 1])
        logical_col(i) = logical(value);
    elseif ischar(value) 
        s_log = lower(strtrim(char(value)));
        if strcmp(s_log, 'true') || strcmp(s_log, '1')
            logical_col(i) = true;
        elseif strcmp(s_log, 'false') || strcmp(s_log, '0')
            logical_col(i) = false;
        else
            logical_ok = false;
        end
    else
        logical_ok = false;
    end

    if isnumeric(value) && isscalar(value)
        numeric_col(i) = double(value);
    elseif islogical(value) && isscalar(value)
        numeric_col(i) = double(value);
    elseif ischar(value) || isstring(value)
        s_num = strtrim(char(value));
        if isempty(s_num)
            numeric_col(i) = NaN;
        else
            num = str2double(s_num);
            if isnan(num) && ~strcmpi(s_num, 'nan')
                numeric_ok = false;
            else
                numeric_col(i) = num;
            end
        end
    else
        numeric_ok = false;
    end
end

if logical_ok
    column_type = 'logical';
    converted_col = logical_col;
elseif numeric_ok
    column_type = 'numeric';
    converted_col = numeric_col;
else
    column_type = 'original';
    converted_col = [];
end
end

