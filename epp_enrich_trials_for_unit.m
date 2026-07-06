function trials = epp_enrich_trials_for_unit(trials, unit_row)
% Flip lateral trial variables by recorded hemisphere (epp_build_rasters).
%
% Reads unit_row.Hemisphere from final_sorting. Perturbation and block are
% set earlier during block alignment in load_unit_combined.

if isempty(trials)
    return;
end

hemisphere_value = get_unit_table_field(unit_row, 'Hemisphere', []);
rec_hemi_sign = parse_recorded_hemisphere(hemisphere_value);

if isnan(rec_hemi_sign)
    error('epp_enrich_trials_for_unit:InvalidHemisphere', ...
        'Unit %s has missing or invalid hemisphere. Expected left (L/l) or right (R/r).', ...
        unit_row.Neuron_ID);
end

target_pos_fields = {'tar_pos', 'nct_pos', 'fix_pos', 'reach_tar_pos_closest', 'saccade_tar_pos_closest', 'hemifield'};
hand_fields = {'demanded_hand', 'used_hand', 'reach_hand'};

for f = target_pos_fields
    if ~isfield(trials, f{1}), continue; end
    IN = [trials.(f{:})];
    if isreal(IN)
        OUT = num2cell(real(IN) * rec_hemi_sign);
    else
        OUT = num2cell(real(IN) * rec_hemi_sign + 1i * imag(IN));
    end
    [trials.(f{:})] = deal(OUT{:});
end

for f = hand_fields
    if ~isfield(trials, f{1}), continue; end
    IN = [trials.(f{:})];
    OUT = IN;
    OUT(IN == 1) = -1; % left hand
    OUT(IN == 2) = 1;  % right hand
    OUT = num2cell(OUT * rec_hemi_sign);
    [trials.(f{:})] = deal(OUT{:});
end
end

function value = get_unit_table_field(unit_row, field_name, default_value)
if isfield(unit_row, field_name) && ~isempty(unit_row.(field_name))
    value = unit_row.(field_name);
else
    value = default_value;
end
end

function rec_hemi_sign = parse_recorded_hemisphere(hemisphere_value)
% L/l -> +1, R/r -> -1 (multiplier for lateral flip).
rec_hemi_sign = NaN;
if isempty(hemisphere_value)
    return;
end
hemisphere_text = lower(hemisphere_value);
last_char = hemisphere_text(end);
if last_char == 'l'
    rec_hemi_sign = 1;
elseif last_char == 'r'
    rec_hemi_sign = -1;
end
end
