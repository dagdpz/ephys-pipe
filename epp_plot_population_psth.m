function epp_plot_population_psth(cfg)
% Stage 5: aggregate per-unit PSTHs into population mean +/- SEM plots.

unit_psths = load_unit_psth_data(cfg.roots.psth);
if isempty(unit_psths)
    warning('epp_plot_population_psth:NoData', 'No unit PSTH files found in %s.', cfg.roots.psth);
    return;
end

n_windows = numel(unit_psths(1).windows);
if n_windows == 0
    return;
end

n_conditions = numel(cfg.CONDITIONS);
pop_data = struct();
pop_data.n_units = numel(unit_psths);
pop_data.Neuron_ID = {unit_psths.Neuron_ID};
pop_data.windows = repmat(struct( ...
    'name', '', 'align_state', [], 'conditions', []), n_windows, 1);

fig = figure('Color', 'w');
for wi = 1:n_windows
    w_ref = unit_psths(1).windows(wi);
    pop_data.windows(wi).name = w_ref.name;
    pop_data.windows(wi).align_state = w_ref.align_state;
    pop_data.windows(wi).conditions = repmat(struct( ...
        'name', '', 'n_units', 0, 'n_trials_total', 0, ...
        'bin_centers_s', [], 'mean_rate_hz', [], 'sem_rate_hz', []), n_conditions, 1);

    ax = subplot(n_windows, 1, wi);
    hold(ax, 'on');

    for ci = 1:n_conditions
        [unit_means, bin_centers, n_trials_vec] = collect_condition_psths(unit_psths, wi, ci);
        if isempty(unit_means)
            continue;
        end

        pop_mean = mean(unit_means, 1, 'omitnan');
        n_units_ci = size(unit_means, 1);
        pop_sem = std(unit_means, 0, 1, 'omitnan') ./ sqrt(max(n_units_ci, 1));

        pop_data.windows(wi).conditions(ci).name = cfg.CONDITIONS(ci).name;
        pop_data.windows(wi).conditions(ci).n_units = n_units_ci;
        pop_data.windows(wi).conditions(ci).n_trials_total = sum(n_trials_vec);
        pop_data.windows(wi).conditions(ci).bin_centers_s = bin_centers;
        pop_data.windows(wi).conditions(ci).mean_rate_hz = pop_mean;
        pop_data.windows(wi).conditions(ci).sem_rate_hz = pop_sem;

        color = normalize_plot_color(cfg.CONDITIONS(ci).color);
        x = bin_centers(:).';
        y = pop_mean(:).';
        sem = pop_sem(:).';
        
        lineProps={'LineWidth', 1.5, 'Color', color};
        shadedErrorBar(x,y,sem,lineProps,0);
    end

    hold(ax, 'off');
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Rate (Hz)');
    title(ax, sprintf('%s (%d)', w_ref.name, w_ref.align_state));
    legend(ax, {cfg.CONDITIONS.name}, 'Interpreter', 'none', 'Location', 'best');
end

saveas(fig, fullfile(cfg.roots.population_psth, 'population_psth.png'));
close(fig);
save(fullfile(cfg.roots.population_psth, 'population_psth.mat'), 'pop_data');
end

function unit_psths = load_unit_psth_data(psth_root)
unit_psths = struct([]);
psth_files = dir(fullfile(psth_root, '*_psth.mat'));
for i = 1:numel(psth_files)
    if strcmp(psth_files(i).name, 'population_psth.mat')
        continue;
    end
    loaded = load(fullfile(psth_root, psth_files(i).name), 'psth_data');
    unit_psths = [unit_psths; loaded.psth_data]; %#ok<AGROW>
end
end

function [unit_means, bin_centers, n_trials_vec] = collect_condition_psths(unit_psths, wi, ci)
unit_means = [];
bin_centers = [];
n_trials_vec = [];

for u = 1:numel(unit_psths)
    if numel(unit_psths(u).windows) < wi
        continue;
    end
    conds = unit_psths(u).windows(wi).conditions;
    if numel(conds) < ci || isempty(conds(ci).mean_rate_hz)
        continue;
    end
    if isempty(bin_centers)
        bin_centers = conds(ci).bin_centers_s(:).';
    end
    unit_means = [unit_means; conds(ci).mean_rate_hz(:).']; %#ok<AGROW>
    n_trials_vec = [n_trials_vec; conds(ci).n_trials]; %#ok<AGROW>
end
end

function color = normalize_plot_color(color_in)
color = color_in(:).';
if numel(color) ~= 3
    color = [0 0 0];
end
if max(color) > 1
    color = color / 255;
end
end
