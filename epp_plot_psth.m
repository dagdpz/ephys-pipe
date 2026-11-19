function epp_plot_psth(cfg, UI)
%EPP_PLOT_PSTH  Stage 4: condition-split PSTH and binned raster for one unit.
%
%   epp_plot_psth(cfg, unit_row)
%
%   Inputs:
%     cfg - struct from epp_load_cfg.
%     unit_row - one element of unit_info. Loads windows + behavior from
%           cfg.roots.raster/<Neuron_ID>_raster.mat (epochs skipped).
%           Uses windows(wi).cond_mask for grouping. Returns without writing
%           if that file is missing.
%           Uses cfg.CONDITIONS for colors and names, and cfg.psth.condition_subplots
%           for axes. condition_subplots empty: all conditions on one axes.
%           A cell of index vectors: one vertical subplot per cell.
%           Condition means are computed here from the aligned raster.
%           A trial that matches several conditions is drawn in each of them.
%
%   Outputs:
%     none. For each unit that has a raster file, writes
%       cfg.roots.psth/<Neuron_ID>_psth.png
%       cfg.roots.psth/<Neuron_ID>_psth.mat  (variable psth_data)
%     psth_data.Neuron_ID, psth_data.windows(wi).name, .align_state,
%     .conditions (name, n_trials, bin_centers_s, mean_rate_hz, sem_rate_hz).
%     The raster, trials, and behavior traces are not written again.
%     Windows are concatenated on one axes (0.1 s gap). X tick labels are
%     time within each window. The y-axis is PSTH (top half), rasters (next
%     quarter), then eye and hand traces (bottom quarter). Y ticks are rate
%     and stay in the PSTH half. The number at the trace height is the
%     position span from each trace center to the edge of its band, snapped
%     to 1, 2, or 5 times a power of ten. Eye and hand traces of one condition
%     are a single line in that condition's color. Raster squares are one
%     per occupied bin. A trial is drawn in every condition it matches.
%
% Stage 4: load aligned rasters, group by condition, plot PSTH + rasters.

Neuron_ID = UI.Neuron_ID;
raster_file = fullfile(cfg.roots.raster, sprintf('%s_raster.mat', Neuron_ID));
if ~isfile(raster_file)
    return
end
% Plot needs traces; stats can load epochs without windows/behavior.
loaded = load(raster_file, 'windows', 'behavior');
windows = loaded.windows;
behavior = loaded.behavior;

n_windows = numel(windows);
n_cond = numel(cfg.CONDITIONS);
conv_kernel = psth_conv_kernel(cfg);
epp_ensure_paths();

fig = figure('units','normalized','outerposition',[0 0 1 1],'color','w');
cond_groups = {1:n_cond};
if ~isempty(cfg.psth.condition_subplots)
    cond_groups = cfg.psth.condition_subplots;
end
n_groups = numel(cond_groups);

raster_trial_height = 0.1;
window_gap_s = 0.1;
label_name = {windows.name};

for gi = 1:n_groups
    group = cond_groups{gi}(:).';

    ax = subplot(n_groups, 1, gi);
    hold(ax, 'on');
    state_separator = 0;
    y_max = 0;
    y_min = 0;
    align_x = zeros(n_windows, 1);
    label_x = zeros(n_windows, 1);
    legend_h = gobjects(numel(cfg.CONDITIONS), 1);
    raster_h = gobjects(0);
    traces_h = gobjects(0);

    for wi = 1:n_windows
        w = windows(wi);
        n_events = size(w.raster, 1);

        core_cols = (w.psth_pad_bins + 1):(size(w.raster, 2) - w.psth_pad_bins);

        % Same seconds-per-pixel in every window. 0.1 s gap is the separator.
        state_shift = state_separator - w.t_start_s;
        state_separator = state_shift + w.t_end_s + window_gap_s;
        align_x(wi) = state_shift;
        label_x(wi) = state_shift + (w.t_start_s + w.t_end_s)/2;
        x_win = state_shift + [w.t_start_s w.t_end_s];

        if n_events == 0 || ~isfield(w, 'cond_mask') || isempty(w.cond_mask)
            continue;
        end

        raster_counts = double(w.raster);

        raster_y = 0;

        for ci = group
            if ci > size(w.cond_mask, 2)
                continue
            end
            keep_mask = w.cond_mask(:, ci);
            if ~any(keep_mask)
                continue
            end
            cond_def = cfg.CONDITIONS(ci);
            color = cond_def.color;
            if max(color) > 1, color = color / 255; end

            % rasters
            rows = find(keep_mask);
            n_rows = numel(rows);

            raster_h(end+1, 1) = line(x_win, [raster_y raster_y], 'Color', color, 'LineWidth', 0.5); %#ok<AGROW>
            xs = [];
            ys = [];
            for r = 1:n_rows
                spike_cols = find(raster_counts(rows(r), core_cols) > 0);
                if isempty(spike_cols)
                    continue
                end
                centers = w.bin_centers(core_cols);
                xs = [xs, centers(spike_cols(:)).' + state_shift]; %#ok<AGROW>
                ys = [ys, (raster_y - r*raster_trial_height) * ones(1, numel(spike_cols))]; %#ok<AGROW>
            end
            if ~isempty(xs)
                raster_h(end+1, 1) = scatter(ax, xs, ys, 10, color, 's', 'filled'); %#ok<AGROW>
            end
            raster_y = raster_y - n_rows*raster_trial_height;
            raster_h(end+1, 1) = line(x_win, [raster_y raster_y], 'Color', color, 'LineWidth', 0.5); %#ok<AGROW>
            if wi == 1
                raster_h(end+1, 1) = text(x_win(1), raster_y + n_rows*raster_trial_height, sprintf('N=%d', n_rows), 'Color', color, 'FontSize', 5, 'VerticalAlignment', 'top'); %#ok<AGROW>
            end
            y_min = min(y_min, raster_y);

            % Smoothed PSTH. Pad columns feed conv, then are dropped.
            psth_edges = w.t_start_s:cfg.psth.bin_size_s:w.t_end_s;
            if numel(psth_edges) < 2 || psth_edges(end) < w.t_end_s
                psth_edges(end+1) = w.t_end_s; %#ok<AGROW>
            end
            psth_centers = psth_edges(1:end-1).' + diff(psth_edges(1:2))/2;
            n_psth_bins = numel(psth_centers);
            raster_to_psth_bin = discretize(w.bin_centers(core_cols), psth_edges);
            cond_rate = raster_counts(keep_mask, :) ./ w.bin_size_s;
            cond_rate = conv2(cond_rate, conv_kernel, 'same');
            cond_rate = cond_rate(:, core_cols);
            cond_rate_psth = NaN(size(cond_rate, 1), n_psth_bins);

            for bi = 1:n_psth_bins
                src_cols = raster_to_psth_bin == bi;
                if any(src_cols)
                    cond_rate_psth(:, bi) = mean(cond_rate(:, src_cols), 2, 'omitnan');
                end
            end
            mean_rate = mean(cond_rate_psth, 1, 'omitnan');
            sem_rate = std(cond_rate_psth, 0, 1, 'omitnan') ./ sqrt(size(cond_rate_psth, 1));

            peak = max(mean_rate + sem_rate, [], 'omitnan');
            if isfinite(peak)
                y_max = max(y_max, peak);
            end

            H = shadedErrorBar(psth_centers(:).' + state_shift, mean_rate(:).', sem_rate(:).', {'color', color, 'linewidth', 1.5}, 1);
            if ~isgraphics(legend_h(ci))
                if isstruct(H) && isfield(H, 'mainLine')
                    legend_h(ci) = H.mainLine;
                elseif isgraphics(H)
                    legend_h(ci) = H;
                end
            end

            % Top-level behavior(wi), parallel to windows — optional if absent.
            if wi > numel(behavior) || isempty(behavior(wi).data)
                continue
            end
            beh = behavior(wi).data;
            ch = behavior(wi).channels;
            i_eye = [find(strcmp(ch, 'x_eye'), 1) find(strcmp(ch, 'y_eye'), 1)];
            i_hnd = [find(strcmp(ch, 'x_hnd'), 1) find(strcmp(ch, 'y_hnd'), 1)];
            xx = w.bin_edges(:);
            xs = [];
            ys = [];
            is_hand = false(0, 1);
            for k = find(keep_mask).'
                if k > size(beh, 1)
                    continue
                end
                segs = {};
                hand_flag = false(0, 1);
                for ii = i_eye
                    segs{end+1} = double(squeeze(beh(k, :, ii))).'; %#ok<AGROW>
                    hand_flag(end+1) = false; %#ok<AGROW>
                end
                for ii = i_hnd
                    segs{end+1} = double(squeeze(beh(k, :, ii))).'; %#ok<AGROW>
                    hand_flag(end+1) = true; %#ok<AGROW>
                end
                for s = 1:numel(segs)
                    v = segs{s};
                    if numel(v) ~= numel(xx)
                        continue
                    end
                    xs = [xs; xx + align_x(wi); NaN]; %#ok<AGROW>
                    ys = [ys; v; NaN]; %#ok<AGROW>
                    is_hand = [is_hand; repmat(hand_flag(s), numel(xx), 1); false]; %#ok<AGROW>
                end
            end
            if isempty(xs)
                continue
            end
            h = line(ax, xs, ys, 'Color', color, 'LineWidth', 0.5, 'HandleVisibility', 'off');
            h.UserData = struct('is_hand', is_hand);
            traces_h(end+1, 1) = h; %#ok<AGROW>
        end
    end

    % PSTH uses [0, y_max] (top half). Scale rasters into the next quarter,
    % then eye/hand into the bottom quarter.
    if ~(y_max > 0)
        y_max = 1;
    end
    raster_span = y_max / 2;
    raster_depth = -y_min;
    if raster_depth > 0
        scale_y_objects(raster_h, raster_span / raster_depth, 0, []);
    end

    trace_top = -raster_span;
    trace_bot = -y_max;
    band_half = (trace_top - trace_bot) / 4;
    eye_center = trace_top - band_half;
    hand_center = trace_bot + band_half;
    pos_peak = 0;
    for hi = 1:numel(traces_h)
        if ~isgraphics(traces_h(hi))
            continue
        end
        y = traces_h(hi).YData;
        y = y(isfinite(y));
        if ~isempty(y)
            pos_peak = max(pos_peak, max(abs(y)));
        end
    end
    trace_span = round_trace_span(pos_peak);
    scale_y_objects(traces_h, band_half / trace_span, eye_center, hand_center);

    if pos_peak > 0
        text(ax, 0, (trace_top + trace_bot)/2, sprintf('%g', trace_span), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
            'Clipping', 'off', 'FontSize', 8, 'Color', 'k');
    end


    % rest of formatting

    %ticks
    xlim(ax, [0, state_separator - window_gap_s]);
    ylim(ax, [-y_max, y_max]);
    yt = get(ax, 'YTick');
    set(ax, 'YTick', yt(yt >= 0));
    xtick = [];
    xticklabel = {};
    tick_step = 0.2;
    for wi = 1:n_windows
        w = windows(wi);
        rel = ceil(w.t_start_s/tick_step)*tick_step : tick_step : floor(w.t_end_s/tick_step)*tick_step;
        rel = round(rel/tick_step)*tick_step;
        if w.t_start_s <= 0 && w.t_end_s >= 0
            rel = unique([rel, 0]);
        end
        xtick = [xtick, align_x(wi) + rel]; %#ok<AGROW>
        xticklabel = [xticklabel, arrayfun(@(t) sprintf('%g', t), rel, 'UniformOutput', false)]; %#ok<AGROW>
    end
    set(ax, 'XTick', xtick, 'XTickLabel', xticklabel);


    % event lines, labels, titles, legend
    for wi = 1:n_windows
        line(ax, [align_x(wi) align_x(wi)], [-y_max y_max], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
        text(ax, label_x(wi), y_max - 0.05*y_max, label_name{wi}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Interpreter', 'none');
    end
    if gi == n_groups
        xlabel(ax, 'Time (s)');
    end
    ylabel(ax, 'Rate (Hz)');
    if gi == 1
        blocks_txt = strtrim(sprintf('%g ', UI.Blocks));
        title(ax, sprintf('%s  %s  ch %g  blocks %s  z %g', Neuron_ID, UI.Target, UI.Channel, blocks_txt, UI.z), 'Interpreter', 'none');
    end
    if any(isgraphics(legend_h))
        legend(ax, legend_h(isgraphics(legend_h)), {cfg.CONDITIONS(isgraphics(legend_h)).name}, 'Interpreter', 'none', 'Location', 'northeast');
    end
    hold(ax, 'off');
end

% saveas re-renders every spike marker. The figure is already drawn.
drawnow;
frame = getframe(fig);
imwrite(frame.cdata, fullfile(cfg.roots.psth, sprintf('%s_psth.png', Neuron_ID)));
close(fig);
end

function conv_kernel = psth_conv_kernel(cfg)
% Same half-width as psth_smooth_pad_bins in epp_build_rasters.
bin_size_s = cfg.raster_bin_size_s;
switch lower(char(cfg.psth.smoothing_kernel))
    case 'gaussian'
        sigma_bins = max(cfg.psth.smoothing_width_s / bin_size_s, eps);
        radius_bins = max(1, ceil(4 * sigma_bins));
        x = -radius_bins:radius_bins;
        conv_kernel = exp(-0.5 * (x / sigma_bins).^2);
        conv_kernel = conv_kernel ./ sum(conv_kernel);
    case 'box'
        box_bins = max(1, round(cfg.psth.smoothing_width_s / bin_size_s));
        conv_kernel = ones(1, box_bins) / box_bins;
    otherwise
        conv_kernel = 1;
end
end

function scale_y_objects(handles, gain, eye_center, hand_center)
% Multiply YData by gain. Traces with UserData.is_hand also get a band offset.
% Text uses Position(2). hand_center empty: pure gain (rasters).
for hi = 1:numel(handles)
    h = handles(hi);
    if ~isgraphics(h)
        continue
    end
    if isprop(h, 'YData')
        y = h.YData(:);
        if isempty(hand_center)
            y = y * gain;
        else
            hand = false(size(y));
            if isstruct(h.UserData) && isfield(h.UserData, 'is_hand') ...
                    && numel(h.UserData.is_hand) == numel(y)
                hand = logical(h.UserData.is_hand(:));
            end
            y(~hand) = y(~hand) * gain + eye_center;
            y(hand) = y(hand) * gain + hand_center;
        end
        h.YData = y.';
    elseif isprop(h, 'Position')
        pos = h.Position;
        pos(2) = pos(2) * gain;
        if ~isempty(hand_center)
            pos(2) = pos(2) + eye_center;
        end
        h.Position = pos;
    end
end
end

function span = round_trace_span(pos_peak)
% Smallest 1, 2, or 5 times a power of ten that covers pos_peak.
if ~(pos_peak > 0) || ~isfinite(pos_peak)
    span = 1;
    return
end
e = floor(log10(pos_peak));
steps = [1 2 5 10] * 10^e;
fit = steps(steps >= pos_peak);
if isempty(fit)
    span = 10^(e + 1);
else
    span = min(fit);
end
end
