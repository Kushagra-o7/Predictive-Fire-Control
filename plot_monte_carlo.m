function fig = plot_monte_carlo(mc, P, OPT)
%PLOT_MONTE_CARLO  Generalization study over random admissible engagements.
%   mc: struct array with fields sc, chk, m

vis = 'on'; if isfield(OPT, 'fig_visible'), vis = OPT.fig_visible; end
M = [mc.m]; C = [mc.chk];
succ = [M.success];

fig = figure('Name', 'Monte-Carlo generalization', 'Color', 'w', 'Visible', vis, ...
             'Position', [120 90 1250 700]);
sgtitle(sprintf(['Monte-Carlo: %d random admissible cases - success %.1f %%, ' ...
    'median intercept %.2f s'], numel(mc), 100*mean(succ), median([M.t_intercept], 'omitnan')), ...
    'FontWeight', 'bold');

subplot(2,3,1);
histogram([M.t_intercept], 0:0.5:P.t_end); grid on;
xlabel('first interception time [s]'); ylabel('cases'); title('Interception time');

subplot(2,3,2);
allmiss = [];
for k = 1:numel(mc), allmiss = [allmiss [mc(k).pulses.miss]]; end %#ok<AGROW>
histogram(allmiss(allmiss < 10), 0:0.25:10); grid on; hold on;
xline(P.R_eng, 'r--', 'R_{eng}');
xlabel('miss distance [m]'); ylabel('pulses');
title(sprintf('All pulses: %.1f %% within R_{eng}', 100*mean(allmiss <= P.R_eng)));

subplot(2,3,3);
scatter([C.a_max], [M.rms_point_m], 25, succ, 'filled'); grid on;
colormap(gca, [0.85 0 0; 0 0.6 0]); clim([0 1]);
xlabel('max target accel [m/s^2]'); ylabel('RMS pointing [m]');
title('Pointing vs maneuver intensity (green = success)');

subplot(2,3,4);
scatter([C.r_min], [M.t_intercept], 25, 'filled'); grid on;
xlabel('min range [m]'); ylabel('interception time [s]'); title('Interception time vs range');

subplot(2,3,5);
bar([mean([M.n_hits]) mean([M.n_pulses])]); grid on;
set(gca, 'XTickLabel', {'hits / case', 'pulses / case'}); title('Mean pulses per case');

subplot(2,3,6);
bo = arrayfun(@(s) ~isempty(s.sc.blackout), mc);
rec = [M.bo_recovery];
histogram(rec(bo), 0:0.1:2); grid on;
xlabel('pointing recovery after radar loss [s]'); ylabel('cases');
title(sprintf('Radar-loss cases: %d', sum(bo)));

if isfield(OPT, 'save_figures') && OPT.save_figures
    save_fig(fig, fullfile(OPT.fig_dir, 'monte_carlo'));
end
end
