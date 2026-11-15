function fig = plot_summary(results, P, OPT, fig_name)
%PLOT_SUMMARY  Cross-scenario comparison of the finalized design.
%   results: struct array with fields id, name, m (metrics from compute_metrics)

if nargin < 4, fig_name = 'summary_cases'; end
vis = 'on'; if isfield(OPT, 'fig_visible'), vis = OPT.fig_visible; end
ids = {results.id};
M = [results.m];
x = 1:numel(results);

fig = figure('Name', 'Scenario summary', 'Color', 'w', 'Visible', vis, 'Position', [100 80 1250 700]);
sgtitle('Finalized design - performance across all deterministic cases', 'FontWeight', 'bold');

subplot(2,3,1);
ti = [M.t_intercept]; ti_plot = ti; ti_plot(isnan(ti)) = P.t_end;
b = bar(x, ti_plot); b.FaceColor = 'flat';
for k = x, if isnan(ti(k)), b.CData(k,:) = [0.85 0 0]; else, b.CData(k,:) = [0 0.45 0.74]; end, end
set(gca, 'XTick', x, 'XTickLabel', ids); grid on;
ylabel('[s]'); title('First interception time (red = none)'); ylim([0 P.t_end]);

subplot(2,3,2);
bar(x, [[M.n_hits]; [M.n_pulses] - [M.n_hits]]', 'stacked');
set(gca, 'XTick', x, 'XTickLabel', ids); grid on;
legend({'hits', 'misses'}, 'Location', 'best'); title('Laser pulses');

subplot(2,3,3);
bar(x, [[M.rms_point_m]; [M.max_point_m]]');
yline(P.R_eng, 'r--', 'R_{eng}');
set(gca, 'XTick', x, 'XTickLabel', ids); grid on;
legend({'RMS', 'max'}, 'Location', 'best'); ylabel('[m]');
title('Pointing error (t \geq 3 s)');

subplot(2,3,4);
bar(x, [[M.miss_min]; [M.miss_mean]]');
yline(P.R_eng, 'r--', 'R_{eng}');
set(gca, 'XTick', x, 'XTickLabel', ids); grid on;
legend({'min', 'mean'}, 'Location', 'best'); ylabel('[m]'); title('Pulse miss distance');

subplot(2,3,5);
bar(x, [[M.rms_pos_err]; [M.rms_pred_err]]');
set(gca, 'XTick', x, 'XTickLabel', ids); grid on;
legend({'estimate (t)', 'prediction (t+T_p)'}, 'Location', 'best'); ylabel('[m]');
title('RMS position error');

subplot(2,3,6);
bar(x, 100*[M.frac_in_R]);
set(gca, 'XTick', x, 'XTickLabel', ids); grid on; ylim([0 100]);
ylabel('[%]'); title('Time with beam inside R_{eng} (t \geq 3 s)');

if isfield(OPT, 'save_figures') && OPT.save_figures
    save_fig(fig, fullfile(OPT.fig_dir, fig_name));
end
end
