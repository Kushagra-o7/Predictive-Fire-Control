function figs = plot_scenario_results(res, m, P, OPT)
%PLOT_SCENARIO_RESULTS  Standard figure set for one engagement.
%
%   Figure 1 - Engagement: geometry with laser pulses, elevation angles,
%              beam-to-target pointing error, tracking errors.
%   Figure 2 - Estimation and drive: estimation errors with 2-sigma bounds,
%              prediction error, motor voltage/current, turret rate,
%              target range/speed/acceleration vs envelope limits.

if nargin < 4, OPT = struct(); end
vis = 'on'; if isfield(OPT, 'fig_visible'), vis = OPT.fig_visible; end
t = res.t; tr = res.truth; sc = res.sc;
d2r = 180/pi;
cTrue = [0 0 0]; cMeas = [0.6 0.6 0.6]; cEst = [0 0.45 0.74];
cPred = [0.85 0.33 0.10]; cTur = [0.47 0.67 0.19]; cHit = [0 0.6 0]; cMiss = [0.85 0 0];
pl = res.pulses;

%% ---------------- Figure 1: engagement ---------------------------------
f1 = figure('Name', [sc.name ' - engagement'], 'Color', 'w', 'Visible', vis, ...
            'Position', [60 60 1250 820]);
sgtitle(sprintf('%s  |  %s', sc.name, summary_line(m)), 'FontWeight', 'bold');

% (a) geometry
subplot(2,2,1); hold on; grid on; box on;
th = linspace(P.theta_min, P.theta_max, 200);
plot(P.r_min*cos(th), P.r_min*sin(th), ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
plot(P.r_max*cos(th), P.r_max*sin(th), ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
v = res.valid == 1;
plot(res.r_meas(v).*cos(res.theta_meas(v)), res.r_meas(v).*sin(res.theta_meas(v)), '.', ...
     'Color', cMeas, 'MarkerSize', 4, 'DisplayName', 'radar meas.');
plot(tr.x, tr.y, '-', 'Color', cTrue, 'LineWidth', 1.2, 'DisplayName', 'true');
plot(res.x_cur(1,:), res.x_cur(4,:), '--', 'Color', cEst, 'LineWidth', 1.0, 'DisplayName', 'estimate');
plot(res.x_pred(1,:), res.x_pred(4,:), '-.', 'Color', cPred, 'DisplayName', 'predicted (t+T_p)');
for k = 1:numel(pl)
    L = 1.05 * hypot(pl(k).target(1), pl(k).target(2));
    c = cMiss; if pl(k).hit, c = cHit; end
    plot([0 L*cos(pl(k).theta)], [0 L*sin(pl(k).theta)], '-', 'Color', [c 0.6], 'HandleVisibility', 'off');
    plot(pl(k).target(1), pl(k).target(2), 'o', 'Color', c, 'MarkerFaceColor', c, 'HandleVisibility', 'off');
end
plot(0, 0, 'k^', 'MarkerFaceColor', 'k', 'DisplayName', 'turret');
axis equal; xlabel('x (ground range) [m]'); ylabel('y (altitude) [m]');
title('Elevation-plane geometry (green = hit, red = miss)');
legend('Location', 'best', 'FontSize', 7);

% (b) elevation angles
subplot(2,2,2); hold on; grid on; box on;
shade_blackout(sc, [-20 95]);
plot(t(v), d2r*res.theta_meas(v), '.', 'Color', cMeas, 'MarkerSize', 4, 'DisplayName', '\theta_{meas}');
plot(t, d2r*tr.theta, '-', 'Color', cTrue, 'LineWidth', 1.2, 'DisplayName', '\theta_{target}');
plot(t, d2r*res.theta_ref, '-', 'Color', cPred, 'DisplayName', '\theta_{ref} (predicted)');
plot(t, d2r*res.theta, '-', 'Color', cTur, 'LineWidth', 1.4, 'DisplayName', '\theta_{turret}');
mark_pulses(pl, d2r*[pl.theta], cHit, cMiss);
yl = d2r*[min([tr.theta res.theta]) max([tr.theta res.theta_ref])];
ylim(yl + [-3 3]); xlim([0 t(end)]);
xlabel('t [s]'); ylabel('elevation [deg]'); title('Elevation angles');
legend('Location', 'best', 'FontSize', 7);

% (c) pointing error
subplot(2,2,3); hold on; grid on; box on;
yl = [-1 1] * max(6, min(30, 1.2*m.max_point_m));
shade_blackout(sc, yl);
patch([0 t(end) t(end) 0], P.R_eng*[-1 -1 1 1], [0 0.6 0], 'FaceAlpha', 0.12, ...
      'EdgeColor', 'none', 'DisplayName', '\pm R_{eng}');
plot(t, res.point_err_m, '-', 'Color', cTur, 'LineWidth', 1.1, 'DisplayName', 'true beam-target distance');
plot(t, res.d_hat, ':', 'Color', cEst, 'DisplayName', 'estimated d_{hat}');
for k = 1:numel(pl)
    c = cMiss; if pl(k).hit, c = cHit; end
    plot(pl(k).t, pl(k).miss_signed, 'o', 'Color', c, 'MarkerFaceColor', c, 'HandleVisibility', 'off');
end
ylim(yl); xlim([0 t(end)]);
xlabel('t [s]'); ylabel('[m]'); title('Pointing error (signed cross-beam distance); markers = pulses');
legend('Location', 'best', 'FontSize', 7);

% (d) angle errors
subplot(2,2,4); hold on; grid on; box on;
e_track = d2r*(res.theta_ref - res.theta);
e_point = d2r*(tr.theta - res.theta);
ss = t >= 2.5;
yl = [-1 1] * max(0.5, 1.2*max(abs([e_track(ss) e_point(ss)])));
shade_blackout(sc, yl);
plot(t, e_track, '-', 'Color', cPred, 'DisplayName', '\theta_{ref} - \theta_{turret} (tracking)');
plot(t, e_point, '-', 'Color', cTur, 'LineWidth', 1.1, 'DisplayName', '\theta_{target} - \theta_{turret} (pointing)');
ylim(yl); xlim([0 t(end)]);
xlabel('t [s]'); ylabel('[deg]'); title('Tracking vs pointing error');
legend('Location', 'best', 'FontSize', 7);

%% ---------------- Figure 2: estimation & drive -------------------------
f2 = figure('Name', [sc.name ' - estimation & drive'], 'Color', 'w', 'Visible', vis, ...
            'Position', [80 40 1250 860]);
sgtitle(sprintf('%s  -  estimation, prediction and drive', sc.name), 'FontWeight', 'bold');
ok = ~isnan(res.x_cur(1,:));

subplot(3,2,1); hold on; grid on; box on;
ex = res.x_cur(1,:) - tr.x;  ey = res.x_cur(4,:) - tr.y;
yl = [-1 1] * max(2, min(40, 1.3*max(abs([ex(t>0.5) ey(t>0.5)]))));
shade_blackout(sc, yl);
plot(t, ex, 'Color', cEst, 'DisplayName', 'x err');
plot(t, ey, 'Color', cPred, 'DisplayName', 'y err');
plot(t,  2*res.pos_std(1,:), ':', 'Color', cEst, 'HandleVisibility', 'off');
plot(t, -2*res.pos_std(1,:), ':', 'Color', cEst, 'HandleVisibility', 'off');
plot(t,  2*res.pos_std(2,:), ':', 'Color', cPred, 'HandleVisibility', 'off');
plot(t, -2*res.pos_std(2,:), ':', 'Color', cPred, 'HandleVisibility', 'off');
ylim(yl); xlim([0 t(end)]); ylabel('[m]'); title('Position estimation error (dotted: \pm2\sigma)');
legend('Location', 'best', 'FontSize', 7);

subplot(3,2,2); hold on; grid on; box on;
evx = res.x_cur(2,:) - tr.vx;  evy = res.x_cur(5,:) - tr.vy;
yl = [-1 1] * max(2, min(60, 1.3*max(abs([evx(t>1) evy(t>1)]))));
shade_blackout(sc, yl);
plot(t, evx, 'Color', cEst, 'DisplayName', 'v_x err');
plot(t, evy, 'Color', cPred, 'DisplayName', 'v_y err');
plot(t,  2*res.vel_std(1,:), ':', 'Color', cEst, 'HandleVisibility', 'off');
plot(t, -2*res.vel_std(1,:), ':', 'Color', cEst, 'HandleVisibility', 'off');
ylim(yl); xlim([0 t(end)]); ylabel('[m/s]'); title('Velocity estimation error');
legend('Location', 'best', 'FontSize', 7);

subplot(3,2,3); hold on; grid on; box on;
tp_ok = ok & t + P.Tp <= t(end);
trp = target_truth(t(tp_ok) + P.Tp, sc);
epx = nan(size(t)); epy = epx;
epx(tp_ok) = res.x_pred(1,tp_ok) - trp.x;  epy(tp_ok) = res.x_pred(4,tp_ok) - trp.y;
yl = [-1 1] * max(2, min(60, 1.3*max(abs([epx(t>1) epy(t>1)]))));
shade_blackout(sc, yl);
plot(t, epx, 'Color', cEst, 'DisplayName', 'x_{pred}(t) - x(t+T_p)');
plot(t, epy, 'Color', cPred, 'DisplayName', 'y_{pred}(t) - y(t+T_p)');
ylim(yl); xlim([0 t(end)]); ylabel('[m]'); title(sprintf('Prediction error, T_p = %.2f s', P.Tp));
legend('Location', 'best', 'FontSize', 7);

subplot(3,2,4); hold on; grid on; box on;
yyaxis left;  plot(t, res.V, 'LineWidth', 1.0); ylabel('V_m [V]');
yyaxis right; plot(t, res.i, 'LineWidth', 1.0); ylabel('i [A]');
xlim([0 t(end)]); title('Motor voltage and current');

subplot(3,2,5); hold on; grid on; box on;
plot(t, d2r*res.theta_dot, 'Color', cTur, 'LineWidth', 1.1, 'DisplayName', 'turret');
dth = gradient(tr.theta, t);
plot(t, d2r*dth, 'k--', 'DisplayName', 'target LOS rate');
xlim([0 t(end)]); xlabel('t [s]'); ylabel('[deg/s]'); title('Angular rates');
legend('Location', 'best', 'FontSize', 7);

subplot(3,2,6); hold on; grid on; box on;
yyaxis left;
plot(t, tr.speed, '-', 'DisplayName', '|v|'); plot(t, 10*tr.accel, '--', 'DisplayName', '10|a|');
yline(P.v_max, ':', 'v_{max}', 'HandleVisibility', 'off'); yline(10*P.a_max, ':', '10 a_{max}', 'HandleVisibility', 'off');
ylabel('[m/s], [0.1 m/s^2]');
yyaxis right;
plot(t, tr.r, '-', 'DisplayName', 'r'); ylabel('r [m]'); ylim([0 2200]);
xlim([0 t(end)]); xlabel('t [s]'); title('Target kinematics vs envelope');
legend('Location', 'best', 'FontSize', 7);

figs = [f1 f2];
if isfield(OPT, 'save_figures') && OPT.save_figures
    save_fig(f1, fullfile(OPT.fig_dir, sprintf('case_%s_engagement', sc.id)));
    save_fig(f2, fullfile(OPT.fig_dir, sprintf('case_%s_estimation', sc.id)));
end
end

%% ------------------------------------------------------------------------
function s = summary_line(m)
if m.success
    s = sprintf('INTERCEPT at %.2f s, %d/%d hits, RMS pointing %.2f m', ...
        m.t_intercept, m.n_hits, m.n_pulses, m.rms_point_m);
else
    s = sprintf('NO INTERCEPT, %d pulses, RMS pointing %.2f m', m.n_pulses, m.rms_point_m);
end
end

function shade_blackout(sc, yl)
if isempty(sc.blackout), return; end
t0 = sc.blackout(1); t1 = t0 + sc.blackout(2);
patch([t0 t1 t1 t0], [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], 'FaceAlpha', 0.15, ...
      'EdgeColor', 'none', 'DisplayName', 'radar loss');
end

function mark_pulses(pl, y, cHit, cMiss)
for k = 1:numel(pl)
    c = cMiss; if pl(k).hit, c = cHit; end
    plot(pl(k).t, y(k), 'p', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 9, 'HandleVisibility', 'off');
end
end
