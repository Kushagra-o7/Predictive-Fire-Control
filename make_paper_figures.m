%MAKE_PAPER_FIGURES  Compact IEEE-column figures and LaTeX table rows for the paper.
%
%   Run main_guardian first (with run_sensitivity = true) so that the files
%   results/results_cases.mat, results_mc.mat, results_comparison.mat and
%   results_sensitivity.mat exist. Then run:
%       >> make_paper_figures
%   Output: paper/figures/*.pdf (vector) and paper/tables/*_rows.tex
%   NOTE: regenerating overwrites the *_rows.tex files - edit captions and
%   table headers in paper/main.tex, not in the generated row files.

clearvars; close all;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here);
rd   = fullfile(here, 'results');
fd   = fullfile(here, 'paper', 'figures');
td   = fullfile(here, 'paper', 'tables');
if ~exist(fd, 'dir'), mkdir(fd); end
if ~exist(td, 'dir'), mkdir(td); end

P = guardian_params();
S = load(fullfile(rd, 'results_cases.mat'));      results = S.results;
M = load(fullfile(rd, 'results_mc.mat'));         mc = M.mc;
C = load(fullfile(rd, 'results_comparison.mat')); cmp = C.comparison; cids = C.case_ids;
E = load(fullfile(rd, 'results_sensitivity.mat')); sens = E.sens;

W  = 3.5;                                   % IEEE column width [in]
FS = 8;                                     % font size
set(groot, 'defaultAxesFontSize', FS, 'defaultTextFontSize', FS, ...
    'defaultAxesFontName', 'Times New Roman', 'defaultTextFontName', 'Times New Roman', ...
    'defaultLineLineWidth', 1.0, 'defaultAxesTickLabelInterpreter', 'tex', ...
    'defaultLegendFontSize', FS-1);
cF = [0 0.447 0.741]; cT = [0.85 0.325 0.098]; cK = [0 0 0]; cG = [0.466 0.674 0.188];
cHit = [0 0.6 0]; cMiss = [0.85 0 0]; grey = [0.55 0.55 0.55];

D = analyze_controller_design(P, struct('fig_visible', 'off', 'save_figures', false));
close(D.fig);
w = D.w; dF = D.design(1); dT = D.design(2);

%% ---------------- Fig: Bode of L(s) -------------------------------------
f = newfig(W, 2.55);
Gp = polyval(P.Gp.num, 1i*w) ./ polyval(P.Gp.den, 1i*w);
subplot(2,1,1); hold on; grid on; box on;
semilogx(w, 20*log10(abs(Gp)), ':', 'Color', grey);
semilogx(w, 20*log10(abs(dF.L)), 'Color', cF);
semilogx(w, 20*log10(abs(dT.L)), '--', 'Color', cT);
set(gca, 'XScale', 'log'); xlim([1e-2 1e3]); ylim([-140 80]); yline(0, 'k-', 'HandleVisibility', 'off');
plot(dF.w_gc, 0, 'o', 'Color', cF, 'MarkerFaceColor', cF, 'HandleVisibility', 'off');
ylabel('Magnitude [dB]');
legend({'G_p(s)', 'L(s): delay-matched (final)', 'L(s): tight PI (comparison)'}, 'Location', 'southwest');
subplot(2,1,2); hold on; grid on; box on;
phF = rad2deg(unwrap(angle(dF.L))); phT = rad2deg(unwrap(angle(dT.L)));
semilogx(w, phF, 'Color', cF); semilogx(w, phT, '--', 'Color', cT);
set(gca, 'XScale', 'log'); xlim([1e-2 1e3]); ylim([-360 -60]); yline(-180, 'k-');
plot(dF.w_gc, -180 + dF.PM_deg, 'o', 'Color', cF, 'MarkerFaceColor', cF);
plot(dF.w_pc, -180, 's', 'Color', cF, 'MarkerFaceColor', cF);
text(0.012, -300, sprintf('final: PM = %.1f^\\circ @ %.2f rad/s, GM = %.1f dB @ %.1f rad/s', ...
    dF.PM_deg, dF.w_gc, dF.GM_dB, dF.w_pc), 'Color', cF);
text(0.012, -335, sprintf('tight: PM = %.1f^\\circ @ %.2f rad/s, GM = %.1f dB @ %.1f rad/s', ...
    dT.PM_deg, dT.w_gc, dT.GM_dB, dT.w_pc), 'Color', cT);
xlabel('\omega [rad/s]'); ylabel('Phase [deg]');
export(f, fd, 'fig_bode');

%% ---------------- Fig: root locus ---------------------------------------
f = newfig(W, 1.7);
kk = [logspace(-3, log10(300), 1500)];
nr = numel(roots(dF.denL));
R = zeros(nr, numel(kk));
for j = 1:numel(kk), R(:, j) = roots(polyadd(dF.denL, kk(j)*dF.numL)); end
ol = roots(dF.denL); cl = dF.cl_poles;
kcrit = 10^(dF.GM_dB/20);
rc = roots(polyadd(dF.denL, kcrit*dF.numL));
for pnl = 1:2
    subplot(1,2,pnl); hold on; grid on; box on;
    plot(real(R).', imag(R).', '.', 'Color', grey, 'MarkerSize', 2);
    plot(real(ol), imag(ol), 'kx', 'MarkerSize', 7, 'LineWidth', 1.2);
    plot(real(cl), imag(cl), 'p', 'Color', cF, 'MarkerFaceColor', cF, 'MarkerSize', 8);
    plot(real(rc(abs(real(rc)) < 1e-2)), imag(rc(abs(real(rc)) < 1e-2)), 'o', 'Color', cT, 'MarkerSize', 6, 'LineWidth', 1.2);
    xline(0, 'k-');
    xlabel('Re(s)'); ylabel('Im(s)');
    if pnl == 1, xlim([-12 4]); ylim([-70 70]); title('(a) full');
    else,        xlim([-6 1]);  ylim([-4 4]);  title('(b) dominant poles');
    end
end
export(f, fd, 'fig_rlocus');

%% ---------------- Fig: delay matching -----------------------------------
f = newfig(W, 2.8);
subplot(2,1,1); hold on; grid on; box on;
semilogx(w, 20*log10(abs(1 - dF.T .* exp(1i*w*P.Tp))), 'Color', cF);
semilogx(w, 20*log10(abs(1 - dF.T)), ':', 'Color', cF);
semilogx(w, 20*log10(abs(1 - dT.T .* exp(1i*w*P.Tp))), '--', 'Color', cT);
semilogx(w, 20*log10(abs(1 - dT.T)), '-.', 'Color', cT);
set(gca, 'XScale', 'log'); xlim([1e-2 1e2]); ylim([-120 20]);
xlabel('\omega [rad/s]'); ylabel('|E/\theta_t| [dB]');
legend({'final, T_p=0.5 s', 'final, no prediction', 'tight, T_p=0.5 s', 'tight, no prediction'}, ...
    'Location', 'southeast');
title('(a) pointing-error transfer 1 - T(s)e^{sT_p}');
subplot(2,1,2); hold on; grid on; box on;
dt = 1e-3; t = 0:dt:6;
th = @(tt) deg2rad(5)*tt + 0.5*deg2rad(1)*tt.^2;
clF = polyadd(dF.denL, dF.numL); clT = polyadd(dT.denL, dT.numL);
plot(t, rad2deg(th(t) - lsim_local(dF.numL, clF, th(t + P.Tp), dt)), 'Color', cF);
plot(t, rad2deg(th(t) - lsim_local(dF.numL, clF, th(t), dt)), ':', 'Color', cF);
plot(t, rad2deg(th(t) - lsim_local(dT.numL, clT, th(t + P.Tp), dt)), '--', 'Color', cT);
plot(t, rad2deg(th(t) - lsim_local(dT.numL, clT, th(t), dt)), '-.', 'Color', cT);
ylim([-4 4]); xlabel('t [s]'); ylabel('e_{point} [deg]');
title('(b) target angle \theta_t = 5t + 0.5t^2 deg, ideal prediction');
export(f, fd, 'fig_delay');

%% ---------------- Fig: step response & voltage --------------------------
f = newfig(W, 2.25);
subplot(2,1,1); hold on; grid on; box on;
yF = lsim_local(dF.numL, clF, ones(size(t)), dt);
yT = lsim_local(dT.numL, clT, ones(size(t)), dt);
plot(t, yF, 'Color', cF); plot(t, yT, '--', 'Color', cT);
xlim([0 4]); ylabel('\theta/\theta_{ref}');
legend({sprintf('final: t_r=%.2f s, t_s=%.2f s, M_p=%.1f%%', dF.step), ...
        sprintf('tight: t_r=%.2f s, t_s=%.2f s, M_p=%.1f%%', dT.step)}, 'Location', 'southeast');
title('(a) closed-loop unit-step response');
subplot(2,1,2); hold on; grid on; box on;
th0 = deg2rad(35);
VF = lsim_local(P.Gc.num, P.Gc.den, th0*(1 - yF), dt);
Pt = P; Pt.ctrl_mode = 'tight'; Pt = design_compensator(Pt);
VT = lsim_local(Pt.Gc.num, Pt.Gc.den, th0*(1 - yT), dt);
plot(t, VF, 'Color', cF); plot(t, VT, '--', 'Color', cT);
xlim([0 4]); xlabel('t [s]'); ylabel('V_m [V]');
title(sprintf('(b) motor voltage, 35^\\circ slew (peak %.0f V vs %.0f V)', max(abs(VF)), max(abs(VT))));
export(f, fd, 'fig_step');

%% ---------------- Fig: scenario B time histories ------------------------
rB = get_case(results, 'B');
f = newfig(W, 2.7);
subplot(2,1,1); hold on; grid on; box on;
plot(rB.t, rad2deg(rB.truth.theta), 'Color', cK);
plot(rB.t, rad2deg(rB.theta_ref), 'Color', cT);
plot(rB.t, rad2deg(rB.theta), 'Color', cG, 'LineWidth', 1.3);
pl = rB.pulses;
plot([pl.t], rad2deg([pl.theta]), 'p', 'Color', cHit, 'MarkerFaceColor', cHit, 'MarkerSize', 6);
ylim([0 60]); xlim([0 10]); ylabel('elevation [deg]');
legend({'\theta_t (true)', '\theta_{ref} (predicted)', '\theta (turret)', 'laser pulse'}, 'Location', 'southeast');
title('(a) scenario B: elevation angles');
subplot(2,1,2); hold on; grid on; box on;
patch([0 10 10 0], P.R_eng*[-1 -1 1 1], cHit, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
plot(rB.t, rB.point_err_m, 'Color', cG);
plot(rB.t, rB.d_hat, ':', 'Color', cF);
for k = 1:numel(pl)
    c = cMiss; if pl(k).hit, c = cHit; end
    plot(pl(k).t, pl(k).miss_signed, 'o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4);
end
ylim([-6 10]); xlim([0 10]); xlabel('t [s]'); ylabel('[m]');
legend({'\pm R_{eng}', 'true distance', 'estimated |d|'}, 'Location', 'north', 'Orientation', 'horizontal');
title('(b) scenario B: pointing error and pulses');
export(f, fd, 'fig_caseB');

%% ---------------- Fig: scenario C (measurement loss) --------------------
rC = get_case(results, 'C');
f = newfig(W, 3.05);
bo = rC.sc.blackout; tb = [bo(1) bo(1)+bo(2)];
subplot(3,1,1); hold on; grid on; box on;
shade(tb, [-30 30]);
patch([0 10 10 0], P.R_eng*[-1 -1 1 1], cHit, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
plot(rC.t, rC.point_err_m, 'Color', cG);
pl = rC.pulses;
for k = 1:numel(pl)
    c = cMiss; if pl(k).hit, c = cHit; end
    plot(pl(k).t, pl(k).miss_signed, 'o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4);
end
ylim([-30 30]); xlim([0 10]); ylabel('[m]'); title('(a) scenario C: pointing error (shaded: radar loss)');
subplot(3,1,2); hold on; grid on; box on;
ex = rC.x_cur(1,:) - rC.truth.x; ey = rC.x_cur(4,:) - rC.truth.y;
shade(tb, [-4 4]);
plot(rC.t, ex, 'Color', cF); plot(rC.t, ey, 'Color', cT);
plot(rC.t, 2*rC.pos_std(1,:), ':', 'Color', cF); plot(rC.t, -2*rC.pos_std(1,:), ':', 'Color', cF);
plot(rC.t, 2*rC.pos_std(2,:), ':', 'Color', cT); plot(rC.t, -2*rC.pos_std(2,:), ':', 'Color', cT);
ylim([-4 4]); xlim([0 10]); ylabel('[m]'); title('(b) position estimation error, \pm2\sigma dotted');
legend({'x', 'y'}, 'Location', 'northeast', 'Orientation', 'horizontal');
subplot(3,1,3); hold on; grid on; box on;
sp = nan(size(rC.t)); ok = ~isnan(rC.sigma_perp); sp(ok) = rC.sigma_perp(ok);
shade(tb, [0 2]);
plot(rC.t, sp, 'Color', cF); yline(P.fire.sigma_max, 'r--');
ylim([0 2]); xlim([0 10]); xlabel('t [s]'); ylabel('\sigma_\perp [m]');
title('(c) cross-beam std used by the fire logic (dashed: gate)');
export(f, fd, 'fig_caseC');

%% ---------------- Fig: estimation & prediction (scenario B) -------------
f = newfig(W, 2.25);
subplot(2,1,1); hold on; grid on; box on;
evx = rB.x_cur(2,:) - rB.truth.vx; evy = rB.x_cur(5,:) - rB.truth.vy;
plot(rB.t, evx, 'Color', cF); plot(rB.t, evy, 'Color', cT);
plot(rB.t, 2*rB.vel_std(1,:), ':', 'Color', cF); plot(rB.t, -2*rB.vel_std(1,:), ':', 'Color', cF);
ylim([-8 8]); xlim([0 10]); ylabel('[m/s]');
legend({'v_x error', 'v_y error', '\pm2\sigma'}, 'Location', 'northeast', 'Orientation', 'horizontal');
title('(a) scenario B: velocity estimation error');
subplot(2,1,2); hold on; grid on; box on;
okp = ~isnan(rB.x_pred(1,:)) & rB.t + P.Tp <= rB.t(end);
trp = target_truth(rB.t(okp) + P.Tp, rB.sc);
epx = nan(size(rB.t)); epy = epx;
epx(okp) = rB.x_pred(1,okp) - trp.x; epy(okp) = rB.x_pred(4,okp) - trp.y;
plot(rB.t, epx, 'Color', cF); plot(rB.t, epy, 'Color', cT);
ylim([-8 8]); xlim([0 10]); xlabel('t [s]'); ylabel('[m]');
legend({'x', 'y'}, 'Location', 'northeast', 'Orientation', 'horizontal');
title('(b) prediction error p_{pred}(t) - p(t+T_p)');
export(f, fd, 'fig_estB');

%% ---------------- Fig: Monte Carlo ---------------------------------------
Mm = [mc.m];
f = newfig(W, 1.6);
subplot(1,2,1); grid on; box on; hold on;
histogram([Mm.t_intercept], 0.75:0.25:3.0, 'FaceColor', cF);
xlabel('first interception [s]'); ylabel('cases');
title(sprintf('(a) %d/%d cases intercepted', sum([Mm.success]), numel(Mm)));
subplot(1,2,2); grid on; box on; hold on;
allmiss = [];
for k = 1:numel(mc), allmiss = [allmiss [mc(k).pulses.miss]]; end %#ok<AGROW>
histogram(allmiss, 0:0.2:3, 'FaceColor', cG);
xline(P.R_eng, 'r--', 'LineWidth', 1.2);
xlabel('pulse miss distance [m]'); ylabel('pulses');
title(sprintf('(b) %d pulses, %.1f%% hits', numel(allmiss), 100*mean(allmiss <= P.R_eng)));
export(f, fd, 'fig_mc');

%% ---------------- Table rows ---------------------------------------------
fid = fopen(fullfile(td, 'cases_rows.tex'), 'w');
for k = 1:numel(results)
    m = results(k).m;
    bo = '--'; if ~isnan(m.bo_recovery), bo = sprintf('%.2f', m.bo_recovery); end
    fprintf(fid, '%s & %.2f & %d/%d & %.2f & %.2f & %.2f & %.2f & %.1f & %s & %.0f \\\\\n', ...
        results(k).id, m.t_intercept, m.n_hits, m.n_pulses, m.miss_mean, m.rms_point_m, ...
        m.rms_pos_err, m.rms_pred_err, 100*m.frac_in_R, bo, m.V_peak);
end
fclose(fid);

fid = fopen(fullfile(td, 'comparison_rows.tex'), 'w');
labels = {'(1) Delay-matched $G_c$, $T_p=0.5$ s (final)', '(2) Tight PI, $T_p=0.5$ s', ...
          '(3) Tight PI, no prediction', '(4) Delay-matched $G_c$, no prediction'};
for v = 1:numel(cmp)
    fprintf(fid, '%s & %d/%d & %.2f & %.2f & %d/%d & %.0f \\\\\n', labels{v}, ...
        sum(cmp(v).success), numel(cids), mean(cmp(v).t_int, 'omitnan'), mean(cmp(v).rms_point), ...
        sum(cmp(v).hits), sum(cmp(v).pulses), max(cmp(v).V_peak));
end
fclose(fid);

fid = fopen(fullfile(td, 'sensitivity_rows.tex'), 'w');
for a = 1:numel(sens.sigma_m)
    fprintf(fid, '%g', sens.sigma_m(a));
    for b = 1:numel(sens.tau_m)
        cellv = sprintf('%.1f / %.2f', 100*sens.inR(a,b), sens.miss(a,b));
        if sens.sigma_m(a) == P.est.sigma_m && sens.tau_m(b) == P.est.tau_m
            cellv = ['\textbf{' cellv '}'];
        end
        fprintf(fid, ' & %s', cellv);
    end
    fprintf(fid, ' \\\\\n');
end
fclose(fid);

% Monte-Carlo and design numbers quoted in the text (for checking)
fprintf('\nMC: success %d/%d, median t_int %.2f, max %.2f, pulses %d, hit %.1f%%, mean RMS point %.2f m\n', ...
    sum([Mm.success]), numel(Mm), median([Mm.t_intercept], 'omitnan'), max([Mm.t_intercept]), ...
    numel(allmiss), 100*mean(allmiss <= P.R_eng), mean([Mm.rms_point_m]));
nbo = sum(arrayfun(@(s) ~isempty(s.sc.blackout), mc));
fprintf('MC: %d of %d cases include radar loss; mean pulse miss %.2f m\n', nbo, numel(mc), mean(allmiss));
fprintf('Figures written to %s\nTables written to %s\n', fd, td);

%% ========================================================================
function f = newfig(w, h)
f = figure('Units', 'inches', 'Position', [1 1 w h], 'Color', 'w', 'Visible', 'off');
end

function export(f, fd, name)
exportgraphics(f, fullfile(fd, [name '.pdf']), 'ContentType', 'vector');
exportgraphics(f, fullfile(fd, [name '.png']), 'Resolution', 300);
close(f);
end

function r = get_case(results, id)
r = results(strcmp({results.id}, id)).res;
r.m = results(strcmp({results.id}, id)).m;
end

function shade(tb, yl)
patch([tb(1) tb(2) tb(2) tb(1)], [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], ...
      'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function c = polyadd(a, b)
n = max(numel(a), numel(b));
c = [zeros(1, n - numel(a)) a] + [zeros(1, n - numel(b)) b];
end

function y = lsim_local(num, den, u, dt)
num = num / den(1); den = den / den(1);
n = numel(den) - 1;
num = [zeros(1, n + 1 - numel(num)) num];
Dd = num(1); r = num(2:end) - Dd*den(2:end);
A = [zeros(n-1,1) eye(n-1); -fliplr(den(2:end))];
B = [zeros(n-1,1); 1];  Cc = fliplr(r);
Ex = expm([A B; zeros(1, n+1)] * dt);
Ad = Ex(1:n,1:n); Bd = Ex(1:n,n+1);
x = zeros(n,1); y = zeros(size(u));
for k = 1:numel(u)
    y(k) = Cc*x + Dd*u(k);
    x = Ad*x + Bd*u(k);
end
end
