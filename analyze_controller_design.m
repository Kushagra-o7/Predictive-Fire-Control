function D = analyze_controller_design(P, OPT)
%ANALYZE_CONTROLLER_DESIGN  Analytical (frequency/time-domain) study of Gc(s).
%
%   D = analyze_controller_design(P, OPT) prints and plots:
%     - plant poles, compensator, open-loop L(s) = Gc Gp
%     - gain/phase margins, crossover, closed-loop poles, damping, bandwidth
%     - error constants and the delay-matching property T(s) exp(sTp) ~ 1
%     - Bode plot, root locus (vs. gain), step response (acquisition), ramp
%       and constant-angular-acceleration pointing error, noise transmission
%   Compares the final 'delay_matched' design with the 'tight' design.
%   Toolbox-free (works without the Control System Toolbox).

if nargin < 2, OPT = struct(); end
vis = 'on'; if isfield(OPT, 'fig_visible'), vis = OPT.fig_visible; end

Pd = P;  Pd.ctrl_mode = 'delay_matched'; Pd = design_compensator(Pd);
Pt = P;  Pt.ctrl_mode = 'tight';         Pt = design_compensator(Pt);
designs = {Pd, Pt};  names = {'delay-matched (final)', 'tight (comparison)'};

fprintf('\n================ CONTROLLER DESIGN ANALYSIS ================\n');
fprintf('Plant  Gp(s) = %g / ( s (%g s^2 + %g s + %g) )\n', P.K_N, P.d2, P.d1, P.d0);
pp = roots(P.Gp.den);
wn_res = abs(pp(abs(pp) > 1)); zeta_res = -real(pp(abs(pp) > 1)) ./ wn_res;
fprintf('  plant poles: %s\n', mat2str(pp.', 4));
fprintf('  electromechanical resonance: wn = %.2f rad/s, zeta = %.3f\n', wn_res(1), zeta_res(1));
fprintf('  low-frequency model: Gp ~ (K_N/d0)/s = %.4g/s  (rad/s per V)\n', P.K_N/P.d0);

w = logspace(-2, 3, 4000);
D = struct();
D.w = w;
for k = 1:2
    Pk = designs{k};
    numL = conv(Pk.Gc.num, P.Gp.num);
    denL = conv(Pk.Gc.den, P.Gp.den);
    L = polyval(numL, 1i*w) ./ polyval(denL, 1i*w);
    T = L ./ (1 + L);  S = 1 ./ (1 + L);
    [GM, PM, wgc, wpc] = margins_local(w, L);
    clpoles = roots(polyadd(denL, numL));
    bw = w(find(abs(T) < 1/sqrt(2), 1));
    % error constants
    nInt = sum(abs(roots(denL)) < 1e-9);
    Kv = Inf; if nInt == 1, Kv = lim_sL(numL, denL); end
    % pointing-error sensitivity with prediction: 1 - T exp(jwTp)
    Ep = abs(1 - T .* exp(1i*w*P.Tp));
    E0 = abs(1 - T);
    d = struct('name', names{k}, 'Gc', Pk.Gc, 'numL', numL, 'denL', denL, 'GM_dB', 20*log10(GM), ...
        'PM_deg', PM, 'w_gc', wgc, 'w_pc', wpc, 'cl_poles', clpoles, 'bw', bw, 'Kv', Kv, ...
        'L', L, 'T', T, 'S', S, 'Ep', Ep, 'E0', E0);
    D.design(k) = d;

    fprintf('\n--- %s:  %s\n', names{k}, Pk.Gc.description);
    fprintf('  gain margin  %.1f dB at %.1f rad/s,  phase margin %.1f deg at %.2f rad/s\n', ...
        d.GM_dB, wpc, PM, wgc);
    fprintf('  closed-loop bandwidth (-3 dB) %.2f rad/s,  Kv = %.3f 1/s (ramp lag 1/Kv = %.3f s)\n', ...
        bw, Kv, 1/Kv);
    for p = clpoles.'
        if imag(p) >= 0
            fprintf('    CL pole %8.3f %+8.3fj   wn = %6.2f  zeta = %.3f\n', real(p), imag(p), ...
                abs(p), -real(p)/abs(p));
        end
    end
    fprintf('  stable: %d\n', all(real(clpoles) < 0));
end

% Delay-matching coefficients (series of 1/T in s) for the final design
Kc = Pd.Gc.Kc; tc = Pd.Gc.tau_c;
c1 = P.d0/(P.K_N*Kc); c2 = (P.d1 + tc*P.d0)/(P.K_N*Kc); c3 = (P.d2 + tc*P.d1)/(P.K_N*Kc);
fprintf('\nDelay matching (final design): 1/T(s) = 1 + c1 s + c2 s^2 + c3 s^3 + ...\n');
fprintf('  c1 = %.4f  (target Tp     = %.4f)\n', c1, P.Tp);
fprintf('  c2 = %.4f  (target Tp^2/2 = %.4f)\n', c2, P.Tp^2/2);
fprintf('  c3 = %.5f  (Tp^3/6 = %.5f) -> residual pointing error ~ %.4f * d^3(theta)/dt^3\n', ...
    c3, P.Tp^3/6, P.Tp^3/6 - c3);
D.c = [c1 c2 c3];

%% ---------------- Figures ------------------------------------------------
fig = figure('Name', 'Controller design analysis', 'Color', 'w', 'Visible', vis, ...
             'Position', [40 40 1300 860]);
sgtitle('Turret elevation loop - analytical design', 'FontWeight', 'bold');
cols = {[0 0.45 0.74], [0.85 0.33 0.10]};

% Bode magnitude / phase of L
subplot(3,3,1); hold on; grid on; box on;
Gp = polyval(P.Gp.num, 1i*w) ./ polyval(P.Gp.den, 1i*w);
semilogx(w, 20*log10(abs(Gp)), 'k:', 'DisplayName', 'G_p');
for k = 1:2, semilogx(w, 20*log10(abs(D.design(k).L)), 'Color', cols{k}, 'DisplayName', ['L: ' names{k}]); end
set(gca, 'XScale', 'log'); yline(0, 'k-', 'HandleVisibility', 'off'); ylim([-150 80]);
ylabel('|L| [dB]'); title('Open-loop magnitude'); legend('FontSize', 7, 'Location', 'southwest');
subplot(3,3,4); hold on; grid on; box on;
for k = 1:2, semilogx(w, rad2deg(unwrap(angle(D.design(k).L))), 'Color', cols{k}); end
set(gca, 'XScale', 'log'); yline(-180, 'k-'); ylim([-360 -60]);
xlabel('\omega [rad/s]'); ylabel('\angle L [deg]');
title(sprintf('PM = %.0f / %.0f deg, GM = %.0f / %.0f dB', D.design(1).PM_deg, ...
    D.design(2).PM_deg, D.design(1).GM_dB, D.design(2).GM_dB));

% Root locus of the final design vs gain multiplier
subplot(3,3,7); hold on; grid on; box on;
kk = logspace(-2, log10(30), 400);
R = zeros(numel(roots(D.design(1).denL)), numel(kk));
for j = 1:numel(kk), R(:,j) = roots(polyadd(D.design(1).denL, kk(j)*D.design(1).numL)); end
plot(real(R), imag(R), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 3);
p0 = D.design(1).cl_poles;
plot(real(p0), imag(p0), 'p', 'Color', cols{1}, 'MarkerFaceColor', cols{1}, 'MarkerSize', 10);
plot(real(roots(D.design(1).denL)), imag(roots(D.design(1).denL)), 'kx', 'MarkerSize', 9);
xlim([-40 5]); ylim([-70 70]); xlabel('Re'); ylabel('Im');
title('Root locus vs gain multiplier (stars = design poles)');

% Closed-loop |T| and sensitivity
subplot(3,3,2); hold on; grid on; box on;
for k = 1:2
    semilogx(w, 20*log10(abs(D.design(k).T)), 'Color', cols{k});
    semilogx(w, 20*log10(abs(D.design(k).S)), '--', 'Color', cols{k});
end
set(gca, 'XScale', 'log'); ylim([-80 10]); xlabel('\omega [rad/s]'); ylabel('[dB]');
title('|T| (solid, noise transmission) and |S| (dashed)');

% Pointing-error transfer with prediction
subplot(3,3,3); hold on; grid on; box on;
for k = 1:2
    semilogx(w, 20*log10(D.design(k).Ep), 'Color', cols{k}, 'DisplayName', ['|1-Te^{sT_p}| ' names{k}]);
end
semilogx(w, 20*log10(D.design(1).E0), ':', 'Color', cols{1}, 'DisplayName', '|1-T| final, no prediction');
set(gca, 'XScale', 'log'); xlim([1e-2 1e2]); ylim([-120 20]);
xlabel('\omega [rad/s]'); ylabel('[dB]'); legend('FontSize', 6, 'Location', 'southeast');
title('Pointing-error transfer \theta_{target} \rightarrow e_{point}');

% Time responses
dt = 1e-3; t = 0:dt:6;
subplot(3,3,5); hold on; grid on; box on;
for k = 1:2
    y = lsim_local(D.design(k).numL, polyadd(D.design(k).denL, D.design(k).numL), ones(size(t)), dt);
    plot(t, y, 'Color', cols{k});
    si = stepinfo_local(t, y);
    fprintf('Step (%s): rise %.2f s, settling(2%%) %.2f s, overshoot %.1f %%\n', names{k}, si(1), si(2), si(3));
    D.design(k).step = si;
end
xlabel('t [s]'); ylabel('\theta / \theta_{ref}'); title('Closed-loop step (acquisition slew)');

% Pointing error for ramp and constant-acceleration target angles,
% reference = exact prediction theta(t + Tp)
subplot(3,3,6); hold on; grid on; box on;
w0 = deg2rad(5);  a0 = deg2rad(1);                      % 5 deg/s, 1 deg/s^2
th_true = @(tt) w0*tt + 0.5*a0*tt.^2;
for k = 1:2
    y = lsim_local(D.design(k).numL, polyadd(D.design(k).denL, D.design(k).numL), th_true(t + P.Tp), dt);
    plot(t, rad2deg(th_true(t) - y), 'Color', cols{k});
end
ylim([-1.5 1.5]); xlabel('t [s]'); ylabel('[deg]');
title('Pointing error, \theta = 5t + 0.5t^2 deg (ideal prediction)');

% Voltage demand for the acquisition slew of the final design
subplot(3,3,8); hold on; grid on; box on;
th0 = deg2rad(35);
y = lsim_local(D.design(1).numL, polyadd(D.design(1).denL, D.design(1).numL), th0*ones(size(t)), dt);
V = lsim_local(Pd.Gc.num, Pd.Gc.den, th0 - y, dt);
plot(t, V, 'Color', cols{1}); xlabel('t [s]'); ylabel('V_m [V]');
title('Voltage for a 35 deg acquisition slew (final)');

% Effect of prediction horizon mismatch
subplot(3,3,9); hold on; grid on; box on;
Tps = [0 0.25 0.5 0.75];
for j = 1:numel(Tps)
    semilogx(w, 20*log10(abs(1 - D.design(1).T .* exp(1i*w*Tps(j)))), 'DisplayName', sprintf('T_p = %.2f', Tps(j)));
end
set(gca, 'XScale', 'log'); xlim([1e-2 1e2]); ylim([-100 20]); xlabel('\omega [rad/s]');
ylabel('[dB]'); legend('FontSize', 7, 'Location', 'southeast');
title('Final loop: pointing error vs horizon used');

D.fig = fig;
if isfield(OPT, 'save_figures') && OPT.save_figures
    save_fig(fig, fullfile(OPT.fig_dir, 'controller_design'));
end
fprintf('=============================================================\n');
end

%% ========================================================================
function c = polyadd(a, b)
n = max(numel(a), numel(b));
c = [zeros(1, n - numel(a)) a] + [zeros(1, n - numel(b)) b];
end

function Kv = lim_sL(num, den)
% lim s->0 of s*num/den for a den with exactly one root at s = 0
den_red = den(1:end-1);              % divide by s (last coefficient is 0)
Kv = polyval(num, 0) / polyval(den_red, 0);
end

function [GM, PM, wgc, wpc] = margins_local(w, L)
mag = abs(L); ph = rad2deg(unwrap(angle(L)));
k = find(mag(1:end-1) >= 1 & mag(2:end) < 1, 1);
if isempty(k), wgc = NaN; PM = Inf;
else
    f = (0 - log(mag(k))) / (log(mag(k+1)) - log(mag(k)));
    wgc = exp(log(w(k)) + f*(log(w(k+1)) - log(w(k))));
    PM = 180 + interp1(w(k:k+1), ph(k:k+1), wgc);
end
k = find(ph(1:end-1) > -180 & ph(2:end) <= -180, 1);
if isempty(k), wpc = NaN; GM = Inf;
else
    wpc = interp1(ph(k:k+1), w(k:k+1), -180);
    GM = 1 / interp1(w(k:k+1), mag(k:k+1), wpc);
end
end

function y = lsim_local(num, den, u, dt)
% Simulate a proper SISO transfer function (ZOH, exact discretization)
num = num / den(1); den = den / den(1);
n = numel(den) - 1;
num = [zeros(1, n + 1 - numel(num)) num];
Dd = num(1); r = num(2:end) - Dd*den(2:end);
A = [zeros(n-1,1) eye(n-1); -fliplr(den(2:end))];
B = [zeros(n-1,1); 1];  C = fliplr(r);
E = expm([A B; zeros(1, n+1)] * dt);
Ad = E(1:n,1:n); Bd = E(1:n,n+1);
x = zeros(n,1); y = zeros(size(u));
for k = 1:numel(u)
    y(k) = C*x + Dd*u(k);
    x = Ad*x + Bd*u(k);
end
end

function si = stepinfo_local(t, y)
yf = y(end);
tr10 = t(find(y >= 0.1*yf, 1)); tr90 = t(find(y >= 0.9*yf, 1));
outside = find(abs(y - yf) > 0.02*abs(yf), 1, 'last');
ts = t(min(outside + 1, numel(t)));
si = [tr90 - tr10, ts, 100*max(0, (max(y) - yf)/yf)];
end
