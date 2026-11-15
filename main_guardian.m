%MAIN_GUARDIAN  Master script - Guardian point-defense, ME319 Assignment 2.
%
%   This file owns the simulation clock. For every test case it
%     1. steps time at OPT.dt (1 kHz) from 0 to t_end,
%     2. calls target_motion_generator  -> radar measurement (20 Hz, noisy, loss)
%     3. calls turret_controller        -> motor voltage + fire command
%     4. models the laser (0->1 edge, 1 s recharge, ray hit test vs TRUE target)
%     5. integrates turret_plant with RK4 (ZOH voltage) + mechanical stops,
%   then computes metrics, plots and summary tables.
%
%   Sections (enable/disable in USER OPTIONS):
%     A. Analytical controller design (Bode, margins, root locus, ...)
%     B. Envelope verification of every test case
%     C. Standardized scenarios A, B, C
%     D. Extended stress cases D - I (different maneuvers / radar loss)
%     E. Monte-Carlo generalization over random admissible engagements
%     F. Architecture comparison (final design vs alternatives)
%     G. Save results (MAT + CSV + figures)
%
%   Run:  >> main_guardian
%   Batch / quick runs: define a struct OPT_OVERRIDE in the workspace before
%   running, e.g.  OPT_OVERRIDE.run_monte_carlo = false; main_guardian

clearvars -except OPT_OVERRIDE; close all; clc;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(here);
clear turret_controller target_motion_generator turret_plant   % fresh persistent state

%% ============================ USER OPTIONS ================================
OPT.dt                  = 1e-3;     % simulation / controller step         [s]
OPT.t_end               = 10;       % engagement duration                  [s]
OPT.theta0              = 0;        % initial turret elevation (stowed)    [rad]
OPT.theta_dot0          = 0;        % initial turret rate                  [rad/s]
OPT.i0                  = 0;        % initial motor current                [A]

OPT.run_design_analysis = true;     % Section A
OPT.run_standard        = true;     % Section C  (A, B, C)
OPT.run_extended        = true;     % Section D  (D ... I)
OPT.run_monte_carlo     = true;     % Section E
OPT.n_monte_carlo       = 50;       %   number of random cases
OPT.run_comparison      = true;     % Section F
OPT.run_sensitivity     = false;    % Section H (slow: estimator tuning study)
OPT.sens_sigma_m        = [4 8 12 16];  % Singer RMS accel candidates   [m/s^2]
OPT.sens_tau_m          = [1 2 4];      % Singer correlation times      [s]
OPT.n_sens_cases        = 20;       %   tuning cases (seeds 1001.., disjoint from MC)
OPT.P_mod               = [];       % optional @(P) ... to modify parameters
OPT.plot_each_case      = true;     % per-case figure sets
OPT.save_figures        = true;
OPT.fig_visible         = 'on';     % 'off' for batch runs
OPT.results_dir         = fullfile(here, 'results');
OPT.fig_dir             = fullfile(OPT.results_dir, 'figures');

if exist('OPT_OVERRIDE', 'var') && isstruct(OPT_OVERRIDE)
    fn = fieldnames(OPT_OVERRIDE);
    for k = 1:numel(fn), OPT.(fn{k}) = OPT_OVERRIDE.(fn{k}); end
end

if OPT.save_figures && ~exist(OPT.fig_dir, 'dir'), mkdir(OPT.fig_dir); end
P = guardian_params();
P.t_end = OPT.t_end;
if ~isempty(OPT.P_mod), P = OPT.P_mod(P); P = design_compensator(P); end

fprintf('Guardian predictive fire-control simulation\n');
fprintf('Compensator: %s\n', P.Gc.description);
fprintf('Estimator  : EKF Singer model, tau_m = %.1f s, sigma_m = %.1f m/s^2\n', ...
        P.est.tau_m, P.est.sigma_m);
fprintf('Prediction : Tp = %.2f s\n', P.Tp);

%% ===================== A. ANALYTICAL CONTROLLER DESIGN ====================
if OPT.run_design_analysis
    design = analyze_controller_design(P, OPT);
end

%% ===================== B. ENVELOPE VERIFICATION ===========================
case_ids = {};
if OPT.run_standard, case_ids = [case_ids, {'A','B','C'}]; end
if OPT.run_extended, case_ids = [case_ids, {'D','E','F','G','H','I'}]; end

fprintf('\n================ TRAJECTORY ENVELOPE CHECK ================\n');
for k = 1:numel(case_ids)
    sc  = scenario_library(case_ids{k});
    chk = check_trajectory_limits(sc, P);
    fprintf('%-24s %s\n', sc.name, chk.summary);
    if ~chk.pass
        warning('Scenario %s violates the operating envelope.', sc.id);
    end
end

%% ===================== C/D. DETERMINISTIC TEST CASES ======================
results = struct('id', {}, 'name', {}, 'res', {}, 'm', {});
for k = 1:numel(case_ids)
    sc  = scenario_library(case_ids{k});
    res = simulate_engagement(sc, P, OPT);
    m   = compute_metrics(res, P);
    results(end+1) = struct('id', sc.id, 'name', sc.name, 'res', res, 'm', m); %#ok<SAGROW>
    fprintf('Case %-22s: %s\n', sc.name, one_line(m));
    if OPT.plot_each_case
        plot_scenario_results(res, m, P, OPT);
    end
end

if ~isempty(results)
    T = metrics_table(results);
    fprintf('\n=================== SUMMARY - FINAL DESIGN ===================\n');
    disp(T);
    plot_summary(results, P, OPT);
end

%% ===================== E. MONTE-CARLO GENERALIZATION ======================
mc = [];
if OPT.run_monte_carlo
    fprintf('\nMonte-Carlo: %d random admissible engagements...\n', OPT.n_monte_carlo);
    mc = struct('sc', {}, 'chk', {}, 'm', {}, 'pulses', {});
    for k = 1:OPT.n_monte_carlo
        sc  = random_scenario(k, P);
        chk = check_trajectory_limits(sc, P);
        res = simulate_engagement(sc, P, OPT);
        m   = compute_metrics(res, P);
        mc(end+1) = struct('sc', sc, 'chk', chk, 'm', m, 'pulses', res.pulses); %#ok<SAGROW>
        if mod(k, 10) == 0, fprintf('  %d / %d done\n', k, OPT.n_monte_carlo); end
    end
    Mm = [mc.m];
    fprintf('  success rate        : %.1f %% (%d / %d)\n', 100*mean([Mm.success]), ...
            sum([Mm.success]), numel(Mm));
    fprintf('  intercept time      : median %.2f s, 90th pct %.2f s, max %.2f s\n', ...
            median([Mm.t_intercept], 'omitnan'), prctile_local([Mm.t_intercept], 90), ...
            max([Mm.t_intercept]));
    fprintf('  hit ratio of pulses : %.1f %%\n', 100*sum([Mm.n_hits])/max(sum([Mm.n_pulses]),1));
    fprintf('  RMS pointing (mean) : %.2f m\n', mean([Mm.rms_point_m]));
    fail = find(~[Mm.success]);
    if ~isempty(fail), fprintf('  failed seeds: %s\n', mat2str(fail)); end
    plot_monte_carlo(mc, P, OPT);
end

%% ===================== F. ARCHITECTURE COMPARISON =========================
% Demonstrates WHY the delay-matched loop is used: the laser is instantaneous,
% so the beam must point at the CURRENT target. A loop that tracks the
% predicted reference tightly points Tp seconds ahead of the target.
comparison = [];
if OPT.run_comparison && ~isempty(case_ids)
    variants = { ...
        'Final: delay-matched, Tp=0.5', struct('ctrl_mode', 'delay_matched', 'Tp', 0.5); ...
        'Tight PI, Tp=0.5 (reference tracked tightly)', struct('ctrl_mode', 'tight', 'Tp', 0.5); ...
        'Tight PI, no prediction (Tp=0)', struct('ctrl_mode', 'tight', 'Tp', 0.0); ...
        'Delay-matched loop, no prediction', struct('ctrl_mode', 'delay_matched', 'Tp', 0.5, 'Tp_ref', 0)};
    fprintf('\n================ ARCHITECTURE COMPARISON ================\n');
    comparison = struct('name', variants(:,1)', 'success', [], 't_int', [], 'rms_point', [], ...
                        'hits', [], 'pulses', [], 'V_peak', []);
    for v = 1:size(variants, 1)
        Pv = P; s = variants{v, 2};
        Pv.ctrl_mode = s.ctrl_mode; Pv.Tp = s.Tp;
        Pv = design_compensator(Pv);              % Gc designed for this Tp
        if isfield(s, 'Tp_ref'), Pv.Tp = s.Tp_ref; end   % horizon actually used
        for k = 1:numel(case_ids)
            sc = scenario_library(case_ids{k});
            res = simulate_engagement(sc, Pv, OPT);
            m = compute_metrics(res, Pv);
            comparison(v).success(k)   = m.success;
            comparison(v).t_int(k)     = m.t_intercept;
            comparison(v).rms_point(k) = m.rms_point_m;
            comparison(v).hits(k)      = m.n_hits;
            comparison(v).pulses(k)    = m.n_pulses;
            comparison(v).V_peak(k)    = m.V_peak;
        end
        fprintf('%-46s success %d/%d, mean t_int %.2f s, mean RMS pointing %6.2f m, hits %d/%d, peak V %.0f\n', ...
            comparison(v).name, sum(comparison(v).success), numel(case_ids), ...
            mean(comparison(v).t_int, 'omitnan'), mean(comparison(v).rms_point), ...
            sum(comparison(v).hits), sum(comparison(v).pulses), max(comparison(v).V_peak));
    end
    plot_comparison(comparison, case_ids, P, OPT);
end

%% ===================== H. ESTIMATOR SENSITIVITY (TUNING) ==================
% Justifies the Singer parameters. Uses tuning seeds 1001.. that are disjoint
% from the Monte-Carlo evaluation seeds 1..N, to avoid tuning on the test set.
sens = [];
if OPT.run_sensitivity
    fprintf('\n================ ESTIMATOR SENSITIVITY ================\n');
    ns = numel(OPT.sens_sigma_m); nt = numel(OPT.sens_tau_m);
    sens.sigma_m = OPT.sens_sigma_m; sens.tau_m = OPT.sens_tau_m;
    sens.success = zeros(ns, nt); sens.miss = zeros(ns, nt);
    sens.inR = zeros(ns, nt); sens.tint = zeros(ns, nt); sens.pred = zeros(ns, nt);
    tune_sc = cell(1, OPT.n_sens_cases);
    for c = 1:OPT.n_sens_cases, tune_sc{c} = random_scenario(1000 + c, P); end
    for a = 1:ns
        for b = 1:nt
            Ps = P; Ps.est.sigma_m = OPT.sens_sigma_m(a); Ps.est.tau_m = OPT.sens_tau_m(b);
            Ms = [];
            for c = 1:OPT.n_sens_cases
                res = simulate_engagement(tune_sc{c}, Ps, OPT);
                Ms = [Ms compute_metrics(res, Ps)]; %#ok<AGROW>
            end
            sens.success(a,b) = mean([Ms.success]);
            sens.miss(a,b)    = mean([Ms.miss_mean], 'omitnan');
            sens.inR(a,b)     = mean([Ms.frac_in_R]);
            sens.tint(a,b)    = median([Ms.t_intercept], 'omitnan');
            sens.pred(a,b)    = mean([Ms.rms_pred_err]);
            fprintf('sigma_m = %4.1f, tau_m = %3.1f : success %5.1f %%, mean miss %.2f m, in-R %5.1f %%, pred RMS %.2f m\n', ...
                OPT.sens_sigma_m(a), OPT.sens_tau_m(b), 100*sens.success(a,b), sens.miss(a,b), ...
                100*sens.inR(a,b), sens.pred(a,b));
        end
    end
    plot_sensitivity(sens, P, OPT);
end

%% ===================== G. SAVE ============================================
if ~exist(OPT.results_dir, 'dir'), mkdir(OPT.results_dir); end
% Each section is saved to its own file, and only when it ran, so a partial
% run never overwrites results of sections that were switched off.
rd = OPT.results_dir;
if ~isempty(results)
    writetable(T, fullfile(rd, 'summary_metrics.csv'));
    save(fullfile(rd, 'results_cases.mat'), 'P', 'OPT', 'results');      % full logs (for paper figures)
end
if ~isempty(mc),         save(fullfile(rd, 'results_mc.mat'),          'P', 'OPT', 'mc');         end
if ~isempty(comparison), save(fullfile(rd, 'results_comparison.mat'),  'P', 'OPT', 'comparison', 'case_ids'); end
if ~isempty(sens),       save(fullfile(rd, 'results_sensitivity.mat'), 'P', 'OPT', 'sens');       end
fprintf('\nResults saved to %s\n', rd);


%% ========================================================================
%%                          LOCAL FUNCTIONS
%% ========================================================================
function res = simulate_engagement(sc, P, OPT)
%SIMULATE_ENGAGEMENT  The simulation clock for one engagement.

dt = OPT.dt;
t  = 0:dt:OPT.t_end;
N  = numel(t);

turret_controller('reset', P);
target_motion_generator('reset');

x = [OPT.theta0; OPT.theta_dot0; OPT.i0];         % plant state [theta; theta_dot; i]

% ---- preallocate logs
res.sc = sc; res.P = P; res.t = t;
res.theta = zeros(1,N); res.theta_dot = zeros(1,N); res.i = zeros(1,N); res.V = zeros(1,N);
res.r_meas = nan(1,N); res.theta_meas = nan(1,N); res.valid = zeros(1,N);
res.theta_ref = nan(1,N); res.x_cur = nan(6,N); res.x_pred = nan(6,N);
res.d_hat = nan(1,N); res.sigma_perp = nan(1,N); res.fire_cmd = zeros(1,N);
res.pos_std = nan(2,N); res.vel_std = nan(2,N);
res.pulses = struct('t', {}, 'theta', {}, 'target', {}, 'miss', {}, 'miss_signed', {}, 'hit', {});
res.ignored = [];

laser_t_last = -Inf;  fire_prev = 0;

for n = 1:N
    tn = t(n);

    % ---- 1. radar
    [r_m, th_m, valid] = target_motion_generator(tn, sc);

    % ---- 2. fire-control computer
    [V, fire] = turret_controller(tn, r_m, th_m, valid, x(1), x(2));
    dbg = turret_controller('debug');

    % ---- 3. laser model (independent of the controller's own logic)
    if fire == 1 && fire_prev == 0                          % 0 -> 1 transition
        if tn - laser_t_last >= P.t_recharge - 1e-9         % recharged?
            tr = target_truth(tn, sc);
            p  = [tr.x; tr.y];
            [miss, miss_signed] = ray_point_distance(x(1), p);
            res.pulses(end+1) = struct('t', tn, 'theta', x(1), 'target', p, ...
                'miss', miss, 'miss_signed', miss_signed, 'hit', miss <= P.R_eng);
            laser_t_last = tn;
        else
            res.ignored(end+1) = tn;                        % too early: ignored
        end
    end
    fire_prev = fire;

    % ---- 4. log
    res.theta(n) = x(1); res.theta_dot(n) = x(2); res.i(n) = x(3); res.V(n) = V;
    res.r_meas(n) = r_m; res.theta_meas(n) = th_m; res.valid(n) = valid;
    res.fire_cmd(n) = fire;
    if dbg.est_init
        res.theta_ref(n) = dbg.theta_ref;
        res.x_cur(:,n) = dbg.x_cur;  res.x_pred(:,n) = dbg.x_pred;
        res.d_hat(n) = dbg.d_hat;    res.sigma_perp(n) = dbg.sigma_perp;
        res.pos_std(:,n) = dbg.pos_std(:); res.vel_std(:,n) = dbg.vel_std(:);
    end

    % ---- 5. plant integration over [tn, tn+dt], voltage held (ZOH)
    if n < N
        x = rk4_step(@(tt, xx) plant_rhs(tt, xx, V), tn, x, dt);
        % mechanical hard stops
        if x(1) > P.theta_max, x(1) = P.theta_max; x(2) = min(x(2), 0); end
        if x(1) < P.theta_min, x(1) = P.theta_min; x(2) = max(x(2), 0); end
    end
end

% truth on the whole grid and true pointing error
res.truth = target_truth(t, sc);
nb = [-sin(res.theta); cos(res.theta)];
res.point_err_m = nb(1,:).*res.truth.x + nb(2,:).*res.truth.y;  % signed cross-beam distance
end

function dx = plant_rhs(t, x, V)
[th_d, th_dd, i_d] = turret_plant(t, V, x(1), x(2), x(3));
dx = [th_d; th_dd; i_d];
end

function x = rk4_step(f, t, x, h)
k1 = f(t, x);
k2 = f(t + h/2, x + h/2*k1);
k3 = f(t + h/2, x + h/2*k2);
k4 = f(t + h,   x + h*k3);
x  = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

function [d, d_signed] = ray_point_distance(theta, p)
% Minimum distance between the ray from the origin along theta and point p.
u = [cos(theta); sin(theta)];  nrm = [-sin(theta); cos(theta)];
d_signed = nrm.' * p;
if u.' * p >= 0, d = abs(d_signed); else, d = norm(p); end
end

function s = one_line(m)
if m.success
    s = sprintf('HIT at %5.2f s | %d/%d pulses hit | min miss %.2f m | RMS point %.2f m', ...
        m.t_intercept, m.n_hits, m.n_pulses, m.miss_min, m.rms_point_m);
else
    s = sprintf('NO HIT        | %d pulses | min miss %.2f m | RMS point %.2f m', ...
        m.n_pulses, m.miss_min, m.rms_point_m);
end
end

function T = metrics_table(results)
M = [results.m];
T = table(string({results.id}'), [M.success]', [M.t_intercept]', [M.n_hits]', [M.n_pulses]', ...
    [M.miss_min]', [M.miss_mean]', [M.t_acquire]', [M.rms_point_m]', [M.max_point_m]', ...
    100*[M.frac_in_R]', [M.rms_track_deg]', [M.rms_pos_err]', [M.rms_vel_err]', [M.rms_pred_err]', ...
    [M.bo_max_point_m]', [M.bo_recovery]', [M.V_peak]', [M.i_peak]', ...
    'VariableNames', {'Case','Success','t_intercept_s','Hits','Pulses','MissMin_m','MissMean_m', ...
    't_acquire_s','RMSpoint_m','MaxPoint_m','PctInR','RMStrack_deg','RMSposErr_m','RMSvelErr_mps', ...
    'RMSpredErr_m','BOmaxPoint_m','BOrecovery_s','Vpeak_V','Ipeak_A'});
end

function q = prctile_local(x, p)
x = sort(x(~isnan(x)));
if isempty(x), q = NaN; return; end
q = x(max(1, ceil(p/100*numel(x))));
end

function plot_sensitivity(sens, P, OPT)
fig = figure('Name', 'Estimator sensitivity', 'Color', 'w', 'Visible', OPT.fig_visible, ...
             'Position', [160 120 1250 380]);
sgtitle('Estimator (Singer) parameter sensitivity on tuning cases', 'FontWeight', 'bold');
data  = {100*sens.success, sens.miss, 100*sens.inR, sens.pred};
ttl   = {'success [%]', 'mean pulse miss [m]', 'time beam in R_{eng} [%]', 'RMS prediction error [m]'};
for k = 1:4
    subplot(1,4,k);
    imagesc(1:numel(sens.tau_m), 1:numel(sens.sigma_m), data{k}); colorbar;
    set(gca, 'XTick', 1:numel(sens.tau_m), 'XTickLabel', sens.tau_m, ...
             'YTick', 1:numel(sens.sigma_m), 'YTickLabel', sens.sigma_m, 'YDir', 'normal');
    xlabel('\tau_m [s]'); ylabel('\sigma_m [m/s^2]'); title(ttl{k});
    for a = 1:numel(sens.sigma_m)
        for b = 1:numel(sens.tau_m)
            text(b, a, sprintf('%.2g', data{k}(a,b)), 'HorizontalAlignment', 'center', ...
                 'Color', 'w', 'FontWeight', 'bold');
        end
    end
    hold on;
    [~, ia] = min(abs(sens.sigma_m - P.est.sigma_m)); [~, ib] = min(abs(sens.tau_m - P.est.tau_m));
    rectangle('Position', [ib-0.5 ia-0.5 1 1], 'EdgeColor', 'r', 'LineWidth', 2);
end
if OPT.save_figures, save_fig(fig, fullfile(OPT.fig_dir, 'estimator_sensitivity')); end
end

function plot_comparison(cmp, case_ids, P, OPT)
fig = figure('Name', 'Architecture comparison', 'Color', 'w', 'Visible', OPT.fig_visible, ...
             'Position', [140 100 1250 450]);
sgtitle('Architecture comparison on deterministic cases', 'FontWeight', 'bold');
x = 1:numel(case_ids);
names = {cmp.name};
subplot(1,3,1);
ti = reshape([cmp.t_int], numel(case_ids), []); ti(isnan(ti)) = P.t_end;
bar(x, ti); set(gca, 'XTick', x, 'XTickLabel', case_ids); grid on; ylim([0 P.t_end]);
ylabel('[s]'); title('First interception (10 s = none)');
subplot(1,3,2);
rp = reshape([cmp.rms_point], numel(case_ids), []);
bar(x, rp); set(gca, 'XTick', x, 'XTickLabel', case_ids, 'YScale', 'log'); grid on;
yline(P.R_eng, 'r--'); ylabel('[m]'); title('RMS pointing error (t \geq 3 s)');
subplot(1,3,3);
h = reshape([cmp.hits], numel(case_ids), []);
bar(x, h); set(gca, 'XTick', x, 'XTickLabel', case_ids); grid on;
ylabel('hits'); title('Successful pulses');
legend(names, 'Location', 'southoutside', 'FontSize', 7);
if OPT.save_figures, save_fig(fig, fullfile(OPT.fig_dir, 'architecture_comparison')); end
end
