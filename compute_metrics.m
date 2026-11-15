function m = compute_metrics(res, P, t_settle)
%COMPUTE_METRICS  Quantitative performance measures for one engagement.
%
%   m = compute_metrics(res, P) uses the log produced by the master simulation.
%   t_settle (default 3 s) excludes the initial slew from steady-state measures.
%
%   Engagement   : success, interception time, pulses, hits, miss distances
%   Pointing     : beam-to-target cross distance (m) and angle error (deg)
%   Tracking     : theta_ref - theta_turret (deliberately ~ Tp*theta_dot for the
%                  delay-matched design)
%   Estimation   : position / velocity estimation errors
%   Prediction   : predicted position (made at t) vs true position at t + Tp
%   Measurement loss : worst pointing error and recovery after the blackout
%   Effort       : peak voltage, current and turret rate

if nargin < 3, t_settle = 3.0; end
t = res.t;  tr = res.truth;

%% Engagement
pl = res.pulses;
m.n_pulses  = numel(pl);
m.n_hits    = sum([pl.hit]);
m.success   = m.n_hits > 0;
m.hit_rate  = m.n_hits / max(m.n_pulses, 1);
m.n_ignored = numel(res.ignored);
if m.success
    m.t_intercept = pl(find([pl.hit], 1)).t;
else
    m.t_intercept = NaN;
end
if m.n_pulses > 0
    miss = [pl.miss];
    m.miss_mean = mean(miss);  m.miss_min = min(miss);  m.miss_max = max(miss);
else
    m.miss_mean = NaN; m.miss_min = NaN; m.miss_max = NaN;
end

%% Pointing (truth-based)
ep = res.point_err_m;                               % signed cross-beam distance
m.t_acquire = first_time(t, abs(ep) <= P.R_eng);
ss = t >= t_settle;
m.rms_point_m   = rms_(ep(ss));
m.max_point_m   = max(abs(ep(ss)));
m.frac_in_R     = mean(abs(ep(ss)) <= P.R_eng);
m.rms_point_deg = rms_(rad2deg(res.theta - tr.theta));

%% Reference tracking
m.rms_track_deg = rms_(rad2deg(res.theta_ref(ss) - res.theta(ss)));

%% Estimation (after the first second)
ok = t >= 1.0 & ~isnan(res.x_cur(1,:));
ex = res.x_cur(1,ok) - tr.x(ok);   ey = res.x_cur(4,ok) - tr.y(ok);
evx = res.x_cur(2,ok) - tr.vx(ok); evy = res.x_cur(5,ok) - tr.vy(ok);
m.rms_pos_err = rms_(hypot(ex, ey));
m.rms_vel_err = rms_(hypot(evx, evy));

%% Prediction: x_pred(t) should equal truth(t + Tp)
tp_ok = t >= 1.0 & t + P.Tp <= t(end) & ~isnan(res.x_pred(1,:));
trp = target_truth(t(tp_ok) + P.Tp, res.sc);
m.rms_pred_err = rms_(hypot(res.x_pred(1,tp_ok) - trp.x, res.x_pred(4,tp_ok) - trp.y));
m.rms_pred_err_deg = rms_(rad2deg(atan2(res.x_pred(4,tp_ok), res.x_pred(1,tp_ok)) - trp.theta));

%% Measurement loss
if ~isempty(res.sc.blackout)
    t0 = res.sc.blackout(1);  t1 = t0 + res.sc.blackout(2);
    win = t >= t0 & t <= t1 + 1.0;
    m.bo_max_point_m = max(abs(ep(win)));
    [~, k1] = min(abs(t - t1));
    m.bo_pos_err_end = hypot(res.x_cur(1,k1) - tr.x(k1), res.x_cur(4,k1) - tr.y(k1));
    after = t >= t1;
    m.bo_recovery = first_time(t(after), abs(ep(after)) <= P.R_eng) - t1;
else
    m.bo_max_point_m = NaN; m.bo_pos_err_end = NaN; m.bo_recovery = NaN;
end

%% Effort
m.V_peak   = max(abs(res.V));
m.i_peak   = max(abs(res.i));
m.rate_peak_dps = rad2deg(max(abs(res.theta_dot)));
m.hit_limit = any(res.theta >= P.theta_max - 1e-6 | res.theta <= P.theta_min + 1e-6);
end

function v = rms_(x)
x = x(~isnan(x));
if isempty(x), v = NaN; else, v = sqrt(mean(x.^2)); end
end

function tf = first_time(t, mask)
k = find(mask, 1);
if isempty(k), tf = NaN; else, tf = t(k); end
end
