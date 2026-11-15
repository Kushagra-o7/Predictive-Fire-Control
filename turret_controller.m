function [V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret)
%TURRET_CONTROLLER  Guardian predictive-feedback fire-control computer.
%
%   [V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, ...
%                                       theta_turret, theta_dot_turret)
%
%   Architecture (Fig. 2 of the assignment):
%     1. Target-state estimation : EKF, Cartesian Singer model per axis,
%                                  state [x vx ax y vy ay], polar measurement.
%                                  During radar loss: open-loop propagation.
%     2. Target prediction       : propagate the estimate to t + Tp (Tp = 0.5 s).
%     3. Elevation reference     : theta_ref = atan2(y_pred, x_pred), clipped to
%                                  the mechanical range [-10, 85] deg.
%     4. Feedback controller     : LTI compensator Gc(s) (design_compensator.m),
%                                  implemented by exact ZOH discretization at the
%                                  actual call interval: V_m = Gc * (theta_ref - theta).
%     5. Fire control            : a 0->1 pulse on fire_cmd ONLY when
%                                  (a) the laser has recharged (>= 1 s since the
%                                      last pulse), AND
%                                  (b) the beam passes within the engagement
%                                      radius of the estimated CURRENT target
%                                      position, with a confidence margin:
%                                      d_hat + k_sigma*sigma_perp <= R_eng.
%
%   Persistent state resets automatically on the first call or whenever time
%   runs backwards (a new run). Additional calls used by the master script:
%     turret_controller('reset')        reset with default parameters
%     turret_controller('reset', P)     reset with a parameter struct P
%     dbg = turret_controller('debug')  internal signals of the last call

persistent S

%% ---------------- Special commands ---------------------------------------
if ischar(t) || isstring(t)
    switch lower(char(t))
        case 'reset'
            if nargin >= 2 && isstruct(r_meas), S = init_state(r_meas);
            else,                                S = init_state(guardian_params());
            end
            V_m = 0; fire_cmd = 0;
        case 'debug'
            if isempty(S), V_m = struct(); else, V_m = S.dbg; end
            fire_cmd = [];
        otherwise
            error('turret_controller:cmd', 'Unknown command "%s".', t);
    end
    return
end

if isempty(S)
    S = init_state(guardian_params());
elseif ~isempty(S.t_prev) && t < S.t_prev - 1e-12
    S = init_state(S.P);                      % new run: keep parameter set
end
P = S.P;

%% ---------------- 0. Timing ----------------------------------------------
if isempty(S.t_prev), dt = 0; else, dt = t - S.t_prev; end
S.t_prev = t;

%% ---------------- 4a. Compensator state update (ZOH on previous error) ---
if dt > 0 && ~isempty(S.Ac)
    if abs(dt - S.dt_cache) > 1e-12
        [S.Ad, S.Bd] = c2d_zoh(S.Ac, S.Bc, dt);
        S.dt_cache = dt;
    end
    S.xc = S.Ad * S.xc + S.Bd * S.e_prev;
end

%% ---------------- 1. Target-state estimation (EKF) -----------------------
% Guard: an identical (r, theta) pair repeated with meas_valid = 1 is a held
% sample, not a new measurement - fusing it again would make the EKF overconfident.
is_new = meas_valid && isfinite(r_meas) && isfinite(theta_meas) && ...
         ~(r_meas == S.z_last(1) && theta_meas == S.z_last(2));
if is_new
    z = [r_meas; theta_meas];
    S.z_last = z;
    if ~S.est_init
        % Initialize from the first fix: position from the polar measurement,
        % velocity / acceleration unknown (large covariance).
        c = cos(theta_meas); s = sin(theta_meas);
        Jp = [c, -r_meas*s; s, r_meas*c];
        Ppos = Jp * S.R * Jp.';
        S.x = [r_meas*c; 0; 0; r_meas*s; 0; 0];
        S.Pk = diag([0 P.est.sigma_v0^2 P.est.sigma_a0^2 0 P.est.sigma_v0^2 P.est.sigma_a0^2]);
        S.Pk([1 4],[1 4]) = Ppos;
        S.est_init = true;
        S.t_init   = t;
    else
        [F, Q] = singer_FQ(t - S.t_est, P.est.tau_m, P.est.sigma_m);
        xm = F * S.x;
        Pm = F * S.Pk * F.' + Q;
        % Polar measurement model h(x) = [r; theta]
        px = xm(1); py = xm(4); rr = hypot(px, py);
        h = [rr; atan2(py, px)];
        H = zeros(2, 6);
        H(1,1) =  px/rr;     H(1,4) = py/rr;
        H(2,1) = -py/rr^2;   H(2,4) = px/rr^2;
        nu = z - h;
        nu(2) = atan2(sin(nu(2)), cos(nu(2)));          % wrap angle residual
        Sk = H * Pm * H.' + S.R;
        K  = Pm * H.' / Sk;
        S.x  = xm + K * nu;
        IKH  = eye(6) - K * H;
        S.Pk = IKH * Pm * IKH.' + K * S.R * K.';        % Joseph form
        S.nis = nu.' / Sk * nu;
    end
    S.t_est = t;
    S.t_last_meas = t;
end

%% ---------------- 2. Prediction & 3. elevation reference ----------------
theta_ref = theta_turret;  e = 0;
x_cur = nan(6,1); x_pred = nan(6,1); Pc = nan(6);
if S.est_init
    dcur = t - S.t_est;                          % time since last update
    [Fc, Qc] = singer_FQ(dcur, P.est.tau_m, P.est.sigma_m);
    x_cur = Fc * S.x;                            % estimate at time t
    Pc    = Fc * S.Pk * Fc.' + Qc;
    Fp    = singer_FQ(dcur + P.Tp, P.est.tau_m, P.est.sigma_m);
    x_pred = Fp * S.x;                           % prediction at t + Tp
    theta_ref = atan2(x_pred(4), x_pred(1));
    theta_ref = min(max(theta_ref, P.theta_min), P.theta_max);
    e = theta_ref - theta_turret;
end

%% ---------------- 4b. Compensator output ---------------------------------
if isempty(S.Ac), V_m = S.Dc * e; else, V_m = S.Cc * S.xc + S.Dc * e; end
V_m = min(max(V_m, -P.V_max), P.V_max);
S.e_prev = e;

%% ---------------- 5. Fire control ----------------------------------------
fire_cmd = 0;
d_hat = NaN; sigma_perp = NaN; in_radius = false;
recharged = (t - S.t_last_fire) >= P.t_recharge + P.fire.recharge_tol;
if S.est_init
    u = [cos(theta_turret); sin(theta_turret)];   % beam direction
    n = [-sin(theta_turret); cos(theta_turret)];  % beam normal
    p = [x_cur(1); x_cur(4)];
    if u.' * p > 0, d_hat = abs(n.' * p); else, d_hat = norm(p); end
    sigma_perp = sqrt(max(n.' * Pc([1 4],[1 4]) * n, 0));
    track_ok   = (t - S.t_init) >= P.fire.min_track_age && sigma_perp <= P.fire.sigma_max;
    in_radius  = track_ok && (d_hat + P.fire.k_sigma * sigma_perp <= P.R_eng);
    if in_radius && recharged && S.fire_prev == 0
        fire_cmd = 1;                             % 0 -> 1 transition = one pulse
        S.t_last_fire = t;
    end
end
S.fire_prev = fire_cmd;                           % back to 0 on the next call

%% ---------------- Debug / logging ----------------------------------------
S.dbg = struct('t', t, 'est_init', S.est_init, 'x_cur', x_cur, 'x_pred', x_pred, ...
    'theta_ref', theta_ref, 'error', e, 'V_m', V_m, 'd_hat', d_hat, ...
    'sigma_perp', sigma_perp, 'in_radius', in_radius, 'recharged', recharged, ...
    'fire_cmd', fire_cmd, 'pos_std', sqrt(max([Pc(1,1) Pc(4,4)], 0)), ...
    'vel_std', sqrt(max([Pc(2,2) Pc(5,5)], 0)), 'time_since_meas', t - S.t_last_meas);
end

%% ========================================================================
function S = init_state(P)
S.P = P;
S.t_prev = [];
S.R = diag([P.sigma_r^2, P.sigma_theta^2]);
% Estimator
S.est_init = false; S.x = zeros(6,1); S.Pk = eye(6);
S.t_est = 0; S.t_init = 0; S.t_last_meas = -Inf; S.nis = NaN; S.z_last = [NaN; NaN];
% Compensator (continuous state space, discretized at run time)
S.Ac = P.Gc.A; S.Bc = P.Gc.B; S.Cc = P.Gc.C; S.Dc = P.Gc.D;
S.xc = zeros(size(S.Ac, 1), 1);
S.Ad = []; S.Bd = []; S.dt_cache = -1;
S.e_prev = 0;
% Fire control
S.t_last_fire = -Inf; S.fire_prev = 0;
S.dbg = struct();
end

function [F, Q] = singer_FQ(dt, tau_m, sigma_m)
% Exact discretization of the Singer model (per axis: [p v a],
% a_dot = -a/tau_m + w, E[w w'] = 2 sigma_m^2/tau_m), assembled for [x; y].
alpha = 1/tau_m;
Ac = [0 1 0; 0 0 1; 0 0 -alpha];
if nargout < 2
    F1 = expm(Ac * dt);
    F = blkdiag(F1, F1);
    return
end
qc = 2 * alpha * sigma_m^2;
G  = [0; 0; 1];
M  = [-Ac, G*qc*G.'; zeros(3), Ac.'] * dt;          % Van Loan
E  = expm(M);
F1 = E(4:6, 4:6).';
Q1 = F1 * E(1:3, 4:6);
Q1 = (Q1 + Q1.') / 2;
F = blkdiag(F1, F1);
Q = blkdiag(Q1, Q1);
end

function [Ad, Bd] = c2d_zoh(A, B, dt)
n = size(A, 1); m = size(B, 2);
E = expm([A B; zeros(m, n + m)] * dt);
Ad = E(1:n, 1:n);
Bd = E(1:n, n+1:n+m);
end
