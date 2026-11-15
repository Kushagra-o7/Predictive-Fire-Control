function P = guardian_params()
%GUARDIAN_PARAMS  Single source of truth for every physical and design parameter.
%
%   P = guardian_params() returns a struct used by the plant, the controller,
%   the target generator and the analysis scripts. Physical values are taken
%   directly from Section 2 of the ME319 Assignment 2 statement. Design values
%   (estimator, compensator, fire logic) are chosen once and are NOT retuned
%   per scenario.

%% ---------------- Physical parameters: elevation drive (Section 2) -------
P.J_turret = 1200;      % turret + cradle inertia about trunnion   [kg m^2]
P.J_m      = 0.20;      % motor rotor inertia                      [kg m^2]
P.N        = 100;       % gear ratio (motor speed / turret speed)  [-]
P.b        = 1500;      % turret-side viscous damping              [N m s/rad]
P.Kt       = 10;        % motor torque constant                    [N m/A]
P.Ke       = 10;        % motor back-emf constant                  [V s/rad]
P.Rm       = 1.0;       % armature resistance                      [Ohm]
P.Lm       = 0.10;      % armature inductance                      [H]

% Motor inertia reflected to the turret side through the N:1 reduction
P.J_eq     = P.J_turret + P.N^2 * P.J_m;          % = 3200 kg m^2

% Mechanical elevation limits
P.theta_min = deg2rad(-10);                        % [rad]
P.theta_max = deg2rad( 85);                        % [rad]

%% ---------------- Plant transfer function --------------------------------
%   Gp(s) = theta(s)/V(s) = K_N / ( s * (d2 s^2 + d1 s + d0) )
%   with K_N = N Kt,  d2 = Lm J_eq,  d1 = Lm b + Rm J_eq,  d0 = Rm b + N^2 Kt Ke
P.K_N = P.N * P.Kt;                                % 1000
P.d2  = P.Lm * P.J_eq;                             % 320
P.d1  = P.Lm * P.b + P.Rm * P.J_eq;                % 3350
P.d0  = P.Rm * P.b + P.N^2 * P.Kt * P.Ke;          % 1 001 500
P.Gp.num = P.K_N;
P.Gp.den = [P.d2 P.d1 P.d0 0];

%% ---------------- Radar (Section 2) --------------------------------------
P.Ts_radar    = 0.05;            % 20 Hz                                 [s]
P.sigma_r     = 1.0;             % range noise std                       [m]
P.sigma_theta = deg2rad(0.02);   % elevation noise std (0.02 deg)        [rad]
P.t_bo_max    = 0.50;            % max measurement-loss duration         [s]
P.t_boo_min   = 1.00;            % earliest measurement-loss onset       [s]

%% ---------------- Laser (Section 2) --------------------------------------
P.t_recharge = 1.0;              % minimum time between pulses           [s]
P.R_eng      = 2.0;              % engagement radius                     [m]

%% ---------------- Target operating envelope (Section 3) ------------------
P.r_min = 800;   P.r_max = 2000; % [m]
P.v_max = 150;                   % [m/s]
P.a_max = 20;                    % [m/s^2]
P.t_end = 10;                    % standardized engagement duration      [s]

%% ---------------- Prediction (Section 4) ---------------------------------
P.Tp = 0.50;                     % prescribed prediction horizon         [s]

%% ---------------- Target-state estimator (design choice) -----------------
% Extended Kalman filter, Cartesian Singer (exponentially-correlated
% acceleration) model per axis; polar (r, theta) radar measurement.
P.est.tau_m    = 1.0;    % maneuver (acceleration) correlation time      [s] (sensitivity study)
P.est.sigma_m  = 8.0;    % RMS maneuver acceleration                     [m/s^2] (sensitivity study)
P.est.sigma_v0 = 100;    % initial velocity std (unknown at first fix)   [m/s]
P.est.sigma_a0 = 10;     % initial acceleration std                      [m/s^2]

%% ---------------- Feedback compensator (design choice) -------------------
% 'delay_matched' : FINAL DESIGN. Gc(s) = Kc/(tau_c s + 1), chosen so that the
%                   closed loop T(s) matches exp(-s*Tp) to 2nd order, i.e. the
%                   loop lag cancels the prediction lead and the beam points
%                   at the CURRENT target position (see design_compensator.m).
% 'tight'         : comparison only. PI + roll-off, type-2 loop that tracks the
%                   predicted reference tightly (beam leads the target).
P.ctrl_mode = 'delay_matched';
P.tight.wc  = 8;         % crossover of the comparison design            [rad/s]
P.V_max     = Inf;       % motor voltage limit (none specified)          [V]

%% ---------------- Fire-control logic (design choice) ---------------------
% Fire only if (a) the laser is recharged and (b) the beam passes within the
% engagement radius of the estimated current target position with margin:
%       d_hat + k_sigma * sigma_perp <= R_eng
P.fire.k_sigma       = 2.0;   % confidence multiplier on cross-beam std    [-]
P.fire.sigma_max     = 0.75;  % do not fire if cross-beam std exceeds this [m]
P.fire.min_track_age = 0.25;  % do not fire before the track has settled   [s]
P.fire.recharge_tol  = 1e-6;  % numerical guard on the recharge check      [s]

%% ---------------- Derived compensator ------------------------------------
P = design_compensator(P);
end
