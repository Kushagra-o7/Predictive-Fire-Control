function [theta_dot_turret, theta_ddot_turret, i_dot] = turret_plant(t, V_m, theta_turret, theta_dot_turret, i)
%TURRET_PLANT  First-principles elevation-drive dynamics (armature-controlled
%   DC motor + 100:1 gear reduction + turret), referred to the turret side.
%
%   [theta_dot_turret, theta_ddot_turret, i_dot] = ...
%        turret_plant(t, V_m, theta_turret, theta_dot_turret, i)
%
%   Electrical (motor):   Lm di/dt = V_m - Rm i - Ke * w_m,     w_m = N * theta_dot
%   Mechanical (turret):  J_eq theta_ddot = N Kt i - b theta_dot
%                         J_eq = J_turret + N^2 J_m   (rotor inertia reflected)
%
%   Gearbox assumed ideal (rigid, lossless, no backlash); turret-side viscous
%   damping only; gravity unbalance neglected (trunnion assumed balanced).
%   The model is linear, so the "small-perturbation" model is exact away from
%   the mechanical stops:
%       Gp(s) = N Kt / ( s [ (Lm s + Rm)(J_eq s + b) + N^2 Kt Ke ] )
%
%   Mechanical hard stops (-10 deg ... +85 deg): when the turret is on a stop
%   and the net torque pushes it further into the stop, the stop reaction
%   cancels the motion (velocity and acceleration outputs set to zero).
%   The integrator in the master file additionally clamps the angle.

persistent P
if isempty(P)
    P = guardian_params();
end

% Electrical dynamics
w_m   = P.N * theta_dot_turret;                       % motor speed  [rad/s]
i_dot = (V_m - P.Rm * i - P.Ke * w_m) / P.Lm;         % [A/s]

% Mechanical dynamics (turret side)
tau_motor_at_turret = P.N * P.Kt * i;                 % [N m]
theta_ddot_turret   = (tau_motor_at_turret - P.b * theta_dot_turret) / P.J_eq;

% Mechanical hard stops
at_upper = theta_turret >= P.theta_max && (theta_dot_turret > 0 || ...
           (theta_dot_turret == 0 && theta_ddot_turret > 0));
at_lower = theta_turret <= P.theta_min && (theta_dot_turret < 0 || ...
           (theta_dot_turret == 0 && theta_ddot_turret < 0));
if at_upper || at_lower
    theta_dot_turret  = 0;
    theta_ddot_turret = 0;
end
end
