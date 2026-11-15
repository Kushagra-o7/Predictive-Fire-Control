function tr = target_truth(t, sc)
%TARGET_TRUTH  Exact (noise-free) target kinematics for scenario struct sc.
%
%   tr = target_truth(t, sc) with t scalar or vector returns a struct with row
%   vectors x, y, vx, vy, ax, ay, r, theta, speed, accel. Used by the
%   measurement generator, the laser hit test and the performance analysis.

t = t(:).';
[tr.x, tr.vx, tr.ax] = axis_motion(t, sc.x0, sc.vx0, sc.man.x);
[tr.y, tr.vy, tr.ay] = axis_motion(t, sc.y0, sc.vy0, sc.man.y);
tr.r     = hypot(tr.x, tr.y);
tr.theta = atan2(tr.y, tr.x);
tr.speed = hypot(tr.vx, tr.vy);
tr.accel = hypot(tr.ax, tr.ay);
end

function [p, v, a] = axis_motion(t, p0, v0, comps)
p = p0 + v0 * t;
v = v0 * ones(size(t));
a = zeros(size(t));
for k = 1:size(comps, 1)
    A = comps(k,1); w = comps(k,2); nper = comps(k,3); t0 = comps(k,4);
    if A == 0 || w == 0, continue; end
    T     = nper * 2*pi / w;                 % maneuver duration (Inf allowed)
    tau   = t - t0;
    tau_c = min(max(tau, 0), T);             % clipped maneuver time
    dv = (A/w) * (1 - cos(w*tau_c));
    dx = (A/w) * (tau_c - sin(w*tau_c)/w) + dv .* max(tau - T, 0);
    active = (tau > 0) & (tau < T);
    p = p + dx;
    v = v + dv;
    a = a + A * sin(w*tau_c) .* active;
end
end
