function chk = check_trajectory_limits(sc, P, t_end)
%CHECK_TRAJECTORY_LIMITS  Verify a scenario against the operating envelope.
%
%   chk = check_trajectory_limits(sc, P) evaluates the true trajectory on a
%   1 ms grid over [0, t_end] and checks (Section 2 & 3):
%       800 m <= r(t) <= 2000 m,  |v| <= 150 m/s,  |a| <= 20 m/s^2,
%       measurement-loss onset >= 1.00 s and duration <= 0.50 s.

if nargin < 3, t_end = P.t_end; end
t  = 0:1e-3:t_end;
tr = target_truth(t, sc);

chk.r_min = min(tr.r);      chk.r_max = max(tr.r);
chk.v_max = max(tr.speed);  chk.a_max = max(tr.accel);
chk.theta_min_deg = rad2deg(min(tr.theta));
chk.theta_max_deg = rad2deg(max(tr.theta));

tol = 1e-9;
chk.range_ok = chk.r_min >= P.r_min - tol && chk.r_max <= P.r_max + tol;
chk.speed_ok = chk.v_max <= P.v_max + tol;
chk.accel_ok = chk.a_max <= P.a_max + tol;
if isempty(sc.blackout)
    chk.blackout_ok = true;
else
    chk.blackout_ok = sc.blackout(1) >= P.t_boo_min - tol && sc.blackout(2) <= P.t_bo_max + tol ...
                      && sc.blackout(2) >= 0;
end
chk.pass = chk.range_ok && chk.speed_ok && chk.accel_ok && chk.blackout_ok;
chk.summary = sprintf(['r = [%6.1f, %6.1f] m  |v|max = %5.1f m/s  |a|max = %4.1f m/s^2  ' ...
    'theta = [%5.1f, %5.1f] deg  -> %s'], chk.r_min, chk.r_max, chk.v_max, chk.a_max, ...
    chk.theta_min_deg, chk.theta_max_deg, ternary(chk.pass, 'PASS', 'FAIL'));
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
