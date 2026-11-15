function P = design_compensator(P)
%DESIGN_COMPENSATOR  Analytical design of the LTI feedback compensator Gc(s).
%
%   P = design_compensator(P) fills P.Gc.num / P.Gc.den (descending powers of s)
%   according to P.ctrl_mode. Call it again whenever P.Tp or P.ctrl_mode changes.
%
%   ---------------------------------------------------------------------------
%   FINAL DESIGN: 'delay_matched'
%   ---------------------------------------------------------------------------
%   Plant:  Gp(s) = K_N / ( s D(s) ),   D(s) = d2 s^2 + d1 s + d0
%   Comp.:  Gc(s) = Kc / (tau_c s + 1)
%
%   The laser has negligible time of flight, so a hit at time t requires
%   theta_turret(t) = theta_target(t). The reference is the PREDICTED angle
%   theta_ref(t) ~ theta_target(t + Tp), i.e. theta_ref = exp(+s Tp) theta.
%   Pointing error:
%       E(s) = [1 - T(s) exp(s Tp)] theta_target(s),   T = L/(1+L), L = Gc Gp.
%   Expanding 1/T(s) = 1 + s(tau_c s + 1) D(s) / (K_N Kc) in powers of s:
%       1/T = 1 + [d0/(K_N Kc)] s + [(d1 + tau_c d0)/(K_N Kc)] s^2 + ...
%   and requiring 1/T = exp(s Tp) = 1 + Tp s + Tp^2 s^2/2 + ... gives
%       Kc    = d0 / (K_N Tp)                (velocity constant Kv = 1/Tp)
%       tau_c = Tp/2 - d1/d0
%   With these two choices the pointing error is zero for constant AND linearly
%   varying target angular rate (constant theta-dot and theta-ddot); the residual
%   error is proportional to the third derivative of the target angle.
%   For Tp = 0.5 s:  Kc = 2003 V/rad, tau_c = 0.2467 s; dominant closed-loop
%   poles ~ Butterworth, wn ~ sqrt(2)/Tp = 2.83 rad/s, zeta ~ 0.707.
%
%   ---------------------------------------------------------------------------
%   COMPARISON DESIGN: 'tight'  (used only in the comparison study)
%   ---------------------------------------------------------------------------
%   Gc(s) = Kp (s + wi)/s * 1/(s/wf + 1),  Kp = wc d0/K_N, wi = wc/4, wf = 3 wc.
%   Type-2 loop: tracks the predicted reference with zero ramp error.

switch lower(P.ctrl_mode)
    case 'delay_matched'
        Kc    = P.d0 / (P.K_N * P.Tp);
        tau_c = P.Tp/2 - P.d1/P.d0;
        if tau_c <= 0
            error('design_compensator:tau', 'Tp too small for delay matching.');
        end
        P.Gc.num = Kc;
        P.Gc.den = [tau_c 1];
        P.Gc.Kc = Kc;  P.Gc.tau_c = tau_c;
        P.Gc.description = sprintf('Gc(s) = %.1f / (%.4f s + 1)', Kc, tau_c);

    case 'tight'
        wc = P.tight.wc;
        Kp = wc * P.d0 / P.K_N;
        wi = wc/4;  wf = 3*wc;
        P.Gc.num = Kp * [1 wi];
        P.Gc.den = conv([1 0], [1/wf 1]);
        P.Gc.Kp = Kp; P.Gc.wi = wi; P.Gc.wf = wf;
        P.Gc.description = sprintf('Gc(s) = %.0f (s+%.2f)/s / (s/%.1f + 1)', Kp, wi, wf);

    otherwise
        error('design_compensator:mode', 'Unknown ctrl_mode "%s".', P.ctrl_mode);
end

% State-space realization (controllable canonical form) used by the
% digital implementation in turret_controller.m
[P.Gc.A, P.Gc.B, P.Gc.C, P.Gc.D] = tf2ss_local(P.Gc.num, P.Gc.den);
end

function [A, B, C, D] = tf2ss_local(num, den)
% Controllable canonical realization of a proper SISO transfer function
% (toolbox-free so the code also runs in Octave / without Control Toolbox).
num = num(:).'; den = den(:).';
num = num / den(1); den = den / den(1);
n = numel(den) - 1;
if numel(num) > n + 1
    error('tf2ss_local:improper', 'Compensator must be proper.');
end
num = [zeros(1, n + 1 - numel(num)) num];
D = num(1);
if n == 0
    A = zeros(0); B = zeros(0, 1); C = zeros(1, 0); return
end
r = num(2:end) - D * den(2:end);          % strictly proper remainder
A = [zeros(n-1, 1) eye(n-1); -fliplr(den(2:end))];
B = [zeros(n-1, 1); 1];
C = fliplr(r);
end
