function sc = scenario_library(id)
%SCENARIO_LIBRARY  Definitions of all deterministic test cases.
%
%   sc = scenario_library(id)   returns the scenario struct for id ('A'...'I')
%   ids = scenario_library('list') returns the list of all ids.
%
%   Geometry: turret/radar at the origin, x = horizontal ground range toward
%   the target [m], y = altitude [m]; elevation theta = atan2(y, x).
%
%   Target motion = constant-velocity baseline + smooth bounded maneuvers.
%   Each maneuver component (one row of man.x or man.y) is
%       [A  w  nper  t0]  ->  a(t) = A sin(w (t - t0))  for 0 <= t - t0 <= nper*2*pi/w
%   (nper = Inf for a continuous weave; integer/half-integer nper keeps the
%   acceleration continuous). Position/velocity are integrated analytically in
%   target_truth.m. Every case is checked by check_trajectory_limits.m.
%
%   blackout = [t_boo  t_bo]  (onset, duration) or [] for none.
%
%   A, B, C : standardized scenarios of Section 3.
%   D ... I : additional stress / evaluation-style cases (different maneuvers
%             and measurement-loss conditions) used to test generalization.

if nargin < 1, id = 'list'; end
if ischar(id) && strcmpi(id, 'list')
    sc = {'A','B','C','D','E','F','G','H','I'};
    return
end
if isnumeric(id)
    ids = scenario_library('list'); id = ids{id};
end

base = struct('id', upper(id), 'name', '', 'description', '', ...
    'x0', 0, 'y0', 0, 'vx0', 0, 'vy0', 0, ...
    'man', struct('x', zeros(0,4), 'y', zeros(0,4)), ...
    'blackout', [], 'noise_seed', 0);
sc = base;

switch upper(id)
    case 'A'
        sc.name = 'A: Nominal';
        sc.description = 'Non-maneuvering, constant-velocity level inbound flight.';
        sc.x0 = 1600; sc.y0 = 900; sc.vx0 = -100; sc.vy0 = 0;
        sc.noise_seed = 101;

    case 'B'
        sc.name = 'B: Maneuvering';
        sc.description = 'Smooth bounded weave in x and y beginning at t = 1 s.';
        sc.x0 = 1600; sc.y0 = 900; sc.vx0 = -100; sc.vy0 = 0;
        sc.man.x = [  8  1.2  Inf  1.0];
        sc.man.y = [ 12  1.6  Inf  1.0];
        sc.noise_seed = 102;

    case 'C'
        sc.name = 'C: Measurement loss';
        sc.description = 'Nominal trajectory; radar unavailable 1.00-1.50 s.';
        sc.x0 = 1600; sc.y0 = 900; sc.vx0 = -100; sc.vy0 = 0;
        sc.blackout = [1.00 0.50];
        sc.noise_seed = 103;

    case 'D'
        sc.name = 'D: Maneuver + loss';
        sc.description = 'Scenario B maneuver with blackout exactly at maneuver onset.';
        sc.x0 = 1600; sc.y0 = 900; sc.vx0 = -100; sc.vy0 = 0;
        sc.man.x = [  8  1.2  Inf  1.0];
        sc.man.y = [ 12  1.6  Inf  1.0];
        sc.blackout = [1.00 0.50];
        sc.noise_seed = 104;

    case 'E'
        sc.name = 'E: High-g jinks';
        sc.description = 'Near a_max: two-period vertical jink + horizontal half-period break.';
        sc.x0 = 1500; sc.y0 = 700; sc.vx0 = -80; sc.vy0 = 10;
        sc.man.x = [ -9  1.0  0.5  1.0];
        sc.man.y = [ 17  2.5  2.0  1.5];
        sc.noise_seed = 105;

    case 'F'
        sc.name = 'F: Climbing crosser';
        sc.description = 'Fast climb (high elevation rate), weave, blackout at 3.0 s.';
        sc.x0 = 1400; sc.y0 = 250; sc.vx0 = -60; sc.vy0 = 85;
        sc.man.x = [  6  1.8  Inf  1.0];
        sc.man.y = [-10  1.3  Inf  2.0];
        sc.blackout = [3.00 0.50];
        sc.noise_seed = 106;

    case 'G'
        sc.name = 'G: Diving attack';
        sc.description = 'Steep dive with vertical weave and blackout during maneuver.';
        sc.x0 = 1200; sc.y0 = 1500; sc.vx0 = -50; sc.vy0 = -90;
        sc.man.x = [ 10  2.0  Inf  1.0];
        sc.man.y = [ 14  1.4  Inf  1.0];
        sc.blackout = [2.20 0.40];
        sc.noise_seed = 107;

    case 'H'
        sc.name = 'H: Near-overhead';
        sc.description = 'Inbound toward overhead; reference reaches the +85 deg limit.';
        sc.x0 = 900; sc.y0 = 1200; sc.vx0 = -85; sc.vy0 = 0;
        sc.man.y = [  8  1.0  Inf  1.0];
        sc.noise_seed = 108;

    case 'I'
        sc.name = 'I: Fast long-range';
        sc.description = 'About 140 m/s from long range, velocity-reversing jinks, late blackout.';
        sc.x0 = 1950; sc.y0 = 300; sc.vx0 = -110; sc.vy0 = 85;
        sc.man.x = [ 12  2.0  2.0  1.0];
        sc.man.y = [-15  2.0  2.0  1.0];
        sc.blackout = [4.00 0.50];
        sc.noise_seed = 109;

    otherwise
        error('scenario_library:id', 'Unknown scenario id "%s".', id);
end
end
