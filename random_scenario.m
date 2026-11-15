function sc = random_scenario(seed, P)
%RANDOM_SCENARIO  Random admissible engagement for Monte-Carlo generalization.
%
%   sc = random_scenario(seed, P) draws (by rejection sampling) a target with a
%   random initial position, inbound velocity, 1-2 smooth maneuver components
%   per axis (random amplitude, frequency, duration and onset) and, with
%   probability 0.7, a random measurement-loss interval (onset in [1, 7] s,
%   duration in [0.1, 0.5] s). Every returned case satisfies
%   check_trajectory_limits. Used to test generalization beyond the
%   deterministic scenarios.

if nargin < 2, P = guardian_params(); end
try
    rs = RandStream('mt19937ar', 'Seed', seed);
    U = @(varargin) rand(rs, varargin{:});
catch
    rand('state', seed); %#ok<RAND>
    U = @(varargin) rand(varargin{:});
end
uni = @(lo, hi) lo + (hi - lo) * U();

for attempt = 1:500
    sc = struct('id', sprintf('R%03d', seed), 'name', sprintf('Random %d', seed), ...
        'description', 'Random admissible case', 'x0', 0, 'y0', 0, 'vx0', 0, 'vy0', 0, ...
        'man', struct('x', zeros(0,4), 'y', zeros(0,4)), 'blackout', [], ...
        'noise_seed', 10000 + seed);

    r0  = uni(1100, 1950);
    th0 = deg2rad(uni(5, 65));
    sc.x0 = r0 * cos(th0);  sc.y0 = r0 * sin(th0);
    spd = uni(50, 140);
    hdg = pi + th0 + deg2rad(uni(-70, 70));         % roughly inbound
    sc.vx0 = spd * cos(hdg); sc.vy0 = spd * sin(hdg);

    a_budget = uni(6, 19);                          % total acceleration budget
    share = U();
    for ax = 'xy'
        if ax == 'x', Aax = a_budget * share; else, Aax = a_budget * (1 - share); end
        nc = 1 + (U() > 0.6);
        comps = zeros(nc, 4);
        for k = 1:nc
            nper_set = [0.5 1 1.5 2 Inf];
            comps(k,:) = [sign(U() - 0.5) * Aax / nc, uni(0.6, 3.0), ...
                          nper_set(1 + floor(U() * numel(nper_set))), uni(1.0, 3.0)];
        end
        sc.man.(ax) = comps;
    end

    if U() < 0.7
        sc.blackout = [uni(1.0, 7.0), uni(0.1, 0.5)];
    end

    chk = check_trajectory_limits(sc, P);
    % additionally keep the target inside the turret's reachable elevation
    % range (with a 5 deg margin) so every case is physically engageable
    reachable = chk.theta_min_deg >= rad2deg(P.theta_min) + 5 && ...
                chk.theta_max_deg <= rad2deg(P.theta_max) - 5;
    if chk.pass && reachable
        return
    end
end
error('random_scenario:fail', 'Could not draw an admissible scenario for seed %d.', seed);
end
