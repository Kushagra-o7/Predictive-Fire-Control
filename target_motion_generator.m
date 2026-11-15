function [r_meas, theta_meas, meas_valid] = target_motion_generator(t, scenario)
%TARGET_MOTION_GENERATOR  Simulated 20 Hz radar measurement of the V-44 target.
%
%   [r_meas, theta_meas, meas_valid] = target_motion_generator(t, scenario)
%
%   t         time [s] (may be called at any rate, e.g. every control step)
%   scenario  scenario id ('A', 'B', 'C', ... see scenario_library), an index,
%             or a full scenario struct (e.g. from random_scenario.m)
%
%   r_meas      measured range [m]            (NaN when meas_valid = 0)
%   theta_meas  measured elevation [rad]      (NaN when meas_valid = 0)
%   meas_valid  1 only on the first call at/after a new radar sample instant
%               t_k = k*Ts_radar that lies outside the measurement-loss
%               interval; 0 otherwise (no new measurement available).
%
%   Noise: zero-mean Gaussian, sigma_r = 1 m, sigma_theta = 0.02 deg. The noise
%   sequence is indexed by the sample number k and seeded per scenario, so a run
%   is exactly reproducible and independent of the calling rate.
%
%   Special call:  target_motion_generator('reset')            restart
%                  target_motion_generator('reset', seed)      restart and
%                  override the noise seed (Monte-Carlo studies).

persistent S

if ischar(t) || isstring(t)
    if strcmpi(char(t), 'reset')
        S = [];
        if nargin >= 2 && ~isempty(scenario)
            S.seed_override = scenario;
        end
    end
    r_meas = []; theta_meas = []; meas_valid = [];
    return
end

if nargin < 2 || isempty(scenario), scenario = 'A'; end
if isstruct(scenario), sc = scenario; else, sc = scenario_library(scenario); end

% (Re)initialize on first use, new scenario, or time running backwards
key = sprintf('%s_%d', sc.id, sc.noise_seed);
need_init = isempty(S) || ~isfield(S, 'key') || ~strcmp(S.key, key) || t < S.t_last - 1e-9;
if need_init
    seed_override = [];
    if ~isempty(S) && isfield(S, 'seed_override'), seed_override = S.seed_override; end
    S = struct();
    S.seed_override = seed_override;
    S.key    = key;
    S.P      = guardian_params();
    S.k_last = -1;
    S.t_last = -Inf;
    seed = sc.noise_seed;
    if ~isempty(seed_override), seed = seed_override; end
    S.noise = gaussian_noise(seed, 4000, 2);        % covers 200 s at 20 Hz
end
S.t_last = t;

r_meas = NaN; theta_meas = NaN; meas_valid = 0;

k = floor(t / S.P.Ts_radar + 1e-9);                   % latest sample index
if k <= S.k_last || k < 0
    return                                            % no new sample yet
end
S.k_last = k;
t_k = k * S.P.Ts_radar;

% Measurement loss (low-observability) interval
if ~isempty(sc.blackout)
    if t_k >= sc.blackout(1) - 1e-9 && t_k < sc.blackout(1) + sc.blackout(2) - 1e-9
        return
    end
end

tr = target_truth(t_k, sc);
idx = mod(k, size(S.noise, 1)) + 1;
r_meas     = tr.r     + S.P.sigma_r     * S.noise(idx, 1);
theta_meas = tr.theta + S.P.sigma_theta * S.noise(idx, 2);
meas_valid = 1;
end

function n = gaussian_noise(seed, rows, cols)
% Reproducible N(0,1) samples without disturbing the global RNG.
try
    s = RandStream('mt19937ar', 'Seed', seed);
    n = randn(s, rows, cols);
catch                                   % Octave fallback
    old = randn('state'); %#ok<RAND>
    randn('state', seed); %#ok<RAND>
    n = randn(rows, cols);
    randn('state', old); %#ok<RAND>
end
end
