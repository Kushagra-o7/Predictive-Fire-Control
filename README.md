# Guardian Point-Defense: Predictive Turret Control (ME319 Assignment 2)

This is the MATLAB implementation of the predictive-feedback fire-control system for the
turret elevation axis: radar → target-state estimation → prediction (Tp = 0.5 s) →
predicted elevation reference → LTI compensator → turret.

## Quick start

```matlab
>> main_guardian          % runs everything, makes figures, writes results/
```

For quicker or batch runs, set options before calling the master script:

```matlab
>> OPT_OVERRIDE = struct('run_monte_carlo', false, 'run_comparison', false);
>> main_guardian
```

Tested on MATLAB R2025b. No toolboxes are needed, because the analysis is written
without the Control System Toolbox.

## Files

| File | Role |
|---|---|
| `main_guardian.m` | **Master script.** Owns the simulation clock (1 kHz). Runs every test case, models the laser, integrates the plant with RK4, then computes metrics, plots and summaries. |
| `target_motion_generator.m` | *Prescribed interface.* 20 Hz noisy radar measurement and the measurement-loss interval. |
| `turret_controller.m` | *Prescribed interface.* EKF estimator, predictor, elevation reference, compensator Gc(s) and fire logic. |
| `turret_plant.m` | *Prescribed interface.* First-principles DC motor, 100:1 gearbox and turret dynamics with hard stops. |
| `guardian_params.m` | Every physical and design parameter, in one place. |
| `design_compensator.m` | Analytical design of Gc(s) (delay matching), plus the comparison design. |
| `analyze_controller_design.m` | Margins, closed-loop poles, Bode plot, root locus, step and ramp responses, delay-matching check. |
| `scenario_library.m` | Test cases A–I. |
| `target_truth.m` | Exact target kinematics. Used for measurements, hit tests and error analysis. |
| `check_trajectory_limits.m` | Checks the range, speed, acceleration and measurement-loss bounds. |
| `random_scenario.m` | Random admissible engagements, for the Monte-Carlo generalization study. |
| `compute_metrics.m` | Interception, accuracy, estimation, prediction, measurement-loss and effort metrics. |
| `plot_scenario_results.m`, `plot_summary.m`, `plot_monte_carlo.m`, `save_fig.m` | Plotting. |
| `make_paper_figures.m` | Compact IEEE-column figures (`paper/figures`) and LaTeX table rows (`paper/tables`). Run it after `main_guardian` with `run_sensitivity = true`. |
| `paper/main.tex` | The IEEE conference paper (IEEEtran). Compile with pdflatex or on Overleaf. |

## Reproducing the paper

```matlab
>> OPT_OVERRIDE = struct('run_sensitivity', true); main_guardian   % about 7 min
>> make_paper_figures
```

Then compile `paper/main.tex`. On Overleaf: upload `paper.zip` and set `main.tex` as the main document.

## Design summary

### 1. Plant (derived in `turret_plant.m` / `guardian_params.m`)
Everything is referred to the turret side. The reflected inertia is J_eq = J_turret + N²J_m = 3200 kg m².

```
Gp(s) = N Kt / ( s [ (Lm s + Rm)(J_eq s + b) + N² Kt Ke ] )
      = 1000 / ( s (320 s² + 3350 s + 1 001 500) )
```

There is a free integrator, plus a lightly damped electromechanical resonance at 55.9 rad/s
(ζ = 0.094). The model is linear, so the small-perturbation model is exact away from the stops.

### 2. Estimator
An extended Kalman filter on the Cartesian state `[x vx ax y vy ay]`, using a Singer model
for each axis (acceleration correlated with time constant τm, RMS σm), with a polar `(r, θ)`
measurement. The Singer model fits a target with bounded, smoothly varying maneuvers. During
radar loss the filter runs open loop and its covariance grows, which automatically blocks
firing on a stale track.

### 3. Prediction and reference
The estimate is propagated to `t + Tp` from the time of the last update, giving
θ_ref = atan2(ŷ(t+Tp), x̂(t+Tp)), clipped to [−10°, 85°].

### 4. Compensator: the delay-matched loop
The laser arrives instantly, so a hit needs **θ_turret(t) = θ_target(t)**. The pointing error is

```
E(s) = [1 − T(s)·e^{sTp}] θ_target(s)
```

Gc(s) = Kc/(τc s + 1) is chosen so that the closed loop T(s) matches e^{−sTp} to second order:

```
Kc = d0/(N Kt Tp) = 2003 V/rad            (Kv = 1/Tp: ramp lag equals the lead)
τc = Tp/2 − d1/d0 = 0.2467 s              (curvature term matched as well)
```

So the **loop lag cancels the prediction lead**, and the beam points at the target's
*current* position. The pointing error is zero for a constant angular rate and a constant
angular acceleration. Results: PM ≈ 65°, GM ≈ 38 dB (the resonance is well suppressed),
dominant poles near a Butterworth pair (ωn ≈ 2.85 rad/s, ζ ≈ 0.71), peak voltage about 3–6×
lower than a high-gain type-2 loop.

**Note:** the tracking error θ_ref − θ_turret is deliberately about Tp·θ̇. The figures plot
both the tracking error and the true pointing error.

### 5. Fire logic (in `turret_controller.m`)
A pulse (a 0→1 edge on `fire_cmd`, held for a single call) is issued **only if both** hold:

1. **Recharged:** at least 1.0 s since the last pulse, and
2. **Inside the engagement radius:** the beam passes within R_eng = 2 m of the
   *estimated current* target position, with a confidence margin:
   `d_hat + k_σ·σ_perp ≤ R_eng`. Here σ_perp is the estimator's standard deviation across the
   beam. The track must also be at least 0.25 s old, and σ_perp must be no more than 0.75 m.

The master script models the laser separately. It ignores 0→1 edges that come before the
recharge time, and it scores each pulse against the **true** target position, using the
distance from the ray to the point.

## Test cases

| Id | Description |
|---|---|
| A | Nominal: constant velocity, level flight |
| B | Maneuvering: smooth weave in x and y from t = 1 s |
| C | Measurement loss: 1.00–1.50 s |
| D | B maneuver plus loss at the maneuver onset |
| E | High-g jinks close to a_max |
| F | Climbing crosser with high elevation rate, loss at 3 s |
| G | Diving attack with weave, loss at 2.2 s |
| H | Inbound toward overhead; the reference reaches the +85° stop |
| I | Fast (≈140 m/s) long-range target with velocity-reversing jinks, loss at 4 s |
| R### | Monte-Carlo random admissible cases (`random_scenario.m`) |

Every case is checked against the envelope: 800 ≤ r ≤ 2000 m, |v| ≤ 150 m/s, |a| ≤ 20 m/s²,
loss onset ≥ 1 s, loss duration ≤ 0.5 s.

## Master-script sections

- **A** Analytical controller design
- **B** Envelope check
- **C/D** Deterministic cases
- **E** Monte Carlo
- **F** Architecture comparison: final design vs a tight PI loop with and without prediction,
  and the delay-matched loop without prediction
- **H** Estimator sensitivity (off by default because it is slow)
- **G** Save the results: `results/summary_metrics.csv`, `results/guardian_results.mat` and
  `results/figures/*.png|.fig`

## Assumptions (state these in the paper)
- 2-D elevation plane, with the turret and radar at the origin; x is ground range, y is altitude.
- Initial turret state θ = 0 (stowed), θ̇ = 0, i = 0.
- Ideal rigid gearbox, damping on the turret side only, balanced trunnion (no gravity torque).
- No motor-voltage limit is given, so `P.V_max = Inf`. Peak voltage and current are reported;
  set a finite `V_max` to study saturation. The compensator has no integrator, so there is no windup.
- Radar noise: σ_r = 1 m, σ_θ = 0.02°. When `meas_valid = 0` the generator returns NaN.
- Controller time step: whatever rate the caller uses. The compensator is discretized exactly
  (ZOH) over the actual call interval, so the controller also works at 20 Hz.
