"""Minimal base version of rk4_mpmath.py for one initial condition from sol_80.txt.

Integrate from 0 to period_count * T using adaptive RK4 and print the final
results. Plot samples are kept in memory; no results or checkpoints are written
to disk. The orbit plot is displayed after integration.
"""

from pathlib import Path
from time import perf_counter

import mpmath as mp


# =============================
# Hyperparameters – edit here
# =============================
SOL80_TXT = Path(__file__).resolve().parent / "Listor" / "sol_80.txt"

HYPERPARAMETERS = {
    "ic_number": 100,  # <<< SELECT IC HERE: row number in sol_80.txt, starting at 1. >>>
    "period_count": 1,  # <<< SELECT THE NUMBER OF PERIODS HERE: end time = period_count * T. >>>
    "dps": 50,
    "h0": "1e-3",  # Initial step size; subsequent steps are chosen adaptively.
    "rtol": "1e-23", # Relative tolerance  
    "atol": "1e-23", # Absolute tolerance 
    "periodicity_threshold": "0.1",
    "safety": "0.9",
    "min_factor": "0.2",
    "max_factor": "2.0",
    "max_steps": 35_000_000,
    "progress_every": 10_000,  # <<< PRINT PROGRESS EVERY N ACCEPTED STEPS. >>>
    "save_every": 3_000,  # <<< KEEP A PLOT SAMPLE IN MEMORY EVERY N ACCEPTED STEPS. >>>
    "view_tilt_degrees": 15,  # <<< VIEW TILT FROM DIRECTLY ABOVE: 0–90 DEGREES. >>>
    "plot_linewidth": 1.26,  # Static orbit line width in points (10% thinner).
}


def validate_hyperparameters(settings):
    for key in ("ic_number", "period_count", "dps", "max_steps",
                "progress_every", "save_every"):
        if type(settings[key]) is not int or settings[key] < 1:
            raise ValueError(f"{key} must be a positive integer.")
    with mp.workdps(settings["dps"]):
        for key in ("h0", "rtol", "atol", "periodicity_threshold",
                    "safety", "min_factor", "max_factor"):
            value = mp.mpf(settings[key])
            if not mp.isfinite(value) or value <= 0:
                raise ValueError(f"{key} must be finite and positive.")
        if not mp.mpf(settings["safety"]) < 1:
            raise ValueError("safety must be less than 1.")
        if not mp.mpf(settings["min_factor"]) < 1 < mp.mpf(settings["max_factor"]):
            raise ValueError("Requires min_factor < 1 < max_factor.")
    if not 0 <= settings["view_tilt_degrees"] <= 90:
        raise ValueError("view_tilt_degrees must be between 0 and 90.")
    if not mp.isfinite(settings["plot_linewidth"]) or settings["plot_linewidth"] <= 0:
        raise ValueError("plot_linewidth must be finite and positive.")


def load_initial_condition(path, ic_number):
    """Read x, y and T from the selected source row without conversion to float."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not 1 <= ic_number <= len(lines):
        raise ValueError(f"The IC number must be between 1 and {len(lines)}.")
    columns = lines[ic_number - 1].split()
    if len(columns) != 4:
        raise ValueError(f"Row {ic_number} must contain x y T T*.")
    x, y, period = [mp.mpf(value) for value in columns[:3]]
    if not all(mp.isfinite(value) for value in (x, y, period)) or period <= 0:
        raise ValueError("x, y and T must be finite, and T must be positive.")
    return x, y, period


# =============================
# Initial state and equations of motion
# =============================
def initial_conditions_mp(x, y):
    # G = 1, three unit masses and zero initial velocities (free fall).
    # State: [x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3].
    zero = mp.mpf("0")
    return [x, y, mp.mpf("-0.5"), zero, mp.mpf("0.5"), zero] + [zero] * 6


def vec_sub(a, b):
    return [a[i] - b[i] for i in range(2)]


def vec_add(a, b):
    return [a[i] + b[i] for i in range(2)]


def vec_norm(a):
    return mp.sqrt(sum(ai**2 for ai in a))


def derivatives_mp(t, state):
    r1, r2, r3 = state[0:2], state[2:4], state[4:6]

    def grav_acc(r_from, r_to):
        diff = vec_sub(r_to, r_from)
        dist = vec_norm(diff)
        return [diff[i] / dist**3 for i in range(2)]

    a1 = vec_add(grav_acc(r1, r2), grav_acc(r1, r3))
    a2 = vec_add(grav_acc(r2, r1), grav_acc(r2, r3))
    a3 = vec_add(grav_acc(r3, r1), grav_acc(r3, r2))
    return state[6:12] + a1 + a2 + a3


# =============================
# RK4 with adaptive step size
# =============================
def rk4_step_mp(f, t, y, h):
    f1 = f(t, y)
    f2 = f(t + h/2, [y[i] + h * f1[i] / 2 for i in range(len(y))])
    f3 = f(t + h/2, [y[i] + h * f2[i] / 2 for i in range(len(y))])
    f4 = f(t + h, [y[i] + h * f3[i] for i in range(len(y))])
    return [
        y[i] + h * (f1[i] + 2*f2[i] + 2*f3[i] + f4[i]) / 6
        for i in range(len(y))
    ]


def solve_adaptive_rk4_mp(f, y0, T, *, rtol, atol, h0, safety,
                          min_factor, max_factor, max_steps,
                          progress_every, save_every):
    t = mp.mpf("0")
    h = mp.mpf(h0)
    y = list(y0)
    steps = 0
    rejected = 0
    stop_reason = "completed"
    # Plot sampling does not affect adaptive integration.
    t_values = [t]
    y_values = [list(y)]

    while t < T:
        if steps >= max_steps:
            stop_reason = "max_steps"
            print(f"Stopped at t = {mp.nstr(t, 12)}: step limit {max_steps} reached.", flush=True)
            break

        # Only shorten the final step to reach T; no fixed integration subintervals.
        reaches_end = t + h >= T
        if reaches_end:
            h = T - t
        if t + h == t or t + h/2 == t:
            raise ArithmeticError("The step size is too small for the selected dps.")

        y_big = rk4_step_mp(f, t, y, h)
        y_half = rk4_step_mp(f, t, y, h / 2)
        y_small = rk4_step_mp(f, t + h / 2, y_half, h / 2)

        # Use the original local error control solely to choose the step size.
        err = max(
            abs(y_small[i] - y_big[i])
            / (atol + rtol * max(abs(y[i]), abs(y_small[i])))
            for i in range(len(y))
        )
        if not mp.isfinite(err):
            raise ArithmeticError("The RK4 step produced a non-finite error estimate.")

        if err <= 1:
            t = T if reaches_end else t + h
            y = y_small
            steps += 1

            if steps % save_every == 0:
                t_values.append(t)
                y_values.append(list(y))

            if steps % progress_every == 0:
                print(
                    f"step {steps}, "
                    f"t = {mp.nstr(t, 12)} / {mp.nstr(T, 12)}, "
                    f"{mp.nstr(100 * t / T, 6)} %, "
                    f"h = {mp.nstr(h, 5)}, "
                    f"err = {mp.nstr(err, 5)}, "
                    f"rejected = {rejected}",
                    flush=True,
                )
        else:
            rejected += 1

        if err == 0:
            factor = max_factor
        else:
            factor = safety * err**(-mp.mpf("1") / 5)
            factor = min(max_factor, max(min_factor, factor))
        h *= factor

    # Always include the initial and final points, even at the step limit.
    if t_values[-1] != t:
        t_values.append(t)
        y_values.append(list(y))
    return {
        "t": t,
        "y": y,
        "steps": steps,
        "rejected": rejected,
        "stop_reason": stop_reason,
        "t_values": t_values,
        "y_values": y_values,
    }


# =============================
# Proximity functions – evaluated only at the final point
# =============================
def pairwise_distances_mp(state):
    r1, r2, r3 = state[0:2], state[2:4], state[4:6]
    return sorted([vec_norm(vec_sub(r1, r2)),
                   vec_norm(vec_sub(r1, r3)),
                   vec_norm(vec_sub(r2, r3))])


def proximity_values_mp(state, initial_state):
    return {
        "shape_error": mp.sqrt(sum(
            (a - b)**2 for a, b in zip(pairwise_distances_mp(state),
                                      pairwise_distances_mp(initial_state)))),
        "position_error": vec_norm([state[i] - initial_state[i] for i in range(6)]),
        "velocity_error": vec_norm([state[i] - initial_state[i] for i in range(6, 12)]),
        "state_error": vec_norm([state[i] - initial_state[i] for i in range(12)]),
    }


def configure_orbit_view(ax, positions, settings):
    """Use equal x/y length scales and the selected camera tilt for planar motion."""
    import numpy as np

    low = positions.min(axis=(0, 1))
    high = positions.max(axis=(0, 1))
    center = (low + high) / 2
    extent = max(float(np.max(high - low)), 1e-12) * 0.6
    ax.set_xlim(center[0] - extent, center[0] + extent)
    ax.set_ylim(center[1] - extent, center[1] + extent)
    ax.set_zlim(-0.2 * extent, 0.2 * extent)
    ax.set_box_aspect((1, 1, 0.2))
    ax.view_init(elev=90 - settings["view_tilt_degrees"], azim=-90)
    ax.set(xlabel="x", ylabel="y")
    ax.set_zticks([])
    ax.zaxis.line.set_visible(False)
    ax.grid(True, alpha=0.2)


def add_periodicity_box(fig, result, settings):
    """Show the final return test, threshold and state error (two significant figures)."""
    values = result.get("proximity_values")
    if values is None:
        values = proximity_values_mp(result["y_values"][-1], result["y_values"][0])
    threshold = result.get("periodicity_threshold", mp.mpf(settings["periodicity_threshold"]))
    if result["stop_reason"] != "completed":
        classification = "not assessed (incomplete)"
    else:
        classification = "periodic" if all(value < threshold for value in values.values()) else "non-periodic"
    text = (
        f"Periodicity (all four): {classification}\n"
        f"Threshold: {mp.nstr(threshold, mp.mp.dps)}\n"
        f"Final state_error: {mp.nstr(values['state_error'], 2, strip_zeros=False)}"
    )
    return fig.text(
        0.04, 0.90, text, ha="left", va="top", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="0.6", alpha=0.95),
    )


def create_orbit_plot(result, settings):
    """Create an in-memory orbit plot with body colours and start/end markers."""
    import numpy as np
    import matplotlib.pyplot as plt

    positions = np.array([
        [float(v) for v in state[:6]] for state in result["y_values"]
    ]).reshape(-1, 3, 2)
    if not np.all(np.isfinite(positions)):
        raise ValueError("Trajectory data cannot be displayed as finite float64 values.")
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    fig.subplots_adjust(bottom=0.1, top=0.79)
    configure_orbit_view(ax, positions, settings)
    add_periodicity_box(fig, result, settings)
    for body, color in enumerate(("#0072B2", "#D55E00", "#009E73")):
        path = positions[:, body, :]
        ax.plot(path[:, 0], path[:, 1], np.zeros(len(path)), color=color,
                linewidth=settings["plot_linewidth"], label=f"Body {body + 1}")
        ax.plot([path[0, 0]], [path[0, 1]], [0], "o", color=color, markersize=8)
        ax.plot([path[-1, 0]], [path[-1, 1]], [0], "x", color=color, markersize=8)
    status = "" if result["stop_reason"] == "completed" else " · incomplete integration"
    fig.suptitle(f"Orbits · IC {settings['ic_number']}{status}", fontsize=15)
    ax.set_title(f"0 ≤ t ≤ {mp.nstr(result['t'], 8)} · tilt {settings['view_tilt_degrees']:g}°")
    ax.legend(loc="upper right")
    fig.text(0.5, 0.05, "○ Start  ·  × End  ·  motion in the plane z = 0", ha="center")
    return fig


def show_orbit_plot(result, settings):
    import matplotlib.pyplot as plt

    print("Displaying the orbit plot.")
    fig = create_orbit_plot(result, settings)
    plt.show(block=True)
    return fig


def main():
    validate_hyperparameters(HYPERPARAMETERS)
    mp.mp.dps = HYPERPARAMETERS["dps"]
    ic_number = HYPERPARAMETERS["ic_number"]
    x, y, T = load_initial_condition(SOL80_TXT, ic_number)
    period_count = HYPERPARAMETERS["period_count"]
    end_time = period_count * T
    initial_state = initial_conditions_mp(x, y)
    threshold = mp.mpf(HYPERPARAMETERS["periodicity_threshold"])

    print(f"IC {ic_number} from {SOL80_TXT.name}")
    print(f"x = {mp.nstr(x, mp.mp.dps)}, y = {mp.nstr(y, mp.mp.dps)}")
    print(f"Number of periods: {period_count}, T = {mp.nstr(T, mp.mp.dps)}")
    print(f"Integrating from 0 to {mp.nstr(end_time, mp.mp.dps)}", flush=True)

    start = perf_counter()
    result = solve_adaptive_rk4_mp(
        derivatives_mp, initial_state, end_time,
        **{key: mp.mpf(HYPERPARAMETERS[key]) for key in
           ("rtol", "atol", "h0", "safety", "min_factor", "max_factor")},
        max_steps=HYPERPARAMETERS["max_steps"],
        progress_every=HYPERPARAMETERS["progress_every"],
        save_every=HYPERPARAMETERS["save_every"],
    )
    print(f"Accepted steps: {result['steps']}")
    print(f"Rejected steps: {result['rejected']}")
    print(f"Elapsed time: {perf_counter() - start:.2f} s")
    t, state, stop_reason = result["t"], result["y"], result["stop_reason"]

    values = proximity_values_mp(state, initial_state)
    result["proximity_values"] = values
    result["periodicity_threshold"] = threshold
    print(f"Final point: t = {mp.nstr(t, mp.mp.dps)}")
    for name, value in values.items():
        print(f"  {name}: {mp.nstr(value, mp.mp.dps)}")
    print(f"Periodicity threshold: {mp.nstr(threshold, mp.mp.dps)}")

    if stop_reason != "completed":
        print("Incomplete integration: the step limit was reached before the end time. Periodicity was not classified.")
    else:
        # Use the same two numerical classifications as in rk4_mpmath.py.
        class_1 = "periodic" if values["state_error"] < threshold else "non-periodic"
        class_2 = "periodic" if all(v < threshold for v in values.values()) else "non-periodic"
        print(f"Class 1 (state_error): {class_1}")
        print(f"Class 2 (all four): {class_2}")
    show_orbit_plot(result, HYPERPARAMETERS)
    return result


if __name__ == "__main__":
    main()
