"""Base version V1.5: adaptive RK4 with an orbit plot and interactive animation.

Integrate from 0 to period_count * T using adaptive RK4 and print the final
results. Plot samples are kept in memory; no results or checkpoints are written
to disk. The orbit plot and animation are displayed after integration.

The animation is based on the point and trail display in 3_Body_Problem.ipynb,
adapted to planar mpmath data and an interactive playback-speed slider.
Run: python3 rk4_mpmath_base_V1.5.py
"""

from pathlib import Path
from time import perf_counter

import mpmath as mp


# =============================
# Hyperparameters – edit here
# =============================
SOL80_TXT = Path(__file__).resolve().parent / "Listor" / "sol_80.txt"

HYPERPARAMETERS = {
    "ic_number": 581,  # <<< SELECT IC HERE: row number in sol_80.txt, starting at 1. >>>
    "period_count": 1,  # <<< SELECT THE NUMBER OF PERIODS HERE: end time = period_count * T. >>>
    "dps": 50,
    "h0": "1e-3",  # Initial step size; subsequent steps are chosen adaptively.
    "rtol": "1e-23",
    "atol": "1e-23",
    "periodicity_threshold": "0.1",
    "safety": "0.9",
    "min_factor": "0.2",
    "max_factor": "2.0",
    "max_steps": 35_000_000,
    "progress_every": 10_000,  # <<< PRINT PROGRESS EVERY N ACCEPTED STEPS. >>>
    "save_every": 3_000,  # <<< KEEP A PLOT SAMPLE IN MEMORY EVERY N ACCEPTED STEPS. >>>
    "view_tilt_degrees": 15,  # <<< VIEW TILT FROM DIRECTLY ABOVE: 0–90 DEGREES. >>>
    "plot_linewidth": 1.26,  # Static orbit line width in points (10% thinner).
    "animation_linewidth": 1.8,  # Animated orbit trail width in points (10% thinner).
    "triangle_linewidth": 2.25,  # Triangle side width in points (10% thinner).
    "animate": True,  # <<< SHOW THE ANIMATION AFTER INTEGRATION. >>>
    "animation_mode": "triangle",  # <<< DISPLAY MODE: "triangle" (coloured sides) or "orbits" (trajectories). >>>
    "animation_fps": 30,  # Animation frame rate.
    "animation_seconds": 30,  # Playback seconds for the complete trajectory at 1x speed.
    "animation_speed": 1.0,  # Initial speed; the slider ranges from 0 (paused) to 5x.
    "animation_trail_fraction": 0.15,  # Initial trail fraction: 0 = bodies only, 1 = full elapsed trail.
}

# Optional manual input. Uncomment the complete block below to override ic_number.
# Keep decimal values in quotes so they are converted directly to mpmath numbers.
# T is the reference period; total integration time is period_count * T.
MANUAL_IC = None
# MANUAL_IC = {
#     "x": "0.25694268581619406691269277542214527477259843827894882545958420075540861203009715",
#     "y": "0.2783746621563187640173069301772599105064839871107602666586336860267930150563534",
#     "T": "0.21573615905021502937878732480336592981156797723875640215934585591298459904348956e1",
# }


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
    for key in ("plot_linewidth", "animation_linewidth", "triangle_linewidth"):
        if not mp.isfinite(settings[key]) or settings[key] <= 0:
            raise ValueError(f"{key} must be finite and positive.")
    if type(settings["animate"]) is not bool:
        raise ValueError("animate must be True or False.")
    if settings["animation_mode"] not in {"triangle", "orbits"}:
        raise ValueError('animation_mode must be "triangle" or "orbits".')
    for key in ("animation_fps", "animation_seconds"):
        value = mp.mpf(settings[key])
        if not mp.isfinite(value) or value <= 0:
            raise ValueError(f"{key} must be finite and positive.")
    if not 0 <= settings["animation_speed"] <= 5:
        raise ValueError("animation_speed must be between 0 and 5.")
    if not 0 <= settings["animation_trail_fraction"] <= 1:
        raise ValueError("animation_trail_fraction must be between 0 and 1.")


def load_initial_condition(path, ic_number):
    """Read x, y and T from the selected source row without conversion to float."""
    return parse_initial_condition(path.read_text(encoding="utf-8"), ic_number)


def parse_initial_condition(source_text, ic_number):
    """Parse a row from local or downloaded source text at the current precision."""
    lines = source_text.splitlines()
    if not 1 <= ic_number <= len(lines):
        raise ValueError(f"The IC number must be between 1 and {len(lines)}.")
    columns = lines[ic_number - 1].split()
    if len(columns) != 4:
        raise ValueError(f"Row {ic_number} must contain x y T T*.")
    x, y, period = [mp.mpf(value) for value in columns[:3]]
    if not all(mp.isfinite(value) for value in (x, y, period)) or period <= 0:
        raise ValueError("x, y and T must be finite, and T must be positive.")
    return x, y, period


def get_initial_condition(settings):
    """Use the catalogue row unless the optional manual block is enabled."""
    if MANUAL_IC is None:
        values = load_initial_condition(SOL80_TXT, settings["ic_number"])
        return (*values, f"IC {settings['ic_number']}")
    if set(MANUAL_IC) != {"x", "y", "T"} or any(
            not isinstance(value, str) for value in MANUAL_IC.values()):
        raise ValueError('MANUAL_IC must contain x, y and T as decimal strings.')
    x, y, period = [mp.mpf(MANUAL_IC[key]) for key in ("x", "y", "T")]
    if not all(mp.isfinite(value) for value in (x, y, period)) or period <= 0:
        raise ValueError("Manual x, y and T must be finite, and T must be positive.")
    return x, y, period, "Manual initial condition"


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
    label = result.get("ic_label", f"IC {settings['ic_number']}")
    fig.suptitle(f"Orbits · {label}{status}", fontsize=15)
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


class OrbitAnimation:
    """Float64 display copies; mpmath results and classifications are never changed."""

    def __init__(self, result, settings):
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation
        from matplotlib.widgets import Button, RadioButtons, Slider

        self.np = np
        self.times = np.array([float(t) for t in result["t_values"]])
        self.positions = np.array([
            [float(v) for v in state[:6]] for state in result["y_values"]
        ]).reshape(-1, 3, 2)
        if not np.all(np.isfinite(self.times)) or not np.all(np.isfinite(self.positions)):
            raise ValueError("Trajectory data cannot be displayed as finite float64 values.")
        if len(self.times) < 2 or self.times[-1] <= self.times[0]:
            raise ValueError("The animation requires at least two distinct times.")
        # Very close mpmath times may round to the same float64 time.
        keep = np.r_[np.diff(self.times) > 0, True]
        self.times, self.positions = self.times[keep], self.positions[keep]
        self.duration = self.times[-1] - self.times[0]
        self.simulation_rate = self.duration / settings["animation_seconds"]
        self.speed = settings["animation_speed"]
        self.mode = settings["animation_mode"]
        self.trail_fraction = settings["animation_trail_fraction"]
        self.trail_duration = self.trail_fraction * self.duration
        self.current_time = self.times[0]
        self.last_tick = perf_counter()

        self.fig = plt.figure(figsize=(8, 9))
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.fig.subplots_adjust(bottom=0.35, top=0.79)
        label = result.get("ic_label", f"IC {settings['ic_number']}")
        self.title_prefix = f"Three bodies · {label}"
        configure_orbit_view(self.ax, self.positions, settings)
        self.periodicity_box = add_periodicity_box(self.fig, result, settings)
        if result["stop_reason"] != "completed":
            self.title_prefix += " · incomplete integration"
        self.fig.suptitle(f"{self.title_prefix} · tilt {settings['view_tilt_degrees']:g}°", fontsize=15)

        colors = ("#0072B2", "#D55E00", "#009E73")
        self.trails, self.points = [], []
        for body, color in enumerate(colors):
            trail, = self.ax.plot([], [], [], color=color, linewidth=settings["animation_linewidth"])
            point, = self.ax.plot([], [], [], "o", color=color, markersize=9,
                                  markeredgecolor="white", label=f"Body {body + 1}")
            self.trails.append(trail)
            self.points.append(point)
        self.triangle_sides = []
        self.pairs = ((0, 1), (1, 2), (2, 0))
        for (first, second), color in zip(self.pairs, colors):
            side, = self.ax.plot([], [], [], color=color, linewidth=settings["triangle_linewidth"],
                                 label=f"Side {first + 1}–{second + 1}")
            self.triangle_sides.append(side)
        self.legend = None
        self.time_label = self.ax.set_title("")
        self.slider = Slider(
            self.fig.add_axes([0.25, 0.22, 0.55, 0.025]), "Speed",
            0, 5, valinit=self.speed, valfmt="%1.2fx",
        )
        self.tilt_slider = Slider(
            self.fig.add_axes([0.25, 0.165, 0.55, 0.025]), "Tilt from z",
            0, 90, valinit=settings["view_tilt_degrees"], valfmt="%1.1f°",
        )
        self.trail_slider = Slider(
            self.fig.add_axes([0.25, 0.11, 0.55, 0.025]), "Trail length",
            0, 1, valinit=self.trail_fraction, valfmt="%1.2f",
        )
        self.mode_control = RadioButtons(
            self.fig.add_axes([0.15, 0.015, 0.32, 0.055]), ("Orbits", "Triangle"),
            active=1 if self.mode == "triangle" else 0,
        )
        self.restart_button = Button(self.fig.add_axes([0.59, 0.02, 0.22, 0.04]), "Restart")
        self.fig.text(0.5, 0.285, "0x = paused  ·  1x = full trajectory in "
                      f"{settings['animation_seconds']:g} s", ha="center", fontsize=10)
        self.fig.text(0.5, 0.083, "Orbits: 0 = bodies only · 1 = full elapsed trail", ha="center", fontsize=9)
        self.slider.on_changed(self.set_speed)
        self.tilt_slider.on_changed(self.set_tilt)
        self.trail_slider.on_changed(self.set_trail_fraction)
        self.restart_button.on_clicked(self.restart)
        self.mode_control.on_clicked(self.set_mode)
        self.set_mode("Triangle" if self.mode == "triangle" else "Orbits")
        self.draw()
        # Keep animation and controls as attributes to prevent garbage collection.
        self.animation = FuncAnimation(
            self.fig, self.update, init_func=self.draw, frames=None,
            interval=1000 / settings["animation_fps"], blit=False,
            cache_frame_data=False,
        )

    def position_at(self, t):
        """Linearly interpolate display positions between saved samples."""
        np = self.np
        index = int(np.clip(np.searchsorted(self.times, t, side="right") - 1,
                            0, len(self.times) - 2))
        fraction = np.clip((t - self.times[index])
                           / (self.times[index + 1] - self.times[index]), 0, 1)
        return ((1 - fraction) * self.positions[index]
                + fraction * self.positions[index + 1])

    def draw(self):
        np = self.np
        t = self.current_time
        start_time = max(self.times[0], t - self.trail_duration)
        left = np.searchsorted(self.times, start_time, side="right")
        right = np.searchsorted(self.times, t, side="left")
        current = self.position_at(t)
        trail = np.concatenate((self.position_at(start_time)[None, :, :],
                                self.positions[left:right], current[None, :, :]))
        for body, (line, point) in enumerate(zip(self.trails, self.points)):
            line.set_data_3d(trail[:, body, 0], trail[:, body, 1], np.zeros(len(trail)))
            line.set_visible(self.mode == "orbits" and self.trail_fraction > 0)
            point.set_data_3d([current[body, 0]], [current[body, 1]], [0])
        self.time_label.set_text(f"t = {t:.6g} / {self.times[-1]:.6g}")
        for side, (first, second) in zip(self.triangle_sides, self.pairs):
            side.set_data_3d([current[first, 0], current[second, 0]],
                             [current[first, 1], current[second, 1]], [0, 0])
        return (*self.trails, *self.points, *self.triangle_sides, self.time_label)

    def set_mode(self, label):
        self.mode = "triangle" if label == "Triangle" else "orbits"
        triangle = self.mode == "triangle"
        for line in self.trails:
            line.set_visible(not triangle and self.trail_fraction > 0)
        for side in self.triangle_sides:
            side.set_visible(triangle)
        if self.legend is not None:
            self.legend.remove()
        self.legend = self.ax.legend(
            handles=self.triangle_sides if triangle else self.points,
            loc="upper right",
        )
        self.fig.canvas.draw_idle()

    def advance(self, elapsed):
        if self.speed > 0:
            # Show the final point for one frame before restarting playback.
            if self.current_time >= self.times[-1]:
                self.current_time = self.times[0]
            else:
                self.current_time = min(self.times[-1], self.current_time
                                        + max(0, elapsed) * self.simulation_rate * self.speed)
        return self.draw()

    def update(self, frame):
        now = perf_counter()
        elapsed = now - self.last_tick
        self.last_tick = now
        return self.advance(elapsed)

    def set_speed(self, value):
        self.update(None)
        self.speed = float(value)
        self.fig.canvas.draw_idle()

    def set_tilt(self, value):
        """Change the viewing angle from the z-axis without changing the motion."""
        self.ax.view_init(elev=90 - float(value), azim=self.ax.azim)
        self.fig._suptitle.set_text(f"{self.title_prefix} · tilt {float(value):g}°")
        self.fig.canvas.draw_idle()

    def set_trail_fraction(self, value):
        """Update the visible time window; never draw future trajectory samples."""
        self.trail_fraction = float(value)
        self.trail_duration = self.trail_fraction * self.duration
        self.draw()
        self.fig.canvas.draw_idle()

    def restart(self, event):
        self.current_time = self.times[0]
        self.last_tick = perf_counter()
        self.draw()
        self.fig.canvas.draw_idle()


def show_animation(result, settings):
    import matplotlib.pyplot as plt

    if len(result["t_values"]) < 2:
        print("No animation: too few saved time samples.")
        show_orbit_plot(result, settings)
        return
    print("Displaying the orbit plot and animation. The speed slider pauses at 0x.")
    print("Select Orbits or Triangle. The triangle sides have three different colours.")
    print("Use Tilt from z to change the view, and Trail length to adjust orbit persistence (0–1).")
    print("Plot detail is controlled by save_every; smaller values retain more samples.")
    orbit_fig = create_orbit_plot(result, settings)
    view = OrbitAnimation(result, settings)
    plt.show(block=True)
    return view


def main(*, show_plots=True):
    validate_hyperparameters(HYPERPARAMETERS)
    mp.mp.dps = HYPERPARAMETERS["dps"]
    x, y, T, label = get_initial_condition(HYPERPARAMETERS)
    period_count = HYPERPARAMETERS["period_count"]
    end_time = period_count * T
    initial_state = initial_conditions_mp(x, y)
    threshold = mp.mpf(HYPERPARAMETERS["periodicity_threshold"])

    print(label)
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
    result["ic_label"] = label
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
    if show_plots:
        if HYPERPARAMETERS["animate"]:
            result["animation_view"] = show_animation(result, HYPERPARAMETERS)
        else:
            show_orbit_plot(result, HYPERPARAMETERS)
    return result


if __name__ == "__main__":
    main()
