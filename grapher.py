import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the '3d' projection)


def animate(U, D, TAR_log, LOS_view, LOS_log, X, tar_static, step=10, show=True):
    """Animate a completed engagement: 3D world view + 2D seeker view.

    U:         list of Top_world positions, one per timestep
    D:         list of Bottom_world positions, one per timestep
    TAR_log:   list of target positions, one per timestep
    LOS_view:  list of world-frame LOS unit vectors (for the boresight quiver)
    LOS_log:   list of [alpha, beta] seeker-plane bearings. These are RAW --
               the FOV gate reduces to alpha^2 + beta^2 <= F_len^2, so with
               F_len=1 they already live in the unit disc drawn as the FOV
               circle. Do not pass atan() of them.
    X:         list of scalar ranges to target, one per timestep
    tar_static: the target's initial position (plotted as a fixed marker)
    step:      plot every Nth timestep
    show:      call plt.show() before returning

    Returns (ani3d, ani2d). KEEP A REFERENCE to the returned animations if
    you pass show=False -- FuncAnimation objects that get garbage collected
    stop rendering, which is what the "Animation was deleted without
    rendering anything" warning means.
    """
    if not U or not TAR_log or not LOS_log:
        print("grapher.animate: no trajectory data to plot (empty buffers)")
        return None, None

    # unpack trajectory
    x, y, z = map(list, zip(*U))
    x2, y2, z2 = map(list, zip(*D))
    tx, ty, tz = map(list, zip(*TAR_log))
    lx, ly, lz = map(list, zip(*LOS_view))

    # LOS data
    los_x = [p[0] for p in LOS_log]
    los_y = [p[1] for p in LOS_log]

    n_frames = max(1, len(x) // step)

    # =========================
    # 3D WORLD VIEW
    # =========================
    fig3d = plt.figure()
    ax3d = fig3d.add_subplot(111, projection='3d')

    line1, = ax3d.plot([], [], [], 'red')
    line2, = ax3d.plot([], [], [], 'green')

    # moving target dot + trajectory line
    target_dot = ax3d.scatter([], [], [], c='red', s=20)
    target_line, = ax3d.plot([], [], [], 'orange', linestyle='--')

    target_static = ax3d.scatter(*tar_static, c='yellow')

    los_line, = ax3d.plot([], [], [], 'blue', alpha=0.5)

    # FIXED limits (include target too)
    ax3d.set_xlim(min(min(x), min(tx)), max(max(x), max(tx)))
    ax3d.set_ylim(min(min(y), min(ty)), max(max(y), max(ty)))
    ax3d.set_zlim(min(min(z), min(tz)), max(max(z), max(tz)))

    ax3d.set_xlabel("X")
    ax3d.set_ylabel("Y")
    ax3d.set_zlabel("Z")

    quiver = None

    def update3d(frame):
        nonlocal quiver
        i = min(frame * step, len(x) - 1)
        # missile trajectories
        line1.set_data(x[:i], y[:i])
        line1.set_3d_properties(z[:i])

        line2.set_data(x2[:i], y2[:i])
        line2.set_3d_properties(z2[:i])

        if quiver:
            quiver.remove()

        quiver = ax3d.quiver(x[i],y[i],z[i],lx[i],ly[i],lz[i],length=800,arrow_length_ratio=1,color="red")

        # moving target (dot)
        target_dot._offsets3d = ([tx[i]], [ty[i]], [tz[i]])

        # target trajectory
        target_line.set_data(tx[:i], ty[:i])
        target_line.set_3d_properties(tz[:i])

        # LOS line (missile → target)
        los_line.set_data([x[i], tx[i]], [y[i], ty[i]])
        los_line.set_3d_properties([z[i], tz[i]])

        return line1, line2, target_dot, target_line, los_line, quiver

    ani3d = animation.FuncAnimation(
        fig3d,
        update3d,
        frames=n_frames,
        interval=100,
        blit=False
    )

    # =========================
    # SEEKER VIEW (2D)
    # =========================
    fig2d, ax2d = plt.subplots()

    value_text = ax2d.text(
        0.02, 0.95, '',
        transform=ax2d.transAxes,
        fontsize=10,
        verticalalignment='top'
    )

    point, = ax2d.plot([], [], 'ro', markersize=5)
    trail, = ax2d.plot([], [], 'r-', alpha=0.4)

    circle = plt.Circle((0, 0), 1, fill=False)   # FOV boundary
    ax2d.add_patch(circle)

    ax2d.set_xlim(-1, 1)
    ax2d.set_ylim(-1, 1)
    ax2d.set_aspect('equal')

    def update2d(frame):
        i = min(frame * step, len(los_x) - 1)

        x_val = los_x[i]
        y_val = los_y[i]

        point.set_data([x_val], [y_val])
        trail.set_data(los_x[:i], los_y[:i])

        value_text.set_text(f"Range: {X[i]:.2f} m")

        return point, trail, value_text

    ani2d = animation.FuncAnimation(
        fig2d,
        update2d,
        frames=max(1, len(los_x) // step),
        interval=100,
        blit=True
    )

    if show:
        plt.show()

    return ani3d, ani2d
