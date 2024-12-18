import math
import time
import matplotlib
import numpy as np
import matplotlib.pyplot as plt  # from matplotlib import pyplot as plt


class CurveSimAnimation:

    def __init__(self, p, bodies, lightcurve):
        sampled_lightcurve = np.take(lightcurve, range(0, p.iterations, p.sampling_rate))  # Use only some of the calculated positions for the animation because it is so slow.
        self.fig, ax_right, ax_left, ax_lightcurve, self.red_dot = CurveSimAnimation.init_plot(p, sampled_lightcurve, lightcurve)  # Adjust constants in section [Plot] of config file to fit your screen.
        for body in bodies:  # Circles represent the bodies in the animation. Set their colors and add them to the matplotlib axis.
            body.circle_right.set_color(body.color)
            body.circle_left.set_color(body.color)
            ax_right.add_patch(body.circle_right)
            ax_left.add_patch(body.circle_left)
        self.render(p, bodies, lightcurve)

    @staticmethod
    def tic_delta(scope):
        """Returns a distance between two tics on an axis so that the total
        number of tics on that axis is between 5 and 10."""
        if scope <= 0:  # no or constant values
            return 1
        delta = 10 ** np.floor(math.log10(scope))
        if scope / delta < 5:
            if scope / delta < 2:
                return delta / 5
            else:
                return delta / 2
        else:
            return delta

    @staticmethod
    def init_plot(p, sampled_lightcurve, lightcurve):
        """Initialize the matplotlib figure containing 3 axis:
        Top left: projection (x,y,z) -> (x,z), order = -y.
        Top right: projection (x,y,z) -> (x,y), order = z.
        Bottom: lightcurve"""
        fig = plt.figure()
        fig.set_figwidth(p.figure_width)
        fig.set_figheight(p.figure_height)
        fig.set_facecolor("black")  # background color outside of ax_left and ax_lightcurve
        buffer = 0
        fig.subplots_adjust(left=buffer, right=1.0 - buffer, bottom=buffer, top=1 - buffer)  # Positions of the subplots edges, as a fraction of the figure width.

        ax_left = plt.subplot2grid(shape=(5, 2), loc=(0, 0), rowspan=4, colspan=1)
        ax_left.set_xlim(-p.xlim, p.xlim)
        ax_left.set_ylim(-p.ylim, p.ylim)
        ax_left.set_aspect('equal')
        ax_left.set_facecolor("black")  # background color
        # ax_left.get_xaxis().set_visible(False)
        # ax_left.get_yaxis().set_visible(False)

        ax_right = plt.subplot2grid(shape=(5, 2), loc=(0, 1), rowspan=4, colspan=1)
        ax_right.set_xlim(-p.xlim, p.xlim)
        ax_right.set_ylim(-p.ylim, p.ylim)
        ax_right.set_aspect('equal')
        ax_right.set_facecolor("black")  # background color

        ax_lightcurve = plt.subplot2grid(shape=(5, 1), loc=(4, 0), rowspan=1, colspan=1)
        ax_lightcurve.set_facecolor("black")  # background color

        ax_lightcurve.tick_params(axis='x', colors='grey')
        xmax = p.iterations * p.dt / p.x_unit_value
        ax_lightcurve.set_xlim(0, xmax)
        xvalues = [x * CurveSimAnimation.tic_delta(xmax) for x in range(round(xmax / CurveSimAnimation.tic_delta(xmax)))]
        xlabels = [f'{round(x + p.start_date, 4)} {p.x_unit_name}' for x in xvalues]
        ax_lightcurve.set_xticks(xvalues, labels=xlabels)

        ax_lightcurve.tick_params(axis='y', colors='grey')
        minl = lightcurve.min(initial=None)
        maxl = lightcurve.max(initial=None)
        if minl == maxl:
            minl *= 0.99
        scope = maxl - minl
        buffer = 0.05 * scope
        ax_lightcurve.set_ylim(minl - buffer, maxl + buffer)

        ticdelta = CurveSimAnimation.tic_delta(maxl - minl)
        yvalues = [1 - y * ticdelta for y in range(round(float((maxl - minl) / ticdelta)))]
        ylabels = [f'{round(100 * y, 10)} %' for y in yvalues]
        ax_lightcurve.set_yticks(yvalues, labels=ylabels)

        time_axis = np.arange(0, round(p.iterations * p.dt), round(p.sampling_rate * p.dt), dtype=float)
        time_axis /= p.x_unit_value
        ax_lightcurve.plot(time_axis, sampled_lightcurve[0:len(time_axis)], color="white")

        red_dot = matplotlib.patches.Ellipse((0, 0), p.iterations * p.dt * p.red_dot_width / p.x_unit_value, scope * p.red_dot_height)  # matplotlib patch
        red_dot.set(zorder=2)  # Dot in front of lightcurve.
        red_dot.set_color((1, 0, 0))  # red
        ax_lightcurve.add_patch(red_dot)
        plt.tight_layout()  # Automatically adjust padding horizontally as well as vertically.
        return fig, ax_right, ax_left, ax_lightcurve, red_dot

    @staticmethod
    def next_frame_old(frame, p, bodies, red_dot, lightcurve):
        """Update patches. Send new circle positions to animation function.
        First parameter comes from iterator frames (a parameter of FuncAnimation).
        The other parameters are given to this function via the parameter fargs of FuncAnimation."""
        for body in bodies:  # left view: projection (x,y,z) -> (x,z), order = -y (y-axis points away from viewer)
            body.circle_left.set(zorder=-body.positions[frame * p.sampling_rate][1])
            body.circle_left.center = body.positions[frame * p.sampling_rate][0] / p.scope_left, body.positions[frame * p.sampling_rate][2] / p.scope_left
        for body in bodies:  # right view: projection (x,y,z) -> (x,-y), order = -z (z-axis points away from viewer)
            body.circle_right.set(zorder=-body.positions[frame * p.sampling_rate][2])
            body.circle_right.center = body.positions[frame * p.sampling_rate][0] / p.scope_right, -body.positions[frame * p.sampling_rate][1] / p.scope_right
        red_dot.center = p.dt * p.sampling_rate * frame / p.x_unit_value, lightcurve[frame * p.sampling_rate]
        if frame >= 10 and frame % int(round(p.frames / 10)) == 0:  # Inform user about program's progress.
            print(f'{round(frame / p.frames * 10) * 10:3d}% ', end="")

    @staticmethod
    def next_frame(frame, p, bodies, red_dot, lightcurve):
        """Update patches. Send new circle positions to animation function.
        First parameter comes from iterator frames (a parameter of FuncAnimation).
        The other parameters are given to this function via the parameter fargs of FuncAnimation."""
        for body in bodies:  # left view (overhead): projection (x,y,z) -> (z,-x), order = -y (y-axis points away from viewer)
            body.circle_left.set(zorder=-body.positions[frame * p.sampling_rate][1])
            body.circle_left.center = body.positions[frame * p.sampling_rate][2] / p.scope_left, -body.positions[frame * p.sampling_rate][0] / p.scope_left
        for body in bodies:  # right view (edge-on): projection (x,y,z) -> (z,-y), order = x (x-axis points to viewer)
            body.circle_right.set(zorder=body.positions[frame * p.sampling_rate][0])
            body.circle_right.center = body.positions[frame * p.sampling_rate][2] / p.scope_right, -body.positions[frame * p.sampling_rate][1] / p.scope_right
        red_dot.center = p.dt * p.sampling_rate * frame / p.x_unit_value, lightcurve[frame * p.sampling_rate]
        if frame >= 10 and frame % int(round(p.frames / 10)) == 0:  # Inform user about program's progress.
            print(f'{round(frame / p.frames * 10) * 10:3d}% ', end="")

    def render(self, p, bodies, lightcurve):
        """Calls next_frame() for each frame and saves the video."""
        print(f'Animating {p.frames:8d} frames:     ', end="")
        tic = time.perf_counter()
        anim = matplotlib.animation.FuncAnimation(self.fig, CurveSimAnimation.next_frame, fargs=(p, bodies, self.red_dot, lightcurve), interval=1000 / p.fps, frames=p.frames, blit=False)
        anim.save(p.video_file, fps=p.fps, metadata={"title": " "}, extra_args=['-vcodec', 'libx264'])  # https://www.ffmpeg.org/libavcodec.html
        toc = time.perf_counter()
        print(f' {toc - tic:7.2f} seconds  ({p.frames / (toc - tic):.0f} frames/second)')
        print(f'{p.video_file} saved.')
