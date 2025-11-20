from colorama import Fore, Style
import math
import matplotlib
import matplotlib.pyplot as plt  # from matplotlib import pyplot as plt
import numpy as np
import shutil
import sys
import time

class CurveSimAnimation:

    def __init__(self, p, bodies, sim_rv, sim_flux, time_s0):
        CurveSimAnimation.check_ffmpeg()
        self.fig, ax_right, ax_left, ax_lightcurve, self.flux_dot = CurveSimAnimation.init_plot(p, sim_rv, sim_flux, time_s0)  # Adjust constants in section [Plot] of config file to fit your screen.
        for body in bodies:  # Circles represent the bodies in the animation. Set their colors and add them to the matplotlib axis.
            body.circle_right.set_color(body.color)
            body.circle_left.set_color(body.color)
            ax_right.add_patch(body.circle_right)
            ax_left.add_patch(body.circle_left)
        self.render(p, bodies, sim_rv, sim_flux, time_s0)

    @staticmethod
    def check_ffmpeg():
        """Checks if ffmpeg is available"""
        if shutil.which('ffmpeg') is None:
            print(f"{Fore.RED}ERROR: ffmpeg is not available. Please install ffmpeg to save the video.")
            print("Visit ffmpeg.org to download an executable version.")
            print(f"Extract the zip file and add the bin directory to your system's PATH environment variable.{Style.RESET_ALL}")
            sys.exit(1)

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
    def init_left_plot(p):
        # left plot (overhead view)
        ax_left = plt.subplot2grid(shape=(5, 2), loc=(0, 0), rowspan=4, colspan=1)
        ax_left.set_xlim(-p.xlim, p.xlim)
        ax_left.set_ylim(-p.ylim, p.ylim)
        ax_left.set_aspect('equal')
        ax_left.set_facecolor("black")  # background color
        return ax_left

    @staticmethod
    def init_right_plot(p):
        # right plot (edge-on view)
        ax_right = plt.subplot2grid(shape=(5, 2), loc=(0, 1), rowspan=4, colspan=1)
        ax_right.set_xlim(-p.xlim, p.xlim)
        ax_right.set_ylim(-p.ylim, p.ylim)
        ax_right.set_aspect('equal')
        ax_right.set_facecolor("black")  # background color
        ax_right.text(.99, .99, "lichtgestalter/CurveSimulator", color='grey', fontsize=10, ha='right', va='top', transform=ax_right.transAxes)
        return ax_right

    @staticmethod
    def init_lightcurve_plot(sim_flux, time_s0, p):
        # lightcurve
        ax_lightcurve = plt.subplot2grid(shape=(5, 1), loc=(4, 0), rowspan=1, colspan=1)
        ax_lightcurve.set_facecolor("black")  # background color
        ax_lightcurve.text(1.00, -0.05, "BJD (TDB)", color='grey', fontsize=10, ha='right', va='bottom', transform=ax_lightcurve.transAxes)

        # lightcurve x-ticks, x-labels
        ax_lightcurve.tick_params(axis='x', colors='grey')
        xmax = time_s0[-1] / p.day
        ax_lightcurve.set_xlim(0, xmax)
        x_listticdelta = CurveSimAnimation.tic_delta(xmax)
        digits = max(0, round(-math.log10(x_listticdelta) + 0.4))  # The labels get as many decimal places as the intervals between the tics.
        xvalues = [x * x_listticdelta for x in range(round(xmax / x_listticdelta))]
        xlabels = [f'{round(x + p.start_date, 4):.{digits}f}' for x in xvalues]
        ax_lightcurve.set_xticks(xvalues, labels=xlabels)

        # lightcurve y-ticks, y-labels
        ax_lightcurve.tick_params(axis='y', colors='grey')
        minl = sim_flux.min(initial=None)
        maxl = sim_flux.max(initial=None)
        if minl == maxl:
            minl *= 0.99
        scope = maxl - minl
        buffer = 0.05 * scope
        ax_lightcurve.set_ylim(minl - buffer, maxl + buffer)
        y_listticdelta = CurveSimAnimation.tic_delta(scope)
        digits = max(0, round(-math.log10(y_listticdelta) + 0.4) - 2)  # The labels get as many decimal places as the intervals between the tics.
        yvalues = [1 - y * y_listticdelta for y in range(round(float((maxl - minl) / y_listticdelta)))]
        ylabels = [f'{round(100 * y, 10):.{digits}f} %' for y in yvalues]
        ax_lightcurve.set_yticks(yvalues, labels=ylabels)

        # lightcurve data (white line)
        ax_lightcurve.plot(time_s0 / p.day, sim_flux, color="white")

        # lightcurve red dot
        flux_dot = matplotlib.patches.Ellipse((0, 0), (time_s0[-1] - time_s0[0]) * p.flux_dot_width / p.day, scope * p.flux_dot_height)  # matplotlib patch
        flux_dot.set(zorder=2)  # Dot in front of lightcurve.
        flux_dot.set_color((1, 0, 0))  # red
        ax_lightcurve.add_patch(flux_dot)
        return ax_lightcurve, flux_dot

    @staticmethod
    def init_rv_curve_plot(sim_rv, time_s0, p):
        # rv_curve
        ax_rv_curve = plt.subplot2grid(shape=(5, 1), loc=(4, 0), rowspan=1, colspan=1)
        ax_rv_curve.set_facecolor("black")  # background color
        ax_rv_curve.text(1.00, -0.05, "BJD (TDB)", color='grey', fontsize=10, ha='right', va='bottom', transform=ax_rv_curve.transAxes)

        # rv_curve x-ticks, x-labels
        ax_rv_curve.tick_params(axis='x', colors='grey')
        xmax = time_s0[-1] / p.day
        ax_rv_curve.set_xlim(0, xmax)
        x_listticdelta = CurveSimAnimation.tic_delta(xmax)
        digits = max(0, round(-math.log10(x_listticdelta) + 0.4))  # The labels get as many decimal places as the intervals between the tics.
        xvalues = [x * x_listticdelta for x in range(round(xmax / x_listticdelta))]
        xlabels = [f'{round(x + p.start_date, 4):.{digits}f}' for x in xvalues]
        ax_rv_curve.set_xticks(xvalues, labels=xlabels)

        # rv_curve y-ticks, y-labels
        ax_rv_curve.tick_params(axis='y', colors='grey')
        minl = sim_rv.min(initial=None)
        maxl = sim_rv.max(initial=None)
        if minl == maxl:
            minl *= 0.99
        scope = maxl - minl
        buffer = 0.05 * scope
        ax_rv_curve.set_ylim(minl - buffer, maxl + buffer)
        y_listticdelta = CurveSimAnimation.tic_delta(scope)
        digits = max(0, round(-math.log10(y_listticdelta) + 0.4) - 2)  # The labels get as many decimal places as the intervals between the tics.
        yvalues = [1 - y * y_listticdelta for y in range(round(float((maxl - minl) / y_listticdelta)))]
        ylabels = [f'{round(100 * y, 10):.{digits}f} %' for y in yvalues]
        ax_rv_curve.set_yticks(yvalues, labels=ylabels)

        # rv_curve data (white line)
        ax_rv_curve.plot(time_s0 / p.day, sim_rv, color="white")

        # rv_curve red dot
        rv_dot = matplotlib.patches.Ellipse((0, 0), (time_s0[-1] - time_s0[0]) * p.rv_dot_width / p.day, scope * p.rv_dot_height)  # matplotlib patch
        rv_dot.set(zorder=2)  # Dot in front of rv_curve.
        rv_dot.set_color((1, 0, 0))  # red
        ax_rv_curve.add_patch(rv_dot)
        return ax_rv_curve, rv_dot

    @staticmethod
    def init_plot(p, sim_rv, sim_flux, time_s0):
        """Initialize the matplotlib figure containing 3 axis:
        Top left: overhead view
        Top right: edge-on view
        Bottom: lightcurve and rv_curve"""
        fig = plt.figure()
        fig.set_figwidth(p.figure_width)
        fig.set_figheight(p.figure_height)
        fig.set_facecolor("black")  # background color outside of ax_left and ax_lightcurve
        buffer = 0
        fig.subplots_adjust(left=buffer, right=1.0 - buffer, bottom=buffer, top=1 - buffer)  # Positions of the subplots edges, as a fraction of the figure width.

        ax_left = CurveSimAnimation.init_left_plot(p)
        ax_right = CurveSimAnimation.init_right_plot(p)
        ax_lightcurve, flux_dot = CurveSimAnimation.init_lightcurve_plot(sim_flux, time_s0, p)
        ax_rv_curve, rv_dot = CurveSimAnimation.init_rv_curve_plot(sim_rv, time_s0, p)
        plt.tight_layout()  # Automatically adjust padding horizontally as well as vertically.

        return fig, ax_right, ax_left, ax_lightcurve, flux_dot

    @staticmethod
    def next_frame(frame, p, bodies, flux_dot, sim_rv, sim_flux, time_s0):
        """Update patches. Send new circle positions to animation function.
        First parameter comes from iterator frames (a parameter of FuncAnimation).
        The other parameters are given to this function via the parameter fargs of FuncAnimation."""
        # for body in bodies:  # left view: projection (x,y,z) -> (z,x), order = y (y-axis points to viewer)
        #     body.circle_left.set(zorder=body.positions[frame * p.sampling_rate][1])
        #     body.circle_left.center = body.positions[frame * p.sampling_rate][2] / p.scope_left, body.positions[frame * p.sampling_rate][0] / p.scope_left
        for body in bodies:  # left view: projection (x,y,z) -> (x,-z), order = y (y-axis points to viewer)
            body.circle_left.set(zorder=body.positions[frame * p.sampling_rate][1])
            body.circle_left.center = body.positions[frame * p.sampling_rate][0] / p.scope_left, -body.positions[frame * p.sampling_rate][2] / p.scope_left
        # for body in bodies:  # right view: projection (x,y,z) -> (z,y), order = -x (x-axis points away from viewer)
        #     body.circle_right.set(zorder=-body.positions[frame * p.sampling_rate][0])
        #     body.circle_right.center = body.positions[frame * p.sampling_rate][2] / p.scope_right, body.positions[frame * p.sampling_rate][1] / p.scope_right
        for body in bodies:  # right view: projection (x,y,z) -> (x,y), order = z (z-axis points to viewer)
            body.circle_right.set(zorder=body.positions[frame * p.sampling_rate][2])
            body.circle_right.center = body.positions[frame * p.sampling_rate][0] / p.scope_right, body.positions[frame * p.sampling_rate][1] / p.scope_right
        flux_dot.center = time_s0[frame * p.sampling_rate] / p.day, sim_flux[frame * p.sampling_rate]
        if frame >= 10 and frame % int(round(p.frames / 10)) == 0:  # Inform user about program's progress.
            print(f'{round(frame / p.frames * 10) * 10:3d}% ', end="")

    def render(self, p, bodies, sim_rv, sim_flux, time_s0):
        """Calls next_frame() for each frame and saves the video."""
        if p.verbose:
            print(f'Animating {p.frames:8d} frames:     ', end="")
            tic = time.perf_counter()
        frames = len(sim_flux) // p.sampling_rate
        anim = matplotlib.animation.FuncAnimation(self.fig, CurveSimAnimation.next_frame, fargs=(p, bodies, self.flux_dot, sim_rv, sim_flux, time_s0), interval=1000 / p.fps, frames=frames, blit=False)
        anim.save(
            p.video_file,
            fps=p.fps,
            metadata={"title": " "},
            extra_args=[
                '-vcodec', 'libx264',
                '-crf', '18',  # Constant Rate Factor (lower value means better quality)
                '-preset', 'slow',  # Preset for better compression
                '-b:v', '30000k'  # Bitrate 5000k (increase as needed)
            ]
        )
        if p.verbose:
            toc = time.perf_counter()
            print(f' {toc - tic:7.2f} seconds  ({p.frames / (toc - tic):.0f} frames/second)')
            print(f'{p.video_file} saved.')
