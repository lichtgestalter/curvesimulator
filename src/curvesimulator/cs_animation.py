from colorama import Fore, Style
import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import patches
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np
import shutil
import subprocess
import sys
import time

class CurveSimAnimation:

    def __init__(self, p, bodies, sim_rv, sim_flux, time_s0):
        CurveSimAnimation.check_ffmpeg()  # is FFmpeg installed?
        self.fig, ax_right, ax_left, ax_lightcurve, self.rv_dot, self.flux_dot = CurveSimAnimation.init_plot(p, sim_rv, sim_flux, time_s0)  # Adjust constants in section [Plot] of config file to fit your screen.
        # rv_dot/flux_dot change every frame. Mark them "animated" so the initial full draw
        # (used to cache the static background for blitting) does NOT bake them in.
        if self.rv_dot is not None:
            self.rv_dot.set_animated(True)
        if self.flux_dot is not None:
            self.flux_dot.set_animated(True)
        for body in bodies:  # Circles represent the bodies in the animation. Set their colors and add them to the matplotlib axis.
            if body.image_file_left is None or body.image_file_right is None:
                body.circle_left.set_color(body.color)
                body.circle_right.set_color(body.color)
                if p.show_left_plot:
                    body.circle_left.set_animated(True)  # Moves every frame -> exclude from cached background.
                    ax_left.add_patch(body.circle_left)
                if p.show_right_plot:
                    body.circle_right.set_animated(True)
                    ax_right.add_patch(body.circle_right)
            else:
                body.image_left = OffsetImage(mpimg.imread(body.image_file_left), zoom=1.0)
                body.image_right = OffsetImage(mpimg.imread(body.image_file_right), zoom=1.0)
                body.ab_left = AnnotationBbox(body.image_left, (0.2, 0.2), frameon=False, xycoords='data')
                body.ab_left.set_animated(True)
                ax_left.add_artist(body.ab_left)
                body.ab_right = AnnotationBbox(body.image_right, (-0.5, -0.5), frameon=False, xycoords='data')
                body.ab_right.set_animated(True)
                ax_right.add_artist(body.ab_right)
        self.render(p, bodies, sim_rv, sim_flux, time_s0)

    @staticmethod
    def check_ffmpeg():
        """Checks if ffmpeg is in PATH"""
        if shutil.which("ffmpeg") is None:
            print(f"{Fore.RED}\nERROR: FFmpeg is not available. Please install FFmpeg to save the video.")
            print("Visit ffmpeg.org to download an executable version.")
            print(f"On Windows, extract the zip file and add the bin directory to your system's PATH environment variable.{Style.RESET_ALL}")
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
    def init_left_plot(p, shape, loc, rowspan, colspan):
        # left plot (overhead view)
        ax_left = plt.subplot2grid(shape=shape, loc=loc, rowspan=rowspan, colspan=colspan)
        ax_left.set_xlim(-p.xlim, p.xlim)
        ax_left.set_ylim(-p.ylim, p.ylim)
        ax_left.set_aspect("equal")
        ax_left.set_facecolor("black")  # background color
        ax_left.set_title(p.left_title, color="xkcd:light grey", fontsize=10, y=0.9)

        scale_bar_end_x = p.xlim * 0.90
        scale_bar_start_x = scale_bar_end_x - p.scale_bar_length_left / p.scope_left
        dy = p.ylim * 0.02
        scale_bar_height = p.ylim * -0.94
        scale_bar_text_height = p.ylim * -0.97
        scale_bar_text_start_x = (scale_bar_start_x + scale_bar_end_x) / 2
        ax_left.hlines(y=scale_bar_height, xmin=scale_bar_start_x, xmax=scale_bar_end_x, color='grey', linewidth=1)
        ax_left.vlines(x=scale_bar_start_x, ymin=scale_bar_height - dy, ymax=scale_bar_height + dy, color='grey', linewidth=1)
        ax_left.vlines(x=scale_bar_end_x, ymin=scale_bar_height - dy, ymax=scale_bar_height + dy, color='grey', linewidth=1)
        ax_left.text(scale_bar_text_start_x, scale_bar_text_height, f"{p.scale_bar_length_left / p.au:.2f} AU", color='grey', fontsize=8, ha='center', va='top')
        return ax_left

    @staticmethod
    def init_right_plot(p, shape, loc, rowspan, colspan):
        # right plot (edge-on view)
        ax_right = plt.subplot2grid(shape=shape, loc=loc, rowspan=rowspan, colspan=colspan)
        ax_right.set_xlim(-p.xlim, p.xlim)
        ax_right.set_ylim(-p.ylim, p.ylim)
        ax_right.set_aspect("equal")
        ax_right.set_facecolor("black")  # background color
        ax_right.set_title(p.right_title, color="xkcd:light grey", fontsize=10, y=0.9)

        scale_bar_end_x = p.xlim * 0.99
        scale_bar_start_x = scale_bar_end_x - p.scale_bar_length_right / p.scope_right
        dy = p.ylim * 0.02
        scale_bar_height = p.ylim * -0.94
        scale_bar_text_height = p.ylim * -0.97
        scale_bar_text_start_x = (scale_bar_start_x + scale_bar_end_x) / 2
        ax_right.hlines(y=scale_bar_height, xmin=scale_bar_start_x, xmax=scale_bar_end_x, color='grey', linewidth=1)
        ax_right.vlines(x=scale_bar_start_x, ymin=scale_bar_height - dy, ymax=scale_bar_height + dy, color='grey', linewidth=1)
        ax_right.vlines(x=scale_bar_end_x, ymin=scale_bar_height - dy, ymax=scale_bar_height + dy, color='grey', linewidth=1)
        ax_right.text(scale_bar_text_start_x, scale_bar_text_height, f"{p.scale_bar_length_right / p.au:.2f} AU", color='grey', fontsize=8, ha='center', va='top')
        return ax_right

    @staticmethod
    def init_lightcurve_plot(sim_flux, time_s0, p, shape, loc, rowspan, colspan):
        ax_lightcurve = plt.subplot2grid(shape=shape, loc=loc, rowspan=rowspan, colspan=colspan)
        ax_lightcurve.set_facecolor("black")  # background color

        # no x-tics/-labels because the rv-plot below uses the same x-tics/-labels

        # lightcurve y-ticks, y-labels
        ax_lightcurve.set_ylabel("Relative Flux", color="grey", fontsize=8)
        ax_lightcurve.tick_params(axis="y", colors="grey", labelsize=8)
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
        ylabels = [f"{round(100 * y, 10):.{digits}f} %" for y in yvalues]
        ax_lightcurve.set_yticks(yvalues, labels=ylabels)

        # lightcurve data (white line)
        x = (time_s0 - p.starts_s0[0]) / p.day
        ax_lightcurve.set_xlim(float(x[0]), float(x[-1]))
        ax_lightcurve.plot(x, sim_flux, color="white")

        # lightcurve red dot
        flux_dot = patches.Ellipse((0, 0), (time_s0[-1] - time_s0[0]) * p.flux_dot_width / p.day, scope * p.flux_dot_height)  # matplotlib patch
        flux_dot.set(zorder=2)  # Dot in front of lightcurve.
        flux_dot.set_color((1, 0, 0))  # red
        ax_lightcurve.add_patch(flux_dot)
        return ax_lightcurve, flux_dot

    @staticmethod
    def init_rv_curve_plot(sim_rv, time_s0, p, shape, loc, rowspan, colspan):
        ax_rv_curve = plt.subplot2grid(shape=shape, loc=loc, rowspan=rowspan, colspan=colspan)
        ax_rv_curve.set_facecolor("black")  # background color
        ax_rv_curve.text(0.97, -0.11, "BJD (TDB)", color="grey", fontsize=8, ha="right", va="bottom", transform=ax_rv_curve.transAxes)

        # rv_curve x-ticks, x-labels
        ax_rv_curve.tick_params(axis="x", colors="grey")
        # Use the same relative x-axis as the lightcurve: days since p.starts_s0[0]
        x = (time_s0 - p.starts_s0[0]) / p.day
        x_listticdelta = CurveSimAnimation.tic_delta(float(x[-1]))
        digits = max(0, round(-math.log10(x_listticdelta) + 0.4))  # The labels get as many decimal places as the intervals between the tics.
        # build tick positions in relative days and corresponding absolute-time labels (BJD)
        n_ticks = max(1, int(round(float(x[-1]) / x_listticdelta)))
        xvalues = [i * x_listticdelta for i in range(n_ticks + 1)]
        xlabels = [f"{round(val + p.epoch + p.starts_s0[0] / p.day, 4):.{digits}f}" for val in xvalues]
        ax_rv_curve.set_xticks(xvalues, labels=xlabels)
        ax_rv_curve.set_xlim(float(x[0]), float(x[-1]))

        # rv_curve y-ticks, y-labels
        ax_rv_curve.set_ylabel("Radial Velocity [m/s]", color="grey", labelpad=22, fontsize=8)
        ax_rv_curve.tick_params(axis="y", colors="grey", labelsize=8)
        minl = sim_rv.min(initial=None)
        maxl = sim_rv.max(initial=None)
        if minl == maxl:
            minl *= 0.99
        scope = maxl - minl
        buffer = 0.05 * scope
        ax_rv_curve.set_ylim(minl - buffer, maxl + buffer)
        y_listticdelta = CurveSimAnimation.tic_delta(scope)
        digits = max(0, round(-math.log10(y_listticdelta) + 0.4) - 2)  # The labels get as many decimal places as the intervals between the tics.
        if maxl > 0 > minl:
            yvalues = [maxl, 0, minl]
        else:
            yvalues = [maxl, minl]
        # yvalues = [maxl - y * y_listticdelta for y in range(round(float((maxl - minl) / y_listticdelta)))]
        ylabels = [f"{round(1 * y, 10):.{digits}f}" for y in yvalues]
        ax_rv_curve.set_yticks(yvalues, labels=ylabels)

        # rv_curve data (white line) using relative days since p.starts_s0[0]
        x = (time_s0 - p.starts_s0[0]) / p.day
        ax_rv_curve.plot(x, sim_rv, color="white")

        # rv_curve green dot
        rv_dot = patches.Ellipse((0, 0), (time_s0[-1] - time_s0[0]) * p.rv_dot_width / p.day, scope * p.rv_dot_height)  # matplotlib patch
        rv_dot.set(zorder=2)  # Dot in front of rv_curve.
        rv_dot.set_color((0, 0.7, 0))  # green
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

        if p.show_left_plot and p.show_right_plot and p.show_lc_plot and p.show_rv_plot:
            ax_left = CurveSimAnimation.init_left_plot(p, shape=(6, 2), loc=(0, 0), rowspan=4, colspan=1)
            ax_right = CurveSimAnimation.init_right_plot(p, shape=(6, 2), loc=(0, 1), rowspan=4, colspan=1)
            ax_lightcurve, flux_dot = CurveSimAnimation.init_lightcurve_plot(sim_flux, time_s0, p, shape=(6, 2), loc=(4, 0), rowspan=1, colspan=2)
            ax_rv_curve, rv_dot = CurveSimAnimation.init_rv_curve_plot(sim_rv, time_s0, p, shape=(6, 2), loc=(5, 0), rowspan=1, colspan=2)
            fig.add_artist(plt.Line2D([0.5, 0.5], [0.4, 0.9], color='xkcd:dark grey', linewidth=1, transform=fig.transFigure))
        # no RV plot
        elif p.show_left_plot and p.show_right_plot and p.show_lc_plot and not p.show_rv_plot:
            ax_left = CurveSimAnimation.init_left_plot(p, shape=(5, 2), loc=(0, 0), rowspan=4, colspan=1)
            ax_right = CurveSimAnimation.init_right_plot(p, shape=(5, 2), loc=(0, 1), rowspan=4, colspan=1)
            ax_lightcurve, flux_dot = CurveSimAnimation.init_lightcurve_plot(sim_flux, time_s0, p, shape=(5, 2), loc=(4, 0), rowspan=1, colspan=2)
            ax_rv_curve, rv_dot = None, None
            fig.add_artist(plt.Line2D([0.5, 0.5], [0.27, 0.9], color='xkcd:dark grey', linewidth=1, transform=fig.transFigure))

        # no light curve plot
        elif p.show_left_plot and p.show_right_plot and not p.show_lc_plot and p.show_rv_plot:
            ax_left = CurveSimAnimation.init_left_plot(p, shape=(5, 2), loc=(0, 0), rowspan=4, colspan=1)
            ax_right = CurveSimAnimation.init_right_plot(p, shape=(5, 2), loc=(0, 1), rowspan=4, colspan=1)
            ax_lightcurve, flux_dot = None, None
            ax_rv_curve, rv_dot = CurveSimAnimation.init_rv_curve_plot(sim_rv, time_s0, p, shape=(5, 2), loc=(4, 0), rowspan=1, colspan=2)
            fig.add_artist(plt.Line2D([0.5, 0.5], [0.27, 0.9], color='xkcd:dark grey', linewidth=1, transform=fig.transFigure))
        # no light curve and no RV plot
        elif p.show_left_plot and p.show_right_plot and not p.show_lc_plot and not p.show_rv_plot:
            ax_left = CurveSimAnimation.init_left_plot(p, shape=(4, 2), loc=(0, 0), rowspan=4, colspan=1)
            ax_right = CurveSimAnimation.init_right_plot(p, shape=(4, 2), loc=(0, 1), rowspan=4, colspan=1)
            ax_lightcurve, flux_dot = None, None
            ax_rv_curve, rv_dot = None, None
            fig.add_artist(plt.Line2D([0.5, 0.5], [0.07, 0.9], color='xkcd:dark grey', linewidth=1, transform=fig.transFigure))

        # right plot only
        elif not p.show_left_plot and p.show_right_plot and not p.show_lc_plot and not p.show_rv_plot:
            ax_left = None
            ax_right = CurveSimAnimation.init_right_plot(p, shape=(1, 1), loc=(0, 0), rowspan=1, colspan=1)
            ax_lightcurve, flux_dot = None, None
            ax_rv_curve, rv_dot = None, None
        # right plot + light curve
        elif not p.show_left_plot and p.show_right_plot and p.show_lc_plot and not p.show_rv_plot:
            ax_left = None
            ax_right = CurveSimAnimation.init_right_plot(p, shape=(5, 1), loc=(0, 0), rowspan=4, colspan=1)
            ax_lightcurve, flux_dot = CurveSimAnimation.init_lightcurve_plot(sim_flux, time_s0, p, shape=(5, 1), loc=(4, 0), rowspan=1, colspan=1)
            ax_rv_curve, rv_dot = None, None
        # right plot + rv curve
        elif not p.show_left_plot and p.show_right_plot and not p.show_lc_plot and p.show_rv_plot:
            ax_left = None
            ax_right = CurveSimAnimation.init_right_plot(p, shape=(5, 1), loc=(0, 0), rowspan=4, colspan=1)
            ax_lightcurve, flux_dot = None, None
            ax_rv_curve, rv_dot = CurveSimAnimation.init_rv_curve_plot(sim_rv, time_s0, p, shape=(5, 1), loc=(4, 0), rowspan=1, colspan=1)
        # right plot + light curve + rv curve
        elif not p.show_left_plot and p.show_right_plot and p.show_lc_plot and p.show_rv_plot:
            ax_left = None
            ax_right = CurveSimAnimation.init_right_plot(p, shape=(6, 1), loc=(0, 0), rowspan=4, colspan=1)
            ax_lightcurve, flux_dot = CurveSimAnimation.init_lightcurve_plot(sim_flux, time_s0, p, shape=(6, 1), loc=(4, 0), rowspan=1, colspan=1)
            ax_rv_curve, rv_dot = CurveSimAnimation.init_rv_curve_plot(sim_rv, time_s0, p, shape=(6, 1), loc=(5, 0), rowspan=1, colspan=1)

        # left plot only
        elif p.show_left_plot and not p.show_right_plot and not p.show_lc_plot and not p.show_rv_plot:
            ax_left = CurveSimAnimation.init_left_plot(p, shape=(1, 1), loc=(0, 0), rowspan=1, colspan=1)
            ax_right = None
            ax_lightcurve, flux_dot = None, None
            ax_rv_curve, rv_dot = None, None
        # left plot + light curve
        elif p.show_left_plot and not p.show_right_plot and p.show_lc_plot and not p.show_rv_plot:
            ax_left = CurveSimAnimation.init_left_plot(p, shape=(5, 1), loc=(0, 0), rowspan=4, colspan=1)
            ax_right = None
            ax_lightcurve, flux_dot = CurveSimAnimation.init_lightcurve_plot(sim_flux, time_s0, p, shape=(5, 1), loc=(4, 0), rowspan=1, colspan=1)
            ax_rv_curve, rv_dot = None, None
        # left plot + rv curve
        elif p.show_left_plot and not p.show_right_plot and not p.show_lc_plot and p.show_rv_plot:
            ax_left = CurveSimAnimation.init_left_plot(p, shape=(5, 1), loc=(0, 0), rowspan=4, colspan=1)
            ax_right = None
            ax_lightcurve, flux_dot = None, None
            ax_rv_curve, rv_dot = CurveSimAnimation.init_rv_curve_plot(sim_rv, time_s0, p, shape=(5, 1), loc=(4, 0), rowspan=1, colspan=1)
        # left plot + light curve + rv curve
        elif p.show_left_plot and not p.show_right_plot and p.show_lc_plot and p.show_rv_plot:
            ax_left = CurveSimAnimation.init_left_plot(p, shape=(6, 1), loc=(0, 0), rowspan=4, colspan=1)
            ax_right = None
            ax_lightcurve, flux_dot = CurveSimAnimation.init_lightcurve_plot(sim_flux, time_s0, p, shape=(6, 1), loc=(4, 0), rowspan=1, colspan=1)
            ax_rv_curve, rv_dot = CurveSimAnimation.init_rv_curve_plot(sim_rv, time_s0, p, shape=(6, 1), loc=(5, 0), rowspan=1, colspan=1)

        # light curve only
        elif not p.show_left_plot and not p.show_right_plot and p.show_lc_plot and not p.show_rv_plot:
            ax_left = None
            ax_right = None
            ax_lightcurve, flux_dot = CurveSimAnimation.init_lightcurve_plot(sim_flux, time_s0, p, shape=(1, 1), loc=(0, 0), rowspan=1, colspan=1)
            ax_rv_curve, rv_dot = None, None
        # rv curve only
        elif not p.show_left_plot and not p.show_right_plot and not p.show_lc_plot and p.show_rv_plot:
            ax_left = None
            ax_right = None
            ax_lightcurve, flux_dot = None, None
            ax_rv_curve, rv_dot = CurveSimAnimation.init_rv_curve_plot(sim_rv, time_s0, p, shape=(1, 1), loc=(0, 0), rowspan=1, colspan=1)
        # light curve + rv curve
        elif not p.show_left_plot and not p.show_right_plot and p.show_lc_plot and p.show_rv_plot:
            ax_left = None
            ax_right = None
            ax_lightcurve, flux_dot = CurveSimAnimation.init_lightcurve_plot(sim_flux, time_s0, p, shape=(2, 1), loc=(0, 0), rowspan=1, colspan=1)
            ax_rv_curve, rv_dot = CurveSimAnimation.init_rv_curve_plot(sim_rv, time_s0, p, shape=(2, 1), loc=(1, 0), rowspan=1, colspan=1)

        else:
            print(f"{Fore.RED}\nERROR: No plot was chosen to be displayed in the video.")
            print("Check the settings of parameters show_left_plot, show_right_plot, show_lc_plot and show_rv_plot.{Style.RESET_ALL}")
            sys.exit(1)

        plt.tight_layout()  # Automatically adjust padding horizontally as well as vertically.
        plt.suptitle(p.main_title, color="white", fontsize=14)
        fig.text(0.99, 0.99, "lichtgestalter/CurveSimulator", color="xkcd:purpley", fontsize=10, ha="right", va="top", transform=fig.transFigure)

        return fig, ax_right, ax_left, ax_lightcurve, rv_dot, flux_dot

    @staticmethod
    def next_frame(frame, p, bodies, rv_dot, flux_dot, sim_rv, sim_flux, time_s0):
        # p.clockwise = False
        if p.clockwise:
            x_direction = 1
        else:
            x_direction = -1
        frame_number = int(frame * p.sampling_rate)
        """Update patches. Send new circle positions to animation function.
        First parameter comes from iterator frames (a parameter of FuncAnimation).
        The other parameters are given to this function via the parameter fargs of FuncAnimation."""
        for body in bodies:  # left view: projection (x,y,z) -> (x,-z), order = y (y-axis points to viewer)
            if body.image_file_left is None or body.image_file_right is None:
                body.circle_left.set(zorder=body.positions[frame_number][1])
                body.circle_left.center = x_direction * body.positions[frame_number][0] / p.scope_left, -body.positions[frame_number][2] / p.scope_left
            else:
                body.ab_left.set_zorder(body.positions[frame_number][1])
                body.ab_left.xybox = (x_direction * body.positions[frame_number][0] / p.scope_left, -body.positions[frame_number][2] / p.scope_left)
        for body in bodies:  # right view: projection (x,y,z) -> (x,y), order = z (z-axis points to viewer)
            if body.image_file_left is None or body.image_file_right is None:
                body.circle_right.set(zorder=body.positions[frame_number][2])
                body.circle_right.center = x_direction * body.positions[frame_number][0] / p.scope_right, body.positions[frame_number][1] / p.scope_right
            else:
                body.ab_right.set_zorder(body.positions[frame_number][2])
                body.ab_right.xybox = (x_direction * body.positions[frame_number][0] / p.scope_right, body.positions[frame_number][1] / p.scope_right)

        # Use relative x (days since p.starts_s0[0]) for both dots so they align with plotted curves
        x_rel = (time_s0[frame_number] - p.starts_s0[0]) / p.day
        if p.show_lc_plot:
            flux_dot.center = x_rel, sim_flux[frame_number]
        if p.show_rv_plot:
            rv_dot.center = x_rel, sim_rv[frame_number]
        # if frame > 10:
        #     bodies[0].circle_left.set_color((1.0, 0.2, 0.2))  # Example code for changing circle color during animation
        if frame >= 10 and frame % int(round(p.frames / 10)) == 0:  # Inform user about program"s progress.
            print(f"{round(frame / p.frames * 10) * 10:3d}% ", end="")

        # Collect the artists that changed this frame. render() uses these to blit
        # only the changed regions instead of redrawing the whole figure.
        artists = [flux_dot] if flux_dot else []
        if rv_dot:
            artists.append(rv_dot)

        # Add AnnotationBboxes (image mode) or circles (non-image mode).
        for body in bodies:
            if hasattr(body, 'ab_left') and body.ab_left:
                artists.append(body.ab_left)
            elif hasattr(body, 'circle_left') and p.show_left_plot:
                artists.append(body.circle_left)
            if hasattr(body, 'ab_right') and body.ab_right:
                artists.append(body.ab_right)
            elif hasattr(body, 'circle_right') and p.show_right_plot:
                artists.append(body.circle_right)

        return artists

    def render(self, p, bodies, sim_rv, sim_flux, time_s0):
        """Calls next_frame() for each frame and saves the video.

        Note: matplotlib's animation.Animation.save() hard-codes blit=False for every frame
        it saves (see its source: `anim._draw_next_frame(d, blit=False)`), so blit=True on
        FuncAnimation has NO effect at all when the animation is written to a video file via
        anim.save() - blit only ever helps interactive on-screen animation (plt.show()).
        To actually benefit from blitting (only redraw the artists that changed instead of
        the whole figure every frame) we drive the rendering manually here and pipe the raw
        pixel data straight into FFmpeg, instead of going through FuncAnimation/anim.save().
        """
        frames = int(len(sim_flux) // p.sampling_rate)
        if p.verbose:
            print(f"Animating {p.frames:8d} frames:     ", end="")
            tic = time.perf_counter()

        fig = self.fig
        canvas = fig.canvas
        canvas.draw()  # Full initial draw of everything NOT marked animated -> the static background.
        renderer = canvas.get_renderer()
        background = canvas.copy_from_bbox(fig.bbox)
        buf = np.asarray(renderer.buffer_rgba())
        height, width = buf.shape[0], buf.shape[1]
        # libx264 with yuv420p output requires even width/height. Odd pixel dimensions
        # (which can happen depending on figure_width/figure_height/dpi) make FFmpeg exit
        # immediately, which in turn causes a BrokenPipeError on the very first stdin.write().
        width -= width % 2
        height -= height % 2
        ffmpeg_cmd = [
            "ffmpeg",
            "-loglevel", "error",   # only print errors (also suppresses the version/configuration banner)
            "-y",
            "-f", "rawvideo",
            "-vcodec", "rawvideo",
            "-s", f"{width}x{height}",
            "-pix_fmt", "rgba",
            "-r", str(p.fps),
            "-i", "-",
            "-an",
            "-vcodec", "libx264",
            "-crf", "18",  # Constant Rate Factor (lower value means better quality)
            "-preset", "slow",  # Preset for better compression
            "-b:v", "30000k",  # Bitrate
            "-pix_fmt", "yuv420p",
            str(p.video_file),
        ]
        proc = subprocess.Popen(ffmpeg_cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

        try:
            for frame in range(frames):
                artists = CurveSimAnimation.next_frame(frame, p, bodies, self.rv_dot, self.flux_dot, sim_rv, sim_flux, time_s0)
                canvas.restore_region(background)
                # Group by axes and sort by zorder, then draw only the changed artists on top of the cached background.
                by_axes = {}
                for artist in artists:
                    if artist is None:
                        continue
                    by_axes.setdefault(artist.axes, []).append(artist)
                for ax, arts in by_axes.items():
                    arts.sort(key=lambda a: a.get_zorder())
                    for a in arts:
                        if ax is not None:
                            ax.draw_artist(a)
                        else:
                            fig.draw_artist(a)
                canvas.blit(fig.bbox)
                frame_buf = np.asarray(renderer.buffer_rgba())[:height, :width]
                try:
                    proc.stdin.write(frame_buf.tobytes())
                except BrokenPipeError:
                    break  # FFmpeg died; stop feeding it and report its stderr below.
        finally:
            try:
                proc.stdin.close()
            except OSError:
                pass
            stderr_output = proc.stderr.read()
            proc.stderr.close()
            proc.wait()
            if proc.returncode != 0:
                print(f"{Fore.RED}\nERROR: FFmpeg exited with code {proc.returncode} while writing {p.video_file}.{Style.RESET_ALL}")
                if stderr_output:
                    print(f"{Fore.RED}FFmpeg output:\n{stderr_output.decode(errors='replace')}{Style.RESET_ALL}")

        if p.verbose:
            toc = time.perf_counter()
            print(f" {toc - tic:7.2f} seconds  ({p.frames / (toc - tic):.0f} frames/second)")
            print(f"{p.video_file} saved.")
