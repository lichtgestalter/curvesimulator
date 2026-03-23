from colorama import Fore, Style
import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import patches, animation
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np
import shutil
import sys
import time
import ast
import configparser
import json
import matplotlib
import matplotlib.animation
import rebound
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, StringVar
import copy
import corner
import emcee
import emcee.autocorr
from functools import wraps
import lmfit
from multiprocessing import Pool
import os
import msvcrt
import random
import bisect
import pandas as pd
import re
import scipy.stats as stats
from multiprocessing import Process, JoinableQueue
import warnings


class CurveSimAnimation:

    def __init__(self, p, bodies, sim_rv, sim_flux, time_s0):
        CurveSimAnimation.check_ffmpeg()  # is ffmpeg installed?
        self.fig, ax_right, ax_left, ax_lightcurve, self.rv_dot, self.flux_dot = CurveSimAnimation.init_plot(p, sim_rv, sim_flux, time_s0)  # Adjust constants in section [Plot] of config file to fit your screen.
        for body in bodies:  # Circles represent the bodies in the animation. Set their colors and add them to the matplotlib axis.
            if body.image_file_left is None or body.image_file_right is None:
                body.circle_left.set_color(body.color)
                body.circle_right.set_color(body.color)
                if p.show_left_plot:
                    ax_left.add_patch(body.circle_left)
                if p.show_right_plot:
                    ax_right.add_patch(body.circle_right)
            else:
                body.image_left = OffsetImage(mpimg.imread(body.image_file_left), zoom=1.0)
                body.image_right = OffsetImage(mpimg.imread(body.image_file_right), zoom=1.0)
                body.ab_left = AnnotationBbox(body.image_left, (0.2, 0.2), frameon=False, xycoords='data')
                ax_left.add_artist(body.ab_left)
                body.ab_right = AnnotationBbox(body.image_right, (-0.5, -0.5), frameon=False, xycoords='data')
                ax_right.add_artist(body.ab_right)
        self.render(p, bodies, sim_rv, sim_flux, time_s0)

    @staticmethod
    def check_ffmpeg():
        """Checks if ffmpeg is in PATH"""
        if shutil.which("ffmpeg") is None:
            print(f"{Fore.RED}\nERROR: ffmpeg is not available. Please install ffmpeg to save the video.")
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
        # lightcurve
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
        # rv_curve
        ax_rv_curve = plt.subplot2grid(shape=shape, loc=loc, rowspan=rowspan, colspan=colspan)
        ax_rv_curve.set_facecolor("black")  # background color
        ax_rv_curve.text(1.00, -0.15, "BJD (TDB)", color="grey", fontsize=10, ha="right", va="bottom", transform=ax_rv_curve.transAxes)

        # rv_curve x-ticks, x-labels
        ax_rv_curve.tick_params(axis="x", colors="grey")
        # Use the same relative x-axis as the lightcurve: days since p.starts_s0[0]
        x = (time_s0 - p.starts_s0[0]) / p.day
        x_listticdelta = CurveSimAnimation.tic_delta(float(x[-1]))
        digits = max(0, round(-math.log10(x_listticdelta) + 0.4))  # The labels get as many decimal places as the intervals between the tics.
        # build tick positions in relative days and corresponding absolute-time labels (BJD)
        n_ticks = max(1, int(round(float(x[-1]) / x_listticdelta)))
        xvalues = [i * x_listticdelta for i in range(n_ticks + 1)]
        xlabels = [f"{round(val + p.start_date + p.starts_s0[0] / p.day, 4):.{digits}f}" for val in xvalues]
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

        # Return artists for compatibility with FuncAnimation (ignored if blit=False)
        # Alle geänderten Artists sammeln
        artists = [flux_dot] if flux_dot else []
        if rv_dot:
            artists.append(rv_dot)

        # AnnotationBbox und Circles hinzufügen
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
        """Calls next_frame() for each frame and saves the video."""
        frames = int(len(sim_flux) // p.sampling_rate)
        if p.verbose:
            print(f"Animating {p.frames:8d} frames:     ", end="")
            tic = time.perf_counter()
        anim = animation.FuncAnimation(self.fig, CurveSimAnimation.next_frame, fargs=(p, bodies, self.rv_dot, self.flux_dot, sim_rv, sim_flux, time_s0), interval=1000 / p.fps, frames=frames, blit=False)
        # anim = animation.FuncAnimation(self.fig, CurveSimAnimation.next_frame_counterclockwise, fargs=(p, bodies, self.rv_dot, self.flux_dot, sim_rv, sim_flux, time_s0), interval=1000 / p.fps, frames=frames, blit=False)
        anim.save(
            p.video_file,
            fps=p.fps,
            metadata={"title": " "},
            extra_args=[
                "-vcodec", "libx264",
                "-crf", "18",  # Constant Rate Factor (lower value means better quality)
                "-preset", "slow",  # Preset for better compression
                "-b:v", "30000k"  # Bitrate 5000k (increase as needed)
            ]
        )
        if p.verbose:
            toc = time.perf_counter()
            print(f" {toc - tic:7.2f} seconds  ({p.frames / (toc - tic):.0f} frames/second)")
            print(f"{p.video_file} saved.")
# import pandas

# from curvesimulator.cs_parameters import CurveSimParameters


class CurveSimBodies(list):

    # noinspection PyUnusedLocal
    def __init__(self, p):
        p.myintegration = False  # debug
        """Initialize instances of physical bodies.
        Read program parameters and properties of the bodies from config file.
        Initialize the circles in the animation (matplotlib patches)"""
        # For ease of use of these constants in the config file are additionally defined here without the prefix "p.".
        super().__init__()  # Call the superclass initializer
        try:
            g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
            r_jup, m_jup, r_earth, m_earth, v_earth = p.r_jup, p.m_jup, p.r_earth, p.m_earth, p.v_earth
            hour, day, year = p.hour, p.day, p.year

        except AttributeError:
            print(f"{Fore.YELLOW}\nWARNING: Section <Astronomical Constants> in the configuration file is incomplete.")
            print(f"See https://github.com/lichtgestalter/curvesimulator/wiki.{Style.RESET_ALL}")
        config = configparser.ConfigParser(inline_comment_prefixes="#")
        config.optionxform = str  # Preserve case of the keys.
        config.read(p.config_file)  # Read config file. (This time the physical objects.)

        # Physical bodies
        for section in config.sections():
            if section not in p.standard_sections:  # section describes a physical object
                file = config.get(section, "file", fallback=None)
                if file is None:
                    kwargs = {
                        "p": p,
                        "primary": config.get(section, "primary", fallback=None),
                        "name": section,
                        "body_type": config.get(section, "body_type", fallback=None),
                        "color": tuple([ast.literal_eval(x) for x in config.get(section, "color", fallback="-1").split(",")]),
                        "image_file_left": config.get(section, "image_file_left", fallback=None),
                        "image_file_right": config.get(section, "image_file_right", fallback=None),
                        "mass": p.read_param(config, section, "mass", fallback="-1"),
                        "radius": p.read_param(config, section, "radius", fallback="-1"),
                        "luminosity": p.read_param(config, section, "luminosity", fallback="0.0"),
                        "limb_darkening_1": p.read_param(config, section, "limb_darkening_1", fallback="None"),
                        "limb_darkening_2": p.read_param(config, section, "limb_darkening_2", fallback="None"),
                        "limb_darkening_parameter_type": config.get(section, "limb_darkening_parameter_type", fallback=None),
                        "startposition": config.get(section, "startposition", fallback=None),
                        "velocity": config.get(section, "velocity", fallback=None),
                        "e": p.read_param(config, section, "e", fallback="-1"),
                        "i": p.read_param(config, section, "i", fallback="-1111"),
                        "P": p.read_param(config, section, "P", fallback="None"),
                        "a": p.read_param(config, section, "a", fallback="None"),
                        "Omega": p.read_param(config, section, "Omega", fallback="None"),
                        "omega": p.read_param(config, section, "omega", fallback="None"),
                        "pomega": p.read_param(config, section, "pomega", fallback="None"),
                        "L": p.read_param(config, section, "L", fallback="None"),
                        "ma": p.read_param(config, section, "ma", fallback="None"),
                        "ea": p.read_param(config, section, "ea", fallback="None"),
                        "nu": p.read_param(config, section, "nu", fallback="None"),
                        "T": p.read_param(config, section, "T", fallback="None"),
                        "t": p.read_param(config, section, "t", fallback="0.0"),
                    }
                    body = CurveSimBody(**kwargs)
                else:
                    body = CurveSimBody.load(file, p)
                self.append(body)
        self.check_body_parameters()
        p.bodynames2bodies(self)
        if p.action == "single_run":
            self.generate_patches(p)

    def __repr__(self):
        names = "CurveSimBodies: "
        for body in self:
            names += body.name + ", "
        return names[:-2]

    def init_rebound(self, p):
        simulation = rebound.Simulation()
        simulation.G = p.g  # gravitational constant
        star_count = sum(1 for body in self if body.body_type == "star")
        if star_count == 1:
            # simulation.integrator = "ias15"
            simulation.integrator = "whfast"
            simulation.dt = p.dt
        if p.verbose:
            print(f"using Rebound integrator {simulation.integrator}:", end="")
        if star_count == 0:
            print(f"{Fore.RED}\nERROR: No body in config file has body type star.{Style.RESET_ALL}")
            sys.exit(1)

        for body in self[0:1]:  # hack debug: works only when the first body is the only star and all other bodies are orbiting this star (no binary, no moons, ...)
            simulation.add(m=body.mass, r=body.radius, hash=body.name)

        for i, body in enumerate(self[1:], start=1):
            kwargs = {}
            kwargs["m"] = body.mass
            kwargs["r"] = body.radius
            kwargs["hash"] = body.name
            kwargs["inc"] = body.i
            kwargs["e"] = body.e
            if body.P is not None:
                kwargs["P"] = body.P
            if body.a is not None:
                kwargs["a"] = body.a
            if body.pomega is not None:
                kwargs["pomega"] = body.pomega
            else:
                if body.Omega is not None:
                    kwargs["Omega"] = body.Omega
                if body.omega is not None:
                    kwargs["omega"] = body.omega
            if body.ma is not None:
                kwargs["M"] = body.ma
            if body.nu is not None:
                kwargs["f"] = body.nu
            if body.ea is not None:
                kwargs["E"] = body.ea
            if body.T is not None:
                kwargs["T"] = body.T
            if body.L is not None:
                kwargs["l"] = body.L

            # primary = CurveSimBodies.get_com_particle(simulation, range(i))
            # kwargs["primary"] = primary

            simulation.add(**kwargs)
        simulation.move_to_com()  # move origin to center of mass before integrating -> better numerical stability
        CurveSimBodies.print_simulation_particles(simulation)

        # if p.action == "single_run":  # obsolete????  does not seem to help for MCMC, but is a good choice when creating a result file including transit times
        #     if p.result_file:
        #         simulation.ri_whfast.safe_mode = 0  # see https://rebound.readthedocs.io/en/latest/ipython_examples/AdvWHFast/
        #         simulation.ri_whfast.corrector = 11  # hopefully more accurate
        return simulation

    @staticmethod
    def get_com_particle(simulation, indices):
        """Returns a temporary particle representing the center of mass
        of the given particle indices."""
        m_tot = 0.0
        x = y = z = 0.0
        vx = vy = vz = 0.0

        for i in indices:
            p = simulation.particles[i]
            m_tot += p.m
            x += p.m * p.x
            y += p.m * p.y
            z += p.m * p.z
            vx += p.m * p.vx
            vy += p.m * p.vy
            vz += p.m * p.vz

        com = rebound.Particle()
        com.m = m_tot
        com.x = x / m_tot
        com.y = y / m_tot
        com.z = z / m_tot
        com.vx = vx / m_tot
        com.vy = vy / m_tot
        com.vz = vz / m_tot

        return com

    @staticmethod
    def print_simulation_particles(simulation):
        print("\n--- Simulation Particles ---")
        for i, particle in enumerate(simulation.particles):
            print(f"\nParticle {i}:")
            print(f"  hash      = {particle.hash}")
            print(f"  mass (m)  = {particle.m}")
            print(f"  radius (r)= {particle.r}")
            print(f"  position  = ({particle.x}, {particle.y}, {particle.z})")
            print(f"  velocity  = ({particle.vx}, {particle.vy}, {particle.vz})")
            try:
                orbit = particle.orbit()
                print(f"Orbital elements (relative to primary if available)")
                print(f"  P         = {orbit.P / (60 * 60 * 24)}")
                print(f"  a         = {orbit.a}")
                print(f"  e         = {orbit.e}")
                print(f"  i         = {math.degrees(orbit.inc)}")
                print(f"  Omega     = {math.degrees(orbit.Omega)}")
                print(f"  omega     = {math.degrees(orbit.omega)}")
                print(f"  ma        = {math.degrees(orbit.M)}")
            except Exception:
                print("  (No orbital elements available)")
            try:
                primary = CurveSimBodies.get_com_particle(simulation, range(i))
                orbit = particle.orbit(primary=primary)
                print(f"Orbital elements (relative to manually computed primary if available)")
                print(f"  P         = {orbit.P / (60 * 60 * 24)}")
                print(f"  a         = {orbit.a}")
                print(f"  e         = {orbit.e}")
                print(f"  i         = {math.degrees(orbit.inc)}")
                print(f"  Omega     = {math.degrees(orbit.Omega)}")
                print(f"  omega     = {math.degrees(orbit.omega)}")
                print(f"  ma        = {math.degrees(orbit.M)}")
            except Exception:
                print("  (No orbital elements available)")
            try:
                primary = CurveSimBodies.get_com_particle(simulation, range(i))
                orbit = particle.orbit(primary=simulation.particles[0])
                print(f"Orbital elements (relative to star if available)")
                print(f"  P         = {orbit.P / (60 * 60 * 24)}")
                print(f"  a         = {orbit.a}")
                print(f"  e         = {orbit.e}")
                print(f"  i         = {math.degrees(orbit.inc)}")
                print(f"  Omega     = {math.degrees(orbit.Omega)}")
                print(f"  omega     = {math.degrees(orbit.omega)}")
                print(f"  ma        = {math.degrees(orbit.M)}")
            except Exception:
                print("  (No orbital elements available)")

    def check_body_parameters(self):
        """Checking parameters of physical bodies in the config file"""
        if len(self) == 0:
            print(f"{Fore.RED}\nERROR in config file: No physical bodies have been specified.")
            sys.exit(1)
        if len(self) == 1:
            print(f"{Fore.RED}\nERROR in config file: Just one physical body has been specified.")
            sys.exit(1)
        for body in self:
            if body.radius <= 0:
                print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid or missing radius.")
                sys.exit(1)
            if body.mass <= 0:
                print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid or missing mass.")
                sys.exit(1)
            if body.luminosity < 0:
                print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid luminosity {body.luminosity=}.")
                sys.exit(1)
            if body.luminosity > 0 and (body.limb_darkening_u1 is None or body.limb_darkening_u2 is None):  # if body.luminosity > 0 and limb darkening parameters are missing
                print(f"{Fore.RED}\nERROR in config file: {body.name} has luminosity but invalid limb darkening parameter {body.limb_darkening=}.")
                sys.exit(1)
            for c in body.color:
                if c < 0 or c > 1 or len(body.color) != 3:
                    print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid or missing color value.")
                    sys.exit(1)
            # if body.velocity is None:
            #     if body.e < 0:
            #         print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid or missing eccentricity e.")
            #         sys.exit(1)
            #     if body.i < -1000:
            #         print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid or missing inclination i.")
            #         sys.exit(1)
            if body.a is not None and body.P is not None:
                print(f"{Fore.RED}\nERROR in config file: Period P and semi-major axis a have been specified for {body.name}.")
                print(f"{Fore.RED}Remove one parameter from the config file.")
                sys.exit(1)
            if body.a is not None and body.a <= 0:
                print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid semi-major axis a.")
                sys.exit(1)
            if body.P is not None and body.P <= 0:
                print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid period P.")
                sys.exit(1)
            anomaly_counter = 0
            anomalies = [body.L, body.ma, body.ea, body.nu, body.T]
            for anomaly in anomalies:
                if anomaly is not None:
                    anomaly_counter += 1
            if anomaly_counter > 1:
                print(f"{Fore.YELLOW}\nWARNING: more than one anomaly (L, ma, ea, nu, T) has been specified in config file for {body.name}.")
                print(f"Check for contradictions and/or remove superflous anomalies.{Style.RESET_ALL}")

    def save(self, prefix="", suffix=""):
        for body in self:
            body.save(prefix, suffix)

    # def calc_primary_body_initial_velocity(self):
    #     """Calculates the initial velocity of the primary body in the star system
    #         from the masses and initial velocities of all other bodies.
    #         The calculation is based on the principles of conservation of momentum
    #         and the center of mass motion"""
    #     assert 0 == self[0].velocity[0] == self[0].velocity[1] == self[0].velocity[2]
    #     for body in self[1:]:
    #         self[0].velocity += body.velocity * body.mass
    #     self[0].velocity /= - self[0].mass

    def total_luminosity(self, stars, iteration, p):
        """Add luminosity of all stars in the system while checking for eclipses.
        Does not yet work correctly for eclipsed eclipses (three or more bodies in line of sight at the same time)."""
        luminosity = 0.0
        for star in stars:
            luminosity += star.luminosity
            for body in self:
                if body != star:  # an object cannot eclipse itself :)
                    eclipsed_area, relative_radius = star.eclipsed_by(body, iteration, p)
                    if eclipsed_area is not None:
                        absolute_depth = star.intensity * eclipsed_area * CurveSimPhysics.limbdarkening(relative_radius, star.limb_darkening_u1, star.limb_darkening_u2) / star.mean_intensity
                        luminosity -= absolute_depth
                        # results["Bodies"][body.name]["Transits"][-1]["impacts_and_depths"][-1].depth = absolute_depth  # this depth is caused by this particular body eclipsing this particular star
        return luminosity

    # @staticmethod
    # def distance_and_direction(body1, body2, iteration, p):
    #     # Calculate distance and direction between 2 bodies:
    #     distance_xyz = body2.positions[iteration - 1] - body1.positions[iteration - 1]
    #     distance = math.sqrt(np.dot(distance_xyz, distance_xyz))
    #     force_total = p.g * body1.mass * body2.mass / distance ** 2  # Use law of gravitation to calculate force acting on body.
    #     x, y, z = distance_xyz[0], distance_xyz[1], distance_xyz[2]
    #     polar_angle = math.acos(z / distance)
    #     azimuth_angle = math.atan2(y, x)
    #     return force_total, azimuth_angle, polar_angle
    #
    # @staticmethod
    # def update_force(azimuth_angle, force, force_total, polar_angle):
    #     # Compute the force of attraction in each direction:
    #     force[0] += math.sin(polar_angle) * math.cos(azimuth_angle) * force_total
    #     force[1] += math.sin(polar_angle) * math.sin(azimuth_angle) * force_total
    #     force[2] += math.cos(polar_angle) * force_total
    #
    # @staticmethod
    # def update_velocity(body1, iteration, force, p):
    #     """https://en.wikipedia.org/wiki/Verlet_integration
    #     https://www.lancaster.ac.uk/staff/drummonn/PHYS281/gravity/"""
    #     if iteration == 1:
    #         body1.acceleration = force / body1.mass
    #     acceleration = force / body1.mass
    #     body1.velocity += (acceleration + body1.acceleration) * 0.5 * p.dt
    #     body1.acceleration = acceleration
    #     return acceleration

    # @staticmethod
    # def update_velocity_euler(body1, force, p):
    #     acceleration = force / body1.mass
    #     body1.velocity += acceleration * p.dt
    #     return acceleration

    # @staticmethod
    # def update_position(body1, iteration, acceleration, p):
    #     """https://en.wikipedia.org/wiki/Verlet_integration
    #     https://www.lancaster.ac.uk/staff/drummonn/PHYS281/gravity/
    #     Verlet integration avoids the numerical problems of the Euler method."""
    #     movement = body1.velocity * p.dt + acceleration * (p.dt ** 2 * 0.5)
    #     body1.positions[iteration] = body1.positions[iteration - 1] + movement

    @staticmethod
    def update_position(body, iteration, rebound_sim):
        particle = rebound_sim.particles[body.name]
        body.positions[iteration] = np.array([particle.x, particle.y, particle.z])

    # @staticmethod
    # def update_position_euler(body1, iteration, acceleration, p):
    #     movement = body1.velocity * p.dt - 0.5 * acceleration * p.dt ** 2
    #     body1.positions[iteration] = body1.positions[iteration - 1] + movement

    @staticmethod
    def progress_bar(iteration, p):
        if p.total_iterations > 5:  # prevent DIV/0 in next line
            if iteration % int(round(p.total_iterations / 10)) == 0:  # Inform user about program"s progress.
                print(f"{round(iteration / p.total_iterations * 10) * 10:3d}% ", end="")
                # print(self.energy(iteration, p))

    def calc_positions_eclipses_luminosity(self, p, time_s0):
        """Calculate distances, forces, accelerations, velocities of the bodies for each iteration.
        The resulting body positions and the lightcurve are stored for later use in the animation."""

        if p.myintegration:  # debug
            simulation = MyIntegration(p, self)
            CurveSimBodies.init_myintegration(self, simulation)
        else:
            simulation = CurveSimBodies.init_rebound(self, p)

        stars = [body for body in self if body.body_type == "star"]
        sim_flux = CurveSimLightcurve(p.total_iterations)  # Initialize lightcurve (essentially a np.ndarray)
        sim_rv = np.full(p.total_iterations, np.nan, dtype=float)
        if not p.myintegration:
            initial_sim_state = CurveSimRebound(simulation)

        for iteration in range(p.total_iterations):
            if p.myintegration:
                if iteration == 0:
                    E0 = simulation.total_energy()
                E = simulation.total_energy()
                rel_error = (E - E0) / abs(E0)
                if iteration % (p.total_iterations // 10) == 0:
                    print(f"Energy drift: {rel_error:.2e}")

            simulation.integrate(time_s0[iteration])
            for body in self:
                CurveSimBodies.update_position(body, iteration, simulation)
            sim_flux[iteration] = self.total_luminosity(stars, iteration, p)  # Update sim_flux.
            sim_rv[iteration] = -simulation.particles[p.rv_body].vz
            if p.verbose:
                CurveSimBodies.progress_bar(iteration, p)
        if not p.myintegration:
            new_sim_state = CurveSimRebound(simulation)
            energy_change = initial_sim_state.sim_check_deltas(new_sim_state)
        else:
            energy_change = None

        lightcurve_max = float(sim_flux.max(initial=None))
        sim_flux /= lightcurve_max  # Normalize flux.
        return sim_rv, sim_flux, self, simulation, energy_change

    def calc_physics(self, p, time_s0):
        """Calculate body positions and the resulting lightcurve."""
        if p.verbose:
            if p.video_file and p.flux_file is None:
                print(f"Generating {p.frames} frames for a {p.frames / p.fps:.0f} seconds long video.")
            print(f"Calculating {p.total_iterations:,} iterations ", end="")
            tic = time.perf_counter()
        sim_rv, sim_flux, bodies, rebound_sim, energy_change = self.calc_positions_eclipses_luminosity(p, time_s0)
        if p.verbose:
            toc = time.perf_counter()
            print(f" {toc - tic:7.3f} seconds  ({p.total_iterations / (toc - tic):.0f} iterations/second)")
            print(f"Log10 of the relative change of energy during simulation: {energy_change:.0f}")
            if energy_change > -6:
                print(f"{Fore.YELLOW}The energy must not change significantly! Consider using a smaller time step (dt).{Style.RESET_ALL}")
        return sim_rv, sim_flux, rebound_sim

    def calc_patch_radii(self, p):
        """If autoscaling is on, this function calculates the radii of the circles (matplotlib patches) of the animation."""
        logs = [math.log10(body.radius) for body in self]  # log10 of all radii
        radii_out = [(p.max_radius - p.min_radius) * (i - min(logs)) / (max(logs) - min(logs)) + p.min_radius for i in logs]  # linear transformation to match the desired minmum and maximum radii
        # print(f"patch radii:", end="  ")
        for body, radius in zip(self, radii_out):
            body.patch_radius = radius

    def generate_patches(self, p):
        """Generates the circles (matplotlib patches) of the animation."""
        if p.autoscaling:
            self.calc_patch_radii(p)
            for body in self:
                body.circle_right = matplotlib.patches.Circle(xy=(0, 0), radius=body.patch_radius)  # Matplotlib patch for right view
                body.circle_left = matplotlib.patches.Circle(xy=(0, 0), radius=body.patch_radius)  # Matplotlib patch for left view
        else:
            for body in self:
                if body.body_type == "planet":
                    extrascale_left, extrascale_right = p.planet_scale_left, p.planet_scale_right  # Scale radius in plot.
                else:
                    extrascale_left, extrascale_right = p.star_scale_left, p.star_scale_right  # It's a star. Scale radius in plot accordingly.
                body.circle_right = matplotlib.patches.Circle((0, 0), radius=body.radius * extrascale_right / p.scope_right)  # Matplotlib patch for right view
                body.circle_left = matplotlib.patches.Circle((0, 0), radius=body.radius * extrascale_left / p.scope_left)  # Matplotlib patch for left view

    # def energy(self, iteration, p):
    #     """Calculates the total energy in the system. This should be constant."""
    #     kinetic_energy = 0
    #     potential_energy = 0
    #     for i, body1 in enumerate(self):
    #         velocity_magnitude = np.linalg.norm(body1.velocity)
    #         kinetic_energy += 0.5 * body1.mass * velocity_magnitude ** 2
    #         for j, body2 in enumerate(self):
    #             if i < j:
    #                 distance = np.linalg.norm(body2.positions[iteration - 1] - body1.positions[iteration - 1])
    #                 if distance > 0:
    #                     potential_energy += body1.mass * body2.mass / distance
    #     return kinetic_energy - p.g * potential_energy

    def find_transits(self, rebound_sim, p, sim_flux, time_s0, time_d):
        print()
        rebound_sim.dt = p.result_dt
        results = CurveSimResults(self)
        for start_index, end_index, dt in zip(p.start_indices[:-1], p.start_indices[1:], p.dts):
            for i in range(start_index, end_index):
                for eclipser in p.eclipsers:
                    for eclipsee in p.eclipsees:
                        eclipser_before_eclipsee = eclipser.positions[i][2] > eclipsee.positions[i][2]
                        transit_between_iterations = (eclipser.positions[i][0] - eclipsee.positions[i][0]) * (eclipser.positions[i - 1][0] - eclipsee.positions[i - 1][0]) <= 0  # transit between i-1 and i?
                        if eclipser_before_eclipsee and transit_between_iterations:
                            tt, impact, depth, close_enough = eclipsee.find_tt(eclipser, i - 1, rebound_sim, p, sim_flux, time_s0, time_d, start_index, end_index, dt)
                            if close_enough:  # eclipser and eclipsee are close enough at actual TT
                                tt_s0 = rebound_sim.t
                                t1 = eclipsee.find_t1234(eclipser, tt_s0, i, rebound_sim, time_s0, start_index, end_index, p, transittimetype="T1")
                                t2 = eclipsee.find_t1234(eclipser, tt_s0, i, rebound_sim, time_s0, start_index, end_index, p, transittimetype="T2")
                                t3 = eclipsee.find_t1234(eclipser, tt_s0, i - 1, rebound_sim, time_s0, start_index, end_index, p, transittimetype="T3")
                                t4 = eclipsee.find_t1234(eclipser, tt_s0, i - 1, rebound_sim, time_s0, start_index, end_index, p, transittimetype="T4")
                                t12, t23, t34, t14 = CurveSimPhysics.calc_transit_intervals(t1, t2, t3, t4)
                                results["Bodies"][eclipser.name]["Transits"].append(Transit(eclipsee))
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["EclipsedBody"] = eclipsee.name
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T1"] = t1
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T2"] = t2
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["TT"] = tt
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T3"] = t3
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T4"] = t4
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T12"] = t12
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T23"] = t23
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T34"] = t34
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T14"] = t14
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["b"] = impact
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["depth"] = depth
        return results

    @staticmethod
    def find_tts(rebound_sim, p, sim_flux, time_s0, time_d):
        tts = []
        rebound_sim.dt = p.result_dt
        for start_index, end_index, dt in zip(p.start_indices[:-1], p.start_indices[1:], p.dts):
            for i in range(start_index, end_index):
                for eclipser in p.eclipsers:
                    for eclipsee in p.eclipsees:
                        eclipser_before_eclipsee = eclipser.positions[i][2] > eclipsee.positions[i][2]
                        transit_between_iterations = (eclipser.positions[i][0] - eclipsee.positions[i][0]) * (eclipser.positions[i - 1][0] - eclipsee.positions[i - 1][0]) <= 0  # transit between i-1 and i?
                        if eclipser_before_eclipsee and transit_between_iterations:
                            tt, b, depth, close_enough = eclipsee.find_tt(eclipser, i - 1, rebound_sim, p, sim_flux, time_s0, time_d, start_index, end_index, dt)
                            if close_enough:
                                tts.append([eclipser.name, eclipsee.name, tt])
        # convert tts into a pandas Dataframe with columns eclipser, eclipsee, tt
        return tts

    def bodies2param_json(self, measured_tt, p):
        result = {}
        result["max_delta"] = float(np.max(np.abs(measured_tt["delta"])))
        result["mean_delta"] = float(np.mean(np.abs(measured_tt["delta"])))
        for i, body in enumerate(self):
            result[body.name] = {}
            for key in p.PARAMS:
                attr = getattr(body, key)
                if attr is not None:
                    if key in p.scale:
                        scale = p.scale[key]
                    else:
                        scale = 1
                    result[body.name][key] = attr * scale
        line = json.dumps(result)
        return line



    def init_myintegration(self, simulation):
        bodies = self

        # --- helper: compute COM of already-added particles ---
        def get_com(indices):
            m_tot = 0.0
            r = np.zeros(3)
            v = np.zeros(3)

            for idx in indices:
                p = simulation._particles_list[idx]
                m_tot += p.m
                r += p.m * np.array([p.x, p.y, p.z])
                v += p.m * np.array([p.vx, p.vy, p.vz])

            return m_tot, r / m_tot, v / m_tot

        simulation._particles_list = []

        # --- first body (star) ---
        star = bodies[0]
        simulation.add(
            name=star.name,
            m=star.mass,
            r=star.radius,
            x=0, y=0, z=0,
            vx=0, vy=0, vz=0
        )

        # simulation._particles_list.append(simulation.particles[star.name])

        # --- remaining bodies ---
        for i, body in enumerate(bodies[1:], start=1):

            m_com, r_com, v_com = get_com(range(i))

            # gravitational parameter (Jacobi!)
            mu = simulation.G * (m_com + body.mass)

            # --- determine semi-major axis ---
            if body.a is not None:
                a = body.a
            elif body.P is not None:
                # Kepler's 3rd law with Jacobi mass
                a = ( (mu * body.P**2) / (4 * np.pi**2) )**(1/3)
            else:
                raise ValueError(f"Body {body.name}: need either a or P")

            e = body.e
            inc = body.i
            Omega = body.Omega or 0.0
            omega = body.omega or 0.0
            M = body.ma or 0.0

            # --- solve Kepler ---
            E = M
            for _ in range(50):
                E -= (E - e*np.sin(E) - M) / (1 - e*np.cos(E))

            # --- true anomaly ---
            nu = 2*np.arctan2(
                np.sqrt(1+e)*np.sin(E/2),
                np.sqrt(1-e)*np.cos(E/2)
            )

            r = a * (1 - e*np.cos(E))

            # orbital plane
            x_orb = r * np.cos(nu)
            y_orb = r * np.sin(nu)

            # velocity in orbital plane
            n = np.sqrt(mu / a**3)
            vx_orb = -a * n * np.sin(E) / (1 - e*np.cos(E))
            vy_orb =  a * n * np.sqrt(1 - e**2) * np.cos(E) / (1 - e*np.cos(E))

            # --- rotation matrices ---
            cosO, sinO = np.cos(Omega), np.sin(Omega)
            cosi, sini = np.cos(inc), np.sin(inc)
            cosw, sinw = np.cos(omega), np.sin(omega)

            def rotate(x, y):
                X = (cosO*cosw - sinO*sinw*cosi)*x + (-cosO*sinw - sinO*cosw*cosi)*y
                Y = (sinO*cosw + cosO*sinw*cosi)*x + (-sinO*sinw + cosO*cosw*cosi)*y
                Z = (sini*sinw)*x + (sini*cosw)*y
                return np.array([X, Y, Z])

            r_vec = rotate(x_orb, y_orb)
            v_vec = rotate(vx_orb, vy_orb)

            # --- convert from Jacobi → inertial ---
            r_vec += r_com
            v_vec += v_com

            simulation.add(
                name=body.name,
                m=body.mass,
                r=body.radius,
                x=r_vec[0], y=r_vec[1], z=r_vec[2],
                vx=v_vec[0], vy=v_vec[1], vz=v_vec[2]
            )

            # simulation._particles_list.append(simulation.particles[body.name])


class MyParticle:
    def __init__(self, name, m, r):
        self.name = name
        self.m = m
        self.r = r
        self.x = self.y = self.z = 0.0
        self.vx = self.vy = self.vz = 0.0


class MyIntegration:
    def __init__(self, p, bodies):
        self.G = p.g
        self.t = 0.0
        self.dt = p.dt
        self.particles = {}          # name → particle
        self._particles_list = []    # ordered list (needed for Jacobi)

        # central mass (Jacobi approx)
        self.star = bodies[0]
        self.M0 = self.star.mass

        # create particles
        for body in bodies:
            self.particles[body.name] = MyParticle(body.name, body.mass, body.radius)

        # store orbital elements
        self.orbits = {}
        for body in bodies[1:]:
            self.orbits[body.name] = {
                "a": body.a,
                "P": body.P,
                "e": body.e,
                "i": body.i,
                "Omega": body.Omega or 0.0,
                "omega": body.omega or 0.0,
                "M0": body.ma or 0.0,
                "t0": 0.0,
                "mu": self.G * (self.M0 + body.mass)
            }

    def add(self, name, m, r, x, y, z, vx, vy, vz):
        p = MyParticle(name, m, r)
        p.x, p.y, p.z = x, y, z
        p.vx, p.vy, p.vz = vx, vy, vz
        self.particles[name] = p
        self._particles_list.append(p)

    def integrate(self, t):
        self.t = t

        # star fixed at origin
        star = self.particles[self.star.name]
        star.x = star.y = star.z = 0.0
        star.vx = star.vy = star.vz = 0.0

        for name, orb in self.orbits.items():
            a = orb["a"]
            e = orb["e"]
            i = orb["i"]
            Omega = orb["Omega"]
            omega = orb["omega"]
            M0 = orb["M0"]
            P = orb["P"]
            mu = orb["mu"]

            # mean motion
            n = 2 * math.pi / P

            # mean anomaly
            M = M0 + n * (t - orb["t0"])

            # solve Kepler: E - e sin E = M
            E = self.solve_kepler(M, e)

            # true anomaly
            nu = 2 * math.atan2(
                math.sqrt(1 + e) * math.sin(E / 2),
                math.sqrt(1 - e) * math.cos(E / 2)
            )

            # distance
            r = a * (1 - e * math.cos(E))

            # orbital plane coords
            x_orb = r * math.cos(nu)
            y_orb = r * math.sin(nu)

            # rotation to inertial frame
            cosO, sinO = math.cos(Omega), math.sin(Omega)
            cosi, sini = math.cos(i), math.sin(i)
            cosw, sinw = math.cos(omega), math.sin(omega)

            x = (cosO * cosw - sinO * sinw * cosi) * x_orb + (-cosO * sinw - sinO * cosw * cosi) * y_orb
            y = (sinO * cosw + cosO * sinw * cosi) * x_orb + (-sinO * sinw + cosO * cosw * cosi) * y_orb
            z = (sini * sinw) * x_orb + (sini * cosw) * y_orb

            p = self.particles[name]
            p.x, p.y, p.z = x, y, z

            # velocity (approx)
            vx = -math.sin(E) * n * a / (1 - e * math.cos(E))
            vy = math.sqrt(1 - e**2) * math.cos(E) * n * a / (1 - e * math.cos(E))

            p.vx, p.vy, p.vz = vx, vy, 0.0

    def solve_kepler(self, M, e, tol=1e-10):
        E = M if e < 0.8 else math.pi
        for _ in range(50):
            dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
            E -= dE
            if abs(dE) < tol:
                break
        return E

    def total_energy(self):
        kinetic = 0.0
        potential = 0.0

        plist = list(self.particles.values())

        # kinetic
        for p in plist:
            v2 = p.vx**2 + p.vy**2 + p.vz**2
            kinetic += 0.5 * p.m * v2

        # potential
        for i in range(len(plist)):
            for j in range(i+1, len(plist)):
                pi = plist[i]
                pj = plist[j]

                dx = pj.x - pi.x
                dy = pj.y - pi.y
                dz = pj.z - pi.z

                r = np.sqrt(dx*dx + dy*dy + dz*dz + 1e-12)
                potential -= self.G * pi.m * pj.m / r

        return kinetic + potential

# from curvesimulator.cs_results import Transit
# from curvesimulator.cs_results import CurveSimResults


def multiple_transit_error():
    print(f"{Fore.RED}\nERROR: Ambiguous transit constellation.")
    print("CurveSimulator can not handle multiple synchronous transits correctly yet.")
    print(f"Please send your config file to CurveSimulator's developers.{Style.RESET_ALL}")
    sys.exit(1)


# noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
class CurveSimBody:
    def __init__(self, p, primary, name, body_type, mass, radius, luminosity, startposition, velocity, P, a, e, i, Omega, omega, pomega, L, ma, ea,
                 nu, T, t, limb_darkening_1, limb_darkening_2, limb_darkening_parameter_type, color, image_file_left, image_file_right):
        """Initialize instance of physical body."""
        # For ease of use of constants in the config file they are additionally defined here without the prefix "p.".
        g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = p.r_jup, p.m_jup, p.r_earth, p.m_earth, p.v_earth
        self.name = name  # name
        self.body_type = body_type  # "star" or "planet"
        self.primary = primary
        self.color = color  # (R, G, B)  each between 0 and 1
        self.image_file_left = image_file_left
        self.image_file_right = image_file_right
        self.mass = mass  # [kg]
        self.radius = radius  # [m]
        self.area_2d = math.pi * radius ** 2  # [m**2]
        self.luminosity = luminosity  # [W]
        self.limb_darkening_u1, self.limb_darkening_u2 = CurveSimPhysics.get_limbdarkening_parameters(limb_darkening_1, limb_darkening_2, limb_darkening_parameter_type)
        self.mean_intensity = CurveSimPhysics.calc_mean_intensity(self.limb_darkening_u1, self.limb_darkening_u2)
        self.intensity = luminosity / self.area_2d  # luminosity per (apparent) area [W/m**2]
        self.positions = np.ndarray((p.max_iterations[0], 3), dtype=float)
        # self.positions = np.ndarray((p.total_iterations, 3), dtype=float)

        self.e = e  # [1] eccentricity
        self.i = i  # [rad] inclination
        self.i_deg = None if i is None else math.degrees(i)  # [deg] inclination

        self.P = P  # [s] period
        self.a = a  # [m] semi-major axis

        self.Omega = Omega  # [rad] longitude of ascending node
        self.Omega_deg = None if Omega is None else math.degrees(Omega)  # [deg] longitude of ascending node
        self.omega = omega  # [rad] argument of periapsis
        self.omega_deg = None if omega is None else math.degrees(omega)  # [deg] argument of periapsis
        self.pomega = pomega  # [rad] longitude of periapsis
        self.pomega_deg = None if pomega is None else math.degrees(pomega)  # [deg] longitude of periapsis

        self.L = L  # [rad] mean longitude
        self.L_deg = None if L is None else math.degrees(L)  # [deg] mean longitude
        self.ma = ma  # [rad] mean anomaly
        self.ma_deg = None if ma is None else math.degrees(ma)  # [deg] mean anomaly
        self.ea = ea  # [rad] eccentric anomaly
        self.ea_deg = None if ea is None else math.degrees(ea)  # [deg] eccentric anomaly
        self.nu = nu  # [rad] true anomaly
        self.nu_deg = None if nu is None else math.degrees(nu)  # [deg] true anomaly. Per definition = 270° at the time of an exoplanet's primary transit.
        self.T = T  # [s] Time of periapsis
        self.t = t  # [s] optional additional time delta. For example time since last time of transit

        self.mu = None  # Gravitational Parameter. Depends on the masses of at least 2 bodies.

        # if not primary and startposition is not None and velocity is not None:  # State vectors are already in config file.
        #     pos = []
        #     for x in startposition.split(","):
        #         pos.append(eval(x))
        #     vel = []
        #     for x in velocity.split(","):
        #         vel.append(eval(x))
        #     if len(pos) != 3:
        #         print(f"{Fore.RED}\nERROR in config file: invalid or missing start position. {pos=}")
        #         sys.exit(1)
        #     if len(vel) != 3:
        #         print(f"{Fore.RED}\nERROR in config file: invalid or missing initial velocity. {vel=}")
        #         sys.exit(1)
        #     self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
        #     self.velocity = np.array(vel, dtype=float)  # [m/s]
        # elif primary:
        #     self.positions[0] = np.array([0.0, 0.0, 0.0], dtype=float)  # [m] initial position
        #     self.velocity = np.array([0.0, 0.0, 0.0], dtype=float)  # [m/s] initial velocity will be updated after all other state vectors have been calculated.
        # else:  # State vectors are not in config file. They will be calculated from Kepler orbit parameters later on after all bodies are initialized.
        #     self.velocity = None

        # Used for calculation of eclipsed area in function eclipsed_by.
        self.d, self.h, self.angle, self.eclipsed_area = 0.0, 0.0, 0.0, 0.0

    def __repr__(self):
        return f"CurveSimBody: {self.name}"
    #
    # def save2(self, p, prefix="", suffix=""):
    #     """Saves attributes of self in a file"""
    #     filename = prefix + self.name + suffix + ".bdy"
    #     print(filename)
    #     print(self.__dict__)
    #     exclude = ["positions", "velocity", "circle_left", "circle_right", "acceleration", "d", "h", "angle", "eclipsed_area", "patch_radius"]
    #     for key in self.__dict__.keys():
    #         if key not in exclude:
    #             pass
    #
    # def load2(filename):
    #     """Read attributes of body from a file"""
    #     body = CurveSimBody()
    #     return body

    def save(self, prefix="", suffix=""):
        """Saves attributes of self in a simple text file (key = repr(value)).
        Excludes large/non-serializable/derived attributes listed in `exclude`.
        """
        filename = "../bodies/" + prefix + self.name + suffix + ".bdy"
        exclude = ["positions", "velocity", "circle_left", "circle_right", "acceleration", "area_2d", "d", "h", "angle", "eclipsed_area", "patch_radius"]
        try:
            with open(filename, "w", encoding="utf-8") as f:
                for key, val in self.__dict__.items():
                    if key in exclude:
                        continue
                    if key.startswith("_"):
                        continue
                    # skip callables and methods
                    if callable(val):
                        continue
                    f.write(f"{key} = {repr(val)}\n")
        except Exception as e:
            print(f"Error saving body to {filename}: {e}")

    @staticmethod
    def load(filename, p):
        """Loads a body from `../bodies/<filename>.bdy`,
        assembles the constructor args in the original __init__ order,
        calls CurveSimBody(*args)"""
        path = "../bodies/" + filename + ".bdy"
        data = {}
        try:
            with open(path, "r", encoding="utf-8") as f:
                for raw in f:
                    line = raw.strip()
                    if not line or line.startswith("#") or "=" not in line:
                        continue
                    key, val_str = line.split("=", 1)
                    key = key.strip()
                    val_str = val_str.strip()
                    data[key] = ast.literal_eval(val_str)
        except Exception as e:
            print(f"Error loading body from {path}: {e}")
            return None

        if "limb_darkening_1" not in data.keys():
            data["limb_darkening_1"] = data.get("limb_darkening_u1")
            data["limb_darkening_2"] = data.get("limb_darkening_u2")
            data["limb_darkening_parameter_type"] = "u"

        missing = []
        for param in p.PARAMS:
            if param not in data.keys():
                data[param] = None
                missing.append(param)
        if missing:
            print(f"{Fore.YELLOW}\nWARNING: Missing parameters in {filename + ".bdy"}: {missing} {Style.RESET_ALL}")


        # Build args in the same order as __init__ signature:
        args = [
            p,
            data.get("primary", False),
            data.get("name"),
            data.get("body_type"),
            data.get("mass"),
            data.get("radius"),
            data.get("luminosity"),
            data.get("startposition"),
            data.get("velocity"),
            data.get("P"),
            data.get("a"),
            data.get("e"),
            data.get("i"),
            data.get("Omega"),
            data.get("omega"),
            data.get("pomega"),
            data.get("L"),
            data.get("ma"),
            data.get("ea"),
            data.get("nu"),
            data.get("T"),
            data.get("t"),
            data.get("limb_darkening_1"),
            data.get("limb_darkening_2"),
            data.get("limb_darkening_parameter_type"),
            data.get("color"),
            data.get("image_file_left"),
            data.get("image_file_right"),
        ]

        try:
            body = CurveSimBody(*args)
        except Exception as e:
            print(f"Error constructing CurveSimBody from {path}: {e}")
            return None

        return body

    # noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
    # def calc_orbit_angles(self):
    #     if self.omega is None:
    #         self.omega = self.pomega - self.Omega
    #     elif self.pomega is None:
    #         self.pomega = self.omega + self.Omega
    #     elif self.Omega is None:
    #         self.Omega = self.pomega - self.omega
    #     else:
    #         error = abs(self.omega - self.pomega + self.Omega)
    #         if error > 0.00001:
    #             print(f"{Fore.RED}\nERROR in config file, body {self.name}:")
    #             print(f"omega, pomega, Omega have been defined in the config file for this body.")
    #             print("This is redundant and in this case contradictory.")
    #             print("Remove one of these parameters from the config file or")
    #             print("make sure that omega - pomega + Omega = 0")
    #             sys.exit(1)
    #
    # def calc_period_or_semi_major_axis(self):
    #     if self.a is None and self.P is None:
    #         print(f"{Fore.RED}\nERROR in config file, body {self.name}:")
    #         print("semi-major axis a or Period P have to be specified in config file.")
    #         sys.exit(1)
    #     elif self.P is None:
    #         self.P = 2 * math.pi * math.sqrt(self.a ** 3 / self.mu)
    #     elif self.a is None:
    #         self.a = ((self.mu * self.P ** 2) / (4 * math.pi ** 2)) ** (1/3)
    #     else:
    #         relative_error = self.P / (2 * math.pi * math.sqrt(self.a ** 3 / self.mu)) - 1
    #         if relative_error > 0.001:
    #             print(f"{Fore.RED}\nERROR in config file, body {self.name}:")
    #             print(f"a and P have been defined in the config file for this body.")
    #             print("This is redundant and in this case contradictory.")
    #             print("Remove one of these parameters from the config file or")
    #             print("make sure that a and P are compatible with Kepler's third law.")
    #             sys.exit(1)
    #
    # def calc_anomalies(self):
    #     """[a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
    #        [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
    #        Numbers in comments refer to numbered formulas in [a] and [b]."""
    #
    #     a, e, L, pomega = self.a, self.e, self.L, self.pomega  # for readability of formulas
    #     ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas
    #
    #     if ma is None and L is not None:
    #         ma = L - pomega
    #         # print("Variant 1: ma-  pomega+  L+, calc ma")
    #     if ea is not None:  # ea provided
    #         nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
    #         ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
    #         # print("Variant 2: ea+, calc nu ma")
    #     else:  # ea not provided
    #         if nu is not None:  # nu provided
    #             ea = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))  # 11a: eccentric anomaly (from true anomaly) [rad]
    #             ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
    #             # print("Variant 3: ea-  nu+, calc ea ma")
    #         else:  # nu, ea not provided
    #             if ma is not None:  # ma provided
    #                 ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
    #                 nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
    #                 # print("Variant 4: ea-  nu-  ma+, calc ea nu")
    #             else:  # nu, ea, ma not provided
    #                 if T is None:  # T not provided
    #                     T = 0.0
    #                     print(f"{self.name}: L, ea, nu, ma, T missing, T set to default value 0.0")
    #                 n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
    #                 ma = n * T  # 1b: Mean anomaly at time of periapsis (from angular motion).
    #                 ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
    #                 nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
    #                 # print("Variant 5: ea-  nu-  ma-  T+, calc n ma ea nu")
    #
    #     n = math.sqrt(mu / a ** 3)  # 12a: mean angular motion
    #     T = ma / n  # Time of periapsis (from mean anomaly and angular motion). Just for completeness.
    #
    #     ma += t * n  # 1b
    #     ma %= 2 * math.pi
    #     ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
    #     nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
    #
    #     self.L, self.ma, self.ea, self.nu, self.T = L, ma, ea, nu, T  # save calculated parameters in body object
    #
    # def keplerian_elements_to_state_vector(self):
    #     """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
    #     Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis.
    #     [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
    #     [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
    #     [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
    #     [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
    #     [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py
    #     Numbers in comments refer to numbered formulas in [a] and [b].
    #     Code based on [c]. Added calculation of eccentric anomaly based on the explanations
    #     in [d] using a stripped down version of [e]."""
    #
    #     self.calc_orbit_angles()  # Omega, omega, pomega
    #     self.calc_period_or_semi_major_axis()  # P, a
    #     self.calc_anomalies()  # L, ma, ea, nu, T
    #     P, a, e, i, Omega, omega, pomega, L = self.P, self.a, self.e, self.i, self.Omega, self.omega, self.pomega, self.L  # for readability of formulas
    #     ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas
    #
    #     r = a * (1 - e * math.cos(ea))  # 4b: radius r
    #     h = math.sqrt(mu * a * (1 - e ** 2))  # 5b: specific angular momentum h
    #     x = r * (math.cos(Omega) * math.cos(omega + nu) - math.sin(Omega) * math.sin(omega + nu) * math.cos(i))  # 6b: position component x
    #     y = r * (math.sin(Omega) * math.cos(omega + nu) + math.cos(Omega) * math.sin(omega + nu) * math.cos(i))  # 6b: position component y
    #     z = r * (math.sin(i) * math.sin(omega + nu))  # 6b: position component z
    #     p = a * (1 - e ** 2)  # 7b: Semi-latus rectum. Used in velocity calculation.
    #     dx = (x * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.cos(Omega) * math.sin(omega + nu) + math.sin(Omega) * math.cos(omega + nu) * math.cos(i))  # 7b: velocity component x
    #     dy = (y * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.sin(Omega) * math.sin(omega + nu) - math.cos(Omega) * math.cos(omega + nu) * math.cos(i))  # 7b: velocity component y
    #     dz = (z * h * e / (r * p)) * math.sin(nu) + (h / r) * (math.cos(omega + nu) * math.sin(i))  # 7b: velocity component z
    #     return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors
    #
    # def calc_state_vector(self, p, bodies):
    #     """Get initial position and velocity of the physical body self."""
    #     self.mu = CurveSimPhysics.gravitational_parameter(bodies, p.g)  # is the same for all bodies in the system, because they are orbiting a common barycenter
    #     if self.velocity is None:  # State vectors are not in config file. So they will be calculated from Kepler orbit parameters instead.
    #         state_vector_function = self.keplerian_elements_to_state_vector
    #         pos, vel, *_ = state_vector_function()
    #         self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
    #         self.velocity = np.array(vel, dtype=float)  # [m/s] initial velocity
    #         self.velocity /= (1 + (self.mass / bodies[0].mass))  # correction because formulas seem to assume a system where all the mass is in one object at the center
    #
    # def state_vector_to_keplerian_elements(self):
    #     """Given the State Vector (position x, y, z and velocity dx, dy, dz) of an exoplanet, calculate its
    #         Kepler Orbit Elements (semi-major axis,  eccentricity, inclination, longitude of ascending node,
    #         argument of periapsis,  true anomaly) with a python function. You may assume that the orbit is
    #         well defined (no edge case, no hyperbole)"""
    #
    #     # Extract position and velocity components
    #     x, y, z = self.positions[0]
    #     dx, dy, dz = self.velocity
    #
    #     # Calculate specific angular momentum
    #     h_vec = np.cross([x, y, z], [dx, dy, dz])
    #     h = np.linalg.norm(h_vec)
    #
    #     # Calculate the semi-major axis
    #     r = np.linalg.norm([x, y, z])
    #     v = np.linalg.norm([dx, dy, dz])
    #     mu = self.mu
    #     a = 1 / (2 / r - v ** 2 / mu)
    #
    #     # Calculate the eccentricity vector and its magnitude
    #     e_vec = (np.cross([dx, dy, dz], h_vec) / mu) - np.array([x, y, z]) / r
    #     e = np.linalg.norm(e_vec)
    #
    #     # Calculate the inclination
    #     i = np.arccos(h_vec[2] / h)
    #
    #     # Calculate the longitude of ascending node
    #     n_vec = np.cross([0, 0, 1], h_vec)
    #     n = np.linalg.norm(n_vec)
    #     if n != 0:
    #         Omega = np.arccos(n_vec[0] / n)
    #         if n_vec[1] < 0:
    #             Omega = 2 * np.pi - Omega
    #     else:
    #         Omega = 0
    #
    #     # Calculate the argument of periapsis
    #     if n != 0:
    #         omega = np.arccos(np.dot(n_vec, e_vec) / (n * e))
    #         if e_vec[2] < 0:
    #             omega = 2 * np.pi - omega
    #     else:
    #         omega = 0
    #
    #     # Calculate the true anomaly
    #     nu = np.arccos(np.dot(e_vec, [x, y, z]) / (e * r))
    #     if np.dot([x, y, z], [dx, dy, dz]) < 0:
    #         nu = 2 * np.pi - nu
    #
    #     # Save calculated parameters in body object
    #     self.a = a
    #     self.e = e
    #     self.i = np.degrees(i)
    #     self.Omega = np.degrees(Omega)
    #     self.omega = np.degrees(omega)
    #     self.nu = np.degrees(nu)
    #
    #     return a, e, np.degrees(i), np.degrees(Omega), np.degrees(omega), np.degrees(nu)

    def full_eclipse(self, other, d):
        if self.radius < other.radius:  # Total eclipse
            area = self.area_2d
            relative_radius = 0
            return area, relative_radius
        else:  # Annular (i.e. ring) eclipse
            area = other.area_2d
            relative_radius = d / self.radius
            return area, relative_radius

    def partial_eclipse(self, other, d):
        # Eclipsed area is the sum of a circle segment of self + a circle segment of other
        # https://de.wikipedia.org/wiki/Kreissegment  https://de.wikipedia.org/wiki/Schnittpunkt#Schnittpunkte_zweier_Kreise
        self.d = (self.radius ** 2 - other.radius ** 2 + d ** 2) / (2 * d)  # Distance of center from self to radical axis
        other.d = (other.radius ** 2 - self.radius ** 2 + d ** 2) / (2 * d)  # Distance of center from other to radical axis
        other.h = other.radius + self.d - d  # Height of circle segment
        self.h = self.radius + other.d - d  # Height of circle segment
        other.angle = 2 * math.acos(1 - other.h / other.radius)  # Angle of circle segment
        self.angle = 2 * math.acos(1 - self.h / self.radius)  # Angle of circle segment
        other.eclipsed_area = other.radius ** 2 * (other.angle - math.sin(other.angle)) / 2  # Area of circle segment
        self.eclipsed_area = self.radius ** 2 * (self.angle - math.sin(self.angle)) / 2  # Area of circle segment
        area = other.eclipsed_area + self.eclipsed_area  # Eclipsed area is sum of two circle segments.
        relative_radius = (self.radius + self.d - other.h) / (2 * self.radius)  # Relative distance between approximated center C of eclipsed area and center of self
        return area, relative_radius

    # def find_tt_old(self, other, iteration, rebound_sim, p, sim_flux, time_s0, time_d, start_index, end_index, dt):
    #     """other eclipses self. Find the exact time of transit (TT).
    #         iteration should be the last one before TT. """
    #     eclipser = rebound_sim.particles[other.name]
    #     eclipsee = rebound_sim.particles[self.name]
    #     rebound_sim.integrate(time_s0[iteration])
    #     dx_left = eclipser.x - eclipsee.x
    #     t_left = rebound_sim.t
    #     rebound_sim.integrate(time_s0[iteration + 1])
    #     t_right = rebound_sim.t
    #     dx_right = eclipser.x - eclipsee.x
    #     interval_extensions = 0
    #     while dx_left * dx_right >= 0:  # dx per definition 0 at TT. If dx_left and dx_right have the same sign due to numeric instability in rebound, enlarge the search interval.
    #         t_left -= dt
    #         t_right += dt
    #         rebound_sim.integrate(t_left)
    #         dx_left = eclipser.x - eclipsee.x
    #         t_left = rebound_sim.t
    #         rebound_sim.integrate(t_right)
    #         t_right = rebound_sim.t
    #         dx_right = eclipser.x - eclipsee.x
    #         interval_extensions += 1
    #     if interval_extensions > 0 and p.verbose:
    #         print(f"{Fore.YELLOW}\nWARNING in function find_tt: Rebound integration results are possibly not accurate enough.")
    #         print(f"Try again with half the overall iteration time step parameter "dt".{Style.RESET_ALL}   ", end="")
    #         print(f"{iteration=}   {interval_extensions=}")
    #     if dx_left * dx_right < 0 and eclipser.z >= eclipsee.z:  # sign of dx changed and eclipser in front of eclipsee
    #         while t_right - t_left > 1e-1:  # bisect until desired precision reached
    #             t_middle = (t_right + t_left) / 2
    #             rebound_sim.integrate(t_middle)
    #             if dx_left * (eclipser.x - eclipsee.x) < 0:  # TT lies between t_left and t_middle
    #                 t_right = rebound_sim.t  # middle is now the new right
    #                 dx_right = eclipser.x - eclipsee.x
    #             else:  # TT lies between t_right and middle
    #                 t_left = rebound_sim.t  # middle is now the new left
    #                 dx_left = eclipser.x - eclipsee.x
    #         tt = rebound_sim.t / p.day + p.start_date
    #         d = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
    #         impact = d / self.radius
    #         close_enough = d <= self.radius + other.radius
    #         depth = 1 - sim_flux.interpolate_max_depth(tt, p, iteration, start_index, end_index, dt, time_d)
    #         return tt, impact, depth, close_enough
    #     else:
    #         print(f"{Fore.RED}\nERROR in function find_tt: Try with a smaller iteration time step dt.")
    #         print(f"If that does not help, please open an issue on https://github.com/lichtgestalter/curvesimulator/issues and provide your config file.{Style.RESET_ALL}")
    #         return -1, -1, -1, False

    def find_tt(self, other, iteration, rebound_sim, p, sim_flux, time_s0, time_d, start_index, end_index, dt):
        """other eclipses self. Find the exact time of transit (TT).
            iteration should be the last one before TT. """
        eclipser = rebound_sim.particles[other.name]
        eclipsee = rebound_sim.particles[self.name]
        rebound_sim.integrate(time_s0[iteration])
        dx_left = eclipser.x - eclipsee.x
        t_left = rebound_sim.t
        rebound_sim.integrate(time_s0[iteration + 1])
        t_right = rebound_sim.t
        dx_right = eclipser.x - eclipsee.x
        interval_extensions = 0
        while dx_left * dx_right >= 0:  # dx per definition 0 at TT. If dx_left and dx_right have the same sign due to numeric instability in rebound, enlarge the search interval.
            t_left -= dt
            t_right += dt
            rebound_sim.integrate(t_left)
            dx_left = eclipser.x - eclipsee.x
            t_left = rebound_sim.t
            rebound_sim.integrate(t_right)
            t_right = rebound_sim.t
            dx_right = eclipser.x - eclipsee.x
            interval_extensions += 1
            if interval_extensions > p.max_interval_extensions:
                if p.verbose:
                    print(f"{Fore.YELLOW}\nWARNING in function find_tt: Maximum acceptable interval extension exceeded.")
                    print(f"This is due to a too large iteration time step parameter <dt>{Style.RESET_ALL}   ", end="")
                    print(f"or due to an unstable star system.{Style.RESET_ALL}   ", end="")
                    print(f"Try again with half the iteration time step parameter <dt>{Style.RESET_ALL}   ", end="")
                    print(f"or choose more plausible start values and more restrictive upper/lower limits for the body parameters{Style.RESET_ALL}  ", end="")
                    print(f"Consider moving the time intervals a bit.{Style.RESET_ALL}   ", end="")
                    print(f"{iteration=}  {time_d[iteration]=} {interval_extensions=}")
                return -1, -1, -1, False
            if iteration - interval_extensions <= start_index or iteration + interval_extensions >= end_index:
                if p.verbose:
                    print(f"{Fore.YELLOW}\nWARNING in function find_tt: Possible TT at the edge of a time interval.")
                    print(f"Consider moving the time intervals a bit.{Style.RESET_ALL}   ", end="")
                    print(f"{iteration=}  {time_d[iteration]=}")
                return -1, -1, -1, False
        if dx_left * dx_right < 0 and eclipser.z >= eclipsee.z:  # sign of dx changed and eclipser in front of eclipsee
            while t_right - t_left > p.transit_precision:  # bisect until desired precision reached
                t_middle = (t_right + t_left) / 2
                rebound_sim.integrate(t_middle)
                if dx_left * (eclipser.x - eclipsee.x) < 0:  # TT lies between t_left and t_middle
                    t_right = rebound_sim.t  # middle is now the new right
                    dx_right = eclipser.x - eclipsee.x
                else:  # TT lies between t_right and middle
                    t_left = rebound_sim.t  # middle is now the new left
                    dx_left = eclipser.x - eclipsee.x
            tt = rebound_sim.t / p.day + p.start_date
            d = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
            impact = d / self.radius
            close_enough = d <= self.radius + other.radius
            if close_enough:
                depth = self.depth_at_tt(other, eclipser, eclipsee)
            else:
                depth = 0
            # print(f"{tt:12.6f};{eclipsee.vz:8.2f}")
            return tt, impact, depth, close_enough
        else:
            if p.verbose:
                print(f"{Fore.YELLOW}\nWARNING in function find_tt: Eclipser not in front of eclipsee at expected TT.")
                print(f"This is due to a too large iteration time step parameter <dt>{Style.RESET_ALL}   ", end="")
                print(f"or due to an unstable star system.{Style.RESET_ALL}   ", end="")
                print(f"Try again with half the iteration time step parameter <dt>{Style.RESET_ALL}   ", end="")
                print(f"or choose more plausible start values and more restrictive upper/lower limits for the body parameters{Style.RESET_ALL}  ", end="")
                print(f"Consider moving the time intervals a bit.{Style.RESET_ALL}   ", end="")
                print(f"{iteration=}  {time_d[iteration]=} {interval_extensions=}")
            return -1, -1, -1, False

    # def find_t1234_old(self, other, iteration, rebound_sim, time_s0, start_index, end_index, p, transittimetype):
    #     """other eclipses self. Find where ingress starts (T1) or egress ends (T4)."""
    #     eclipser = rebound_sim.particles[other.name]
    #     eclipsee = rebound_sim.particles[self.name]
    #     if transittimetype in ["T1", "T4"]:
    #         d_max = self.radius + other.radius
    #     else:
    #         d_max = abs(self.radius - other.radius)
    #     iteration_delta = 0
    #     d = -1
    #     step = -1 if transittimetype in ["T1", "T2"] else 1
    #     # T1/T2: go backwards from iteration (this should be the one right _after_ TT) to find the iteration before the eclipse starts
    #     # T3/T4: go forward from iteration (this should be the one right _before_ TT) to find the iteration after the eclipse ends
    #     while d < d_max:
    #         if iteration + iteration_delta >= end_index or iteration + iteration_delta < start_index:
    #             return None  # incomplete transit at start or end of current simulation interval
    #         iteration_delta += step
    #         d = CurveSimPhysics.distance_2d(other, self, iteration + iteration_delta)
    #     rebound_sim.integrate((time_s0[iteration + iteration_delta]))
    #     d_old = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
    #     t_old = rebound_sim.t
    #     rebound_sim.integrate(time_s0[iteration])
    #     t_new = rebound_sim.t
    #     d_new = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
    #     if transittimetype not in ["T1", "T2"]:
    #         t_new, t_old = t_old, t_new
    #     if d_old > d_max > d_new:  # T1 or T2  or T3 or T4 lies between t_old and t_new
    #         while t_new - t_old > 1e-1:  # bisect until desired precision reached
    #             rebound_sim.integrate((t_new + t_old) / 2)
    #             in_eclipse = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee) < d_max
    #             if transittimetype in ["T1", "T2"]:
    #                 if in_eclipse: # T1 or T2 lies between t_old and (t_new + t_old) / 2
    #                     t_new = rebound_sim.t
    #                 else:
    #                     t_old = rebound_sim.t
    #             else:
    #                 if in_eclipse: # T3 or T4 lies between t_new and (t_new + t_old) / 2
    #                     t_old = rebound_sim.t
    #                 else:
    #                     t_new = rebound_sim.t
    #         return rebound_sim.t / p.day + p.start_date
    #     else:  # grazing transit (or rebound inaccuracy?)
    #         return None

    def find_t1234(self, other, tt_s0, iteration, rebound_sim, time_s0, start_index, end_index, p, transittimetype):
        """other eclipses self. Find where ingress starts (T1) or egress ends (T4)."""
        eclipser = rebound_sim.particles[other.name]
        eclipsee = rebound_sim.particles[self.name]
        if transittimetype in ["T1", "T4"]:
            d_event = self.radius + other.radius  # distance at T1, T4
        else:
            d_event = abs(self.radius - other.radius)  # distance at T2, T3
        iteration_delta = 0
        d = -1
        step = -1 if transittimetype in ["T1", "T2"] else 1
        # T1/T2: go backwards from iteration (this should be the one right _after_ TT) to find the iteration before the eclipse starts
        # T3/T4: go forward from iteration (this should be the one right _before_ TT) to find the iteration after the eclipse ends
        while d < d_event:
            if iteration + iteration_delta >= end_index or iteration + iteration_delta < start_index:
                return None  # incomplete transit at start or end of current simulation interval
            iteration_delta += step
            d = CurveSimPhysics.distance_2d_body(other, self, iteration + iteration_delta)
        rebound_sim.integrate((time_s0[iteration + iteration_delta]))
        d_old = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
        t_old = rebound_sim.t
        rebound_sim.integrate(tt_s0)
        t_new = rebound_sim.t
        d_new = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
        if transittimetype not in ["T1", "T2"]:
            t_new, t_old = t_old, t_new  # T1 or T2  or T3 or T4 lies between t_old and t_new
        while t_new - t_old > p.transit_precision:  # bisect until desired precision reached
            rebound_sim.integrate((t_new + t_old) / 2)
            d = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
            in_eclipse = d < d_event
            if transittimetype in ["T1", "T2"]:
                if in_eclipse:  # T1 or T2 lies between t_old and (t_new + t_old) / 2
                    t_new = rebound_sim.t
                else:
                    t_old = rebound_sim.t
            else:
                if in_eclipse:  # T3 or T4 lies between t_new and (t_new + t_old) / 2
                    t_old = rebound_sim.t
                else:
                    t_new = rebound_sim.t
        if abs(rebound_sim.t - tt_s0) < p.transit_precision:
            return None
        else:
            return rebound_sim.t / p.day + p.start_date

    # def eclipsed_by(self, other, iteration, p):
    #     """Returns area, relative_radius
    #     area: Area of self which is eclipsed by other.
    #     relative_radius: The distance of the approximated center of the eclipsed area from the center of self as a percentage of self.radius (used for limb darkening)."""
    #     # if other.positions[iteration][0] < self.positions[iteration][0]:  # Is other nearer to viewpoint than self? (i.e. its position has a smaller x-coordinate)
    #     if other.positions[iteration][2] > self.positions[iteration][2]:  # Is other nearer to viewpoint than self? (i.e. its position has a larger z-coordinate)
    #         d = CurveSimPhysics.distance_2d_body(other, self, iteration)
    #         if d < self.radius + other.radius:  # Does other eclipse self?
    #             if d <= abs(self.radius - other.radius):  # Annular (i.e. ring) eclipse or total eclipse
    #                 area, relative_radius = self.full_eclipse(other, d)
    #             else:  # Partial eclipse
    #                 area, relative_radius = self.partial_eclipse(other, d)
    #             return area, relative_radius
    #         else:  # No eclipse because, seen from viewer, the bodies are not close enough to each other
    #             return None, None
    #     else:  # other cannot eclipse self, because self is nearer to viewer than other
    #         return None, None

    def eclipsed_by(self, other, iteration, p):
        """Returns area, relative_radius
        area: Area of self which is eclipsed by other.
        relative_radius: The distance of the approximated center of the eclipsed area from the center of self as a percentage of self.radius (used for limb darkening)."""
        # if other.positions[iteration][0] < self.positions[iteration][0]:  # Is other nearer to viewpoint than self? (i.e. its position has a smaller x-coordinate)
        if other.positions[iteration][2] < self.positions[iteration][2]:  # Is other nearer to viewpoint than self? (i.e. its position has a larger z-coordinate)
            # print(f"{iteration=:5d}: {other.name} cannot eclipse {self.name}, because it is further away.")
            return None, None  # other cannot eclipse self, because self is nearer to viewer than other
        if abs(other.positions[iteration][0] - self.positions[iteration][0]) > self.radius + other.radius:
            # print(f"{iteration=:5d}: {other.name} cannot eclipse {self.name}, because their distance in x-direction is too large.")
            return None, None  # difference in x-coordinate too large
        if abs(other.positions[iteration][1] - self.positions[iteration][1]) > self.radius + other.radius:
            # print(f"{iteration=:5d}: {other.name} cannot eclipse {self.name}, because their distance in y-direction is too large.")
            return None, None  # difference in y-coordinate too large
        # print(f"{iteration=:5d}: YAY!!!")
        d = CurveSimPhysics.distance_2d_body(other, self, iteration)
        if d < self.radius + other.radius:  # Does other eclipse self?
            if d <= abs(self.radius - other.radius):  # Annular (i.e. ring) eclipse or total eclipse
                area, relative_radius = self.full_eclipse(other, d)
            else:  # Partial eclipse
                area, relative_radius = self.partial_eclipse(other, d)
            return area, relative_radius
        else:  # No eclipse because, seen from viewer, the bodies are not close enough to each other
            return None, None

    def eclipsed_by_at_tt(self, other, eclipser, eclipsee):
        """ self, other: body
        eclipser, eclipsee: Rebound Particle
        eclipser is the Rebound Particle of self
        eclipsee is the Rebound Particle of other
        Returns area, relative_radius
        area: Area of self which is eclipsed by eclipser.
        relative_radius: The distance of the approximated center
        of the eclipsed area from the center of eclipsee as
        a percentage of eclipsee.radius (used for limb darkening)."""
        if eclipser.z < eclipsee.z:  # Is eclipser nearer to viewpoint than eclipsee? (i.e. its position has a larger z-coordinate)
            return None, None  # eclipser cannot eclipse eclipsee, because eclipsee is nearer to viewer than eclipser
        d = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
        if d < self.radius + other.radius:  # Does other eclipse self?
            if d <= abs(self.radius - other.radius):  # Annular (i.e. ring) eclipse or total eclipse
                area, relative_radius = self.full_eclipse(other, d)
            else:  # Partial eclipse
                area, relative_radius = self.partial_eclipse(other, d)
            return area, relative_radius
        else:  # No eclipse because, seen from viewer, the bodies are not close enough to each eclipser
            return None, None

    def depth_at_tt(self, other, eclipser, eclipsee):
        eclipsed_area, relative_radius = self.eclipsed_by_at_tt(other, eclipser, eclipsee)
        if eclipsed_area is not None:
            limbdarkening = CurveSimPhysics.limbdarkening(relative_radius, self.limb_darkening_u1, self.limb_darkening_u2)
            relative_depth = (self.intensity * eclipsed_area * limbdarkening / self.mean_intensity) / self.luminosity
            return relative_depth
        return None

    def calc_frames_per_orbit(self, p):
        """Calculates for each body how many video frames are needed to complete one orbit.
           ffmpeg (or the video display program?) tends to omit the last few frames.
           Therefore add a handful of extra frames."""
        if self.P is not None:
            return self.P / (p.dt * p.sampling_rate)
        else:
            return None

    def print_particle(self, simulation):
        particle = simulation.particles[self.name]
        print(f"\n{self.name}:")
        print(f"x:     {particle.x:14.6e}")
        print(f"y:     {particle.y:14.6e}")
        print(f"z:     {particle.z:14.6e}")
        print(f"vx:    {particle.vx:14.6e}")
        print(f"vy:    {particle.vy:14.6e}")
        print(f"vz:    {particle.vz:14.6e}")
        try:
            orbit = particle.orbit()
        except ValueError:
            orbit = False
        if orbit:
            print(f"a:     {particle.a:14.6e}")
            print(f"P:     {particle.P:14.6e}")
            print(f"e:     {particle.e:14.6e}")
            print(f"i:     {particle.inc:14.6e}")
            print(f"Omega: {particle.Omega:14.6e}")
            print(f"omega: {particle.omega:14.6e}")
            print(f"ma:    {particle.M:14.6e}")


class CurveSimLightcurve(np.ndarray):

    def __new__(cls, shape, dtype=float):
        obj = np.zeros(shape, dtype=dtype).view(cls)
        return obj

    def __str__(self):
        return f"CurveSimLightcurve: max={max(self)*100:.4f}%, min={min(self)*100:.4f}%, len={len(self)}"

    def interpolate_max_depth(self, tt, p, iteration, start_index, end_index, dt, time_d):
        """
        Interpolates the "self" value at a given "tt" using cubic interpolation
        (Catmull-Rom like) based on surrounding "iteration" points.

        Args:
            self: sim_flux (Lightcurve)
            tt: The time value for which to interpolate [BJD]
            p: CurveSimulator parameters
            iteration: index of the simulation right before TT (iteration < iteration_tt < iteration + 1).
            start_index: index of the first iteration in the current interval (parameters "starts" and "ends")
            end_index: index + 1 of the last iteration in the current interval (parameters "starts" and "ends")

        Returns:
            The interpolated value at tt, or 1 if interpolation indices are out of bounds. (flux = 1 means depth = 0)
        """
        # if not (1 <= iteration < len(self) - 2):  # Ensure indices are within bounds
        if not (start_index < iteration < end_index - 2):  # Ensure indices are within bounds
            if p.verbose:
                print(f"{Fore.YELLOW}\nWARNING: Function interpolate_max_depth: Interpolation indices out of bounds at {iteration=}")
                print(f"Depth of this transit has been stored in result file as 0.")
                print(f"Try to move the intervals (parameters <starts> and <ends>) a bit.{Style.RESET_ALL}")
            return 1
        # iteration_tt = (tt - (time_s0[iteration] / p.day + p.start_date)) / (dt /p.day) + iteration
        iteration_tt = (tt - time_d[iteration]) / (dt /p.day) + iteration
        P0 = self[iteration - 1]  # f_im1
        P1 = self[iteration]  # f_i
        P2 = self[iteration + 1]  # f_ip1
        P3 = self[iteration + 2]  # f_ip2
        alpha = iteration_tt - iteration  # Calculate the normalized position (alpha or t) within the segment [iteration, iteration + 1]
        alpha = np.clip(alpha, 0.0, 1.0)  # Due to floating point arithmetic, it might be slightly outside [0, 1), so clamp it.
        alpha2 = alpha * alpha
        alpha3 = alpha2 * alpha
        interpolated_value = 0.5 * (
                (2 * P1) +
                (-P0 + P2) * alpha +
                (2 * P0 - 5 * P1 + 4 * P2 - P3) * alpha2 +
                (-P0 + 3 * P1 - 3 * P2 + P3) * alpha3
        )
        return interpolated_value

    def save_sim_flux(self, p, time_d):
        noisy_flux = self + np.random.normal(0, p.sim_flux_err, self.shape)
        flux_err = np.full(self.shape, p.sim_flux_err)
        data = np.column_stack((time_d, noisy_flux, flux_err))
        np.savetxt(p.sim_flux_file, data, delimiter=",", header="time,flux,flux_err", comments="")
        if p.verbose:
            print(f"Saved simulated flux to {p.sim_flux_file} including white noise with standard deviation {p.sim_flux_err}")



class FittingGUI:
    """
    A Tkinter application featuring self.number_of_params interactive data control parameters
    and a real-time Matplotlib plot.
    Activation of an parameter is handled via a mouse click on the parameter.
    """
    def __init__(self, root, p, bodies, time_s0, time_d, measured_tt):
        self.bodies = bodies
        self.time_s0 = time_s0
        self.time_d = time_d
        self.measured_tt = measured_tt

        self.residuals_tt_sum_squared, self.measured_tt = self.run_simulation(p)

        self.root = root
        self.root.title("CurveSimulator Manual Fitter")
        self.root.state("zoomed")  # Maximize the window

        # Internal state for the self.number_of_params parameters: value, delta, and Tkinter variables
        self.parameters = []
        self.parameter_frames = []
        self.active_parameter_index = 0  # First parameter is active by default
        self.param_rows = 2
        self.params_per_row = (len(p.fitting_parameters) + 1) // 2
        self.number_of_params = self.param_rows * self.params_per_row
        self.number_of_params = len(p.fitting_parameters)

        self._initialize_grid_weights()
        self._initialize_parameter_data(p)
        self._create_widgets(p)
        self._setup_bindings(p)

        self.update_plot(self.measured_tt)  # Initial plot
        self.activate_parameter(self.active_parameter_index)  # Initial highlighting

    def _initialize_grid_weights(self):
        """Set grid weights for main window to achieve proportional sizing."""
        # Main layout: 3 rows (8%, 8%, 84%) and self.params_per_row columns
        total_weight = 100
        parameter_row_weight = math.ceil(total_weight * 0.08)  # Calculate a suitable weight for the parameter rows (8% of total height)
        self.root.grid_rowconfigure(0, weight=parameter_row_weight)  # Row 0 (Parameters 1-self.params_per_row): 8%
        self.root.grid_rowconfigure(1, weight=parameter_row_weight)  # Row 1 (Parameters 7-self.number_of_params): 8%
        self.root.grid_rowconfigure(2, weight=total_weight - 2 * parameter_row_weight)  # Row 2 (Plot): 84%
        for i in range(self.params_per_row):
            self.root.grid_columnconfigure(i, weight=1)  # self.params_per_row columns for the self.params_per_row parameters per row

    def _initialize_parameter_data(self, p):
        """Initialize internal data model and Tkinter variables for all self.number_of_params parameters."""
        for fp in p.fitting_parameters:
            value = fp.startvalue * fp.scale
            delta = 10 ** round((math.log10(abs(value / 100))), 0) if value != 0 else 1
            fmt_value = FittingGUI._get_format(value)
            fmt_delta = FittingGUI._get_format(delta)
            parameter_data = {
                "value": value,
                "delta": delta,
                "value_var": StringVar(value=f"{value:{fmt_value}}"),
                "delta_var": StringVar(value=f"{delta:{fmt_delta}}")
            }
            self.parameters.append(parameter_data)

    def _create_widgets(self, p):
        """
        Create and place all self.number_of_params control parameters and the Matplotlib plot.
        The activation is triggered by a mouse click on the parameter frame or its contents.
        """
        # 1. Control Parameters (2 rows x self.params_per_row columns)
        for i, fp in enumerate(p.fitting_parameters):
            parameter = self.parameters[i]
            frame = ttk.Frame(self.root, padding="5 5 5 5", relief=tk.RIDGE)  # Create a Frame for the parameter (for grouping and highlighting)
            row = i // self.params_per_row  # Calculate grid position: Row 0/1, Column 0-5
            col = i % self.params_per_row

            frame.grid(row=row, column=col, sticky="nsew", padx=5, pady=5)
            self.parameter_frames.append(frame)
            frame.bind("<Button-1>", lambda event, idx=i: self.activate_parameter(idx))  # Bind mouse click event to the frame to activate the parameter

            for j in range(3):  # Configure frame's internal grid (3 rows: Title, Value, Delta)
                frame.grid_rowconfigure(j, weight=1)
            frame.grid_columnconfigure(0, weight=1)
            frame.grid_columnconfigure(1, weight=2)
            # A. Title Label (Row 0)
            title_label = ttk.Label(frame, text=fp.long_body_parameter_name, font=("Inter", 10, "bold"))
            title_label.grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 5))
            # B. Value Entry (Row 1)
            value_label = ttk.Label(frame, text="Value:", font=("Inter", 9))
            value_label.grid(row=1, column=0, sticky="w")
            value_entry = ttk.Entry(frame, textvariable=parameter["value_var"], width=10, state="readonly")
            value_entry.grid(row=1, column=1, sticky="ew", padx=(5, 0))
            # C. Delta Entry (Row 2)
            delta_label = ttk.Label(frame, text="Delta:", font=("Inter", 9))
            delta_label.grid(row=2, column=0, sticky="w")
            delta_entry = ttk.Entry(frame, textvariable=parameter["delta_var"], width=10, state="readonly")
            delta_entry.grid(row=2, column=1, sticky="ew", padx=(5, 0))
            for widget in (title_label, value_label, value_entry, delta_label, delta_entry):  # Bind child widgets to the same activation command so clicking on them works too
                widget.bind("<Button-1>", lambda event, idx=i: self.activate_parameter(idx))

        # Matplotlib Plot (Row 2)
        plot_frame = ttk.Frame(self.root, padding="10 10 10 10", relief=tk.SUNKEN)
        plot_frame.grid(row=2, column=0, columnspan=self.params_per_row, sticky="nsew", padx=5, pady=5)
        plot_frame.grid_rowconfigure(0, weight=1)
        plot_frame.grid_columnconfigure(0, weight=1)

        # # Create Matplotlib Figure
        self.unique_eclipsers = self.measured_tt["eclipser"].unique()
        n_eclipsers = len(self.unique_eclipsers)
        self.fig, self.axes = plt.subplots(n_eclipsers, figsize=(10, 3.5 * n_eclipsers), sharex=True)
        if n_eclipsers == 1:
            self.axes = [self.axes]
        abs_min = abs(min(self.measured_tt["delta"].min(), 0))
        abs_max = abs(max(self.measured_tt["delta"].max(), 0))
        ylim = (-max(abs_min, abs_max), max(abs_min, abs_max))

        for ax, eclipser in zip(self.axes, self.unique_eclipsers):
            df = self.measured_tt[self.measured_tt["eclipser"] == eclipser]
            ax.plot(df["tt"], df["delta"], marker="o", linestyle="-", color="blue", alpha=0.7)
            ax.axhline(0, color="gray", linestyle="dashed", linewidth=1)
            ax.set_ylabel(f"TT Delta [days]")
            ax.set_title(f"Eclipser: {eclipser}")
            ax.tick_params(labelbottom=True)
            ax.set_ylim(ylim)
        self.axes[-1].set_xlabel("Transit Time [BJD]")
        self.fig.suptitle(f"TT Delta...", fontsize=14)
        plt.tight_layout(rect=(0, 0, 1, 0.97))

        # self.fig = Figure(figsize=(5, 4), dpi=100)
        # self.ax = self.fig.add_subplot(111)
        # self.ax.set_title("Parameter 1 (X) vs. Parameter 2 (Y) Plot")
        # self.ax.set_xlabel("Parameter 1 Value")
        # self.ax.set_ylabel("Parameter 2 Value")
        # self.scatter, = self.ax.plot([], [], "o", color="royalblue", markersize=10)
        # self.ax.grid(True, linestyle="--", alpha=0.6)

        # Embed the figure into Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.grid(row=0, column=0, sticky="nsew")

    def _setup_bindings(self, p):
        """Bind some keys to the root window, passing "p" to the handler."""
        for key in ["<Up>", "<Down>", "<Left>", "<Right>"]:
        # for key in ["<Up>", "<Down>", "<Left>", "<Right>", "<Escape>"]:
            self.root.bind(key, lambda event, p=p: self.handle_key(p, event))

    @staticmethod
    def _get_format(number):
        if 1e-6 < abs(number) < 1e6:
            return ".10f"
        elif number == 0:
            return ".0f"
        else:
            return ".6e"

    def activate_parameter(self, index):
        """Sets the currently active parameter, visually highlights it, and ensures key focus."""
        if self.active_parameter_index is not None and self.active_parameter_index < len(self.parameter_frames):
            self.parameter_frames[self.active_parameter_index]["style"] = "TFrame"  # Reset color of the previously active frame
        self.active_parameter_index = index
        self.root.style = ttk.Style()  # Define a style for the active frame and apply it
        self.root.style.configure("Active.TFrame", background="#e0ffe0", relief="solid", borderwidth=2)  # Use a solid border with a light green background for clear highlighting
        self.parameter_frames[self.active_parameter_index]["style"] = "Active.TFrame"
        self.root.focus_set()  # Crucially, set focus back to the root window so arrow key bindings work

    def handle_key(self, p, event):
        """Processes arrow key presses to modify the active parameter's value or delta."""
        if self.active_parameter_index is None:
            return
        parameter = self.parameters[self.active_parameter_index]
        key = event.keysym
        if key == "Up":
            parameter["delta"] *= 10
        elif key == "Down":
            parameter["delta"] /= 10
        elif key == "Left":
            parameter["value"] -= parameter["delta"]
        elif key == "Right":
            parameter["value"] += parameter["delta"]
        # elif key == "Escape":
        #     parameter["value"] = 1.0
        #     parameter["delta"] = 1.0
        self.update_entry_fields(p)
        self.update_fitting_parameters(p)
        self.residuals_tt_sum_squared, self.measured_tt = self.run_simulation(p)
        max_delta = max(np.abs(self.measured_tt["delta"]))
        mean_delta = np.mean(np.abs(self.measured_tt["delta"]))
        print(f"\n{max_delta=:2.4f}   {mean_delta=:2.4f}    [days] ")
        self.update_plot(self.measured_tt)

    def run_simulation(self, p):
        param_references = [(fp.body_index, fp.parameter_name) for fp in p.fitting_parameters]
        for (body_index, parameter_name), fp in zip(param_references, p.fitting_parameters):
            self.bodies[body_index].__dict__[parameter_name] = fp.startvalue
        sim_rv, sim_flux, rebound_sim = self.bodies.calc_physics(p, self.time_s0)  # run simulation
        return CurveSimMCMC.match_transit_times(self.measured_tt, p, rebound_sim, sim_flux, self.time_d, self.time_s0)

    def update_entry_fields(self, p):
        """Update the Tkinter entry variables from the internal data model."""
        for parameter, fp in zip(self.parameters, p.fitting_parameters):
            fmt_value = FittingGUI._get_format(parameter["value"])
            fmt_delta = FittingGUI._get_format(parameter["delta"])
            parameter["value_var"].set(f"{parameter["value"]:{fmt_value}}")
            parameter["delta_var"].set(f"{parameter["delta"]:{fmt_delta}}")

    def update_fitting_parameters(self, p):
        for parameter, fp in zip(self.parameters, p.fitting_parameters):
            fp.startvalue = parameter["value"] / fp.scale

    # def update_plot(self, measured_tt):
    #     abs_min = abs(min(measured_tt["delta"].min(), 0))
    #     abs_max = abs(max(measured_tt["delta"].max(), 0))
    #     ylim = (-max(abs_min, abs_max), max(abs_min, abs_max))
    #
    #     for ax, eclipser in zip(self.axes, self.unique_eclipsers):
    #         df = measured_tt[measured_tt["eclipser"] == eclipser]
    #         ax.plot(df["tt"], df["delta"], marker="o", linestyle="-", color="blue", alpha=0.7)
    #         ax.set_ylim(ylim)
    #
    #     self.canvas.draw_idle()  # Redraw the canvas

    def update_plot(self, measured_tt):
        colors = ["blue", "black", "gray", "lightgrey"]
        if not hasattr(self, "plot_lines"):
            self.plot_lines = [[] for _ in self.axes]

        abs_min = abs(min(measured_tt["delta"].min(), 0))
        abs_max = abs(max(measured_tt["delta"].max(), 0))
        ylim = (-max(abs_min, abs_max), max(abs_min, abs_max))

        for ax_idx, (ax, eclipser) in enumerate(zip(self.axes, self.unique_eclipsers)):
            df = measured_tt[measured_tt["eclipser"] == eclipser]
            line, = ax.plot(df["tt"], df["delta"], marker="o", linestyle="-", color=colors[0], alpha=0.7)
            self.plot_lines[ax_idx].insert(0, line)
            if len(self.plot_lines[ax_idx]) > 4:
                old_line = self.plot_lines[ax_idx].pop()
                old_line.remove()
            for i, l in enumerate(self.plot_lines[ax_idx]):
                l.set_color(colors[i])
            ax.set_ylim(ylim)

        self.canvas.draw_idle()



    # def update_plot(self):
    #     """Update the Matplotlib plot with the new coordinates (Parameter1.value, Parameter2.value)."""
    #     x_val = self.parameters[0]["value"]  # Parameter 1 value
    #     y_val = self.parameters[1]["value"]  # Parameter 2 value
    #
    #     self.scatter.set_data([x_val], [y_val])  # Update scatter plot data
    #
    #     # Adjust plot limits dynamically (10% padding based on current values)
    #     # Consider a base range to prevent initial plot from being too small
    #     base_range = 10
    #     all_values = [a["value"] for a in self.parameters]
    #
    #     # Calculate overall min/max value across all self.number_of_params parameters
    #     min_v = min(all_values)
    #     max_v = max(all_values)
    #
    #     # Determine plot limits based on the overall range plus padding
    #     padding = (max_v - min_v) * 0.1
    #     if padding < base_range:
    #         padding = base_range
    #
    #     x_min = min_v - padding
    #     x_max = max_v + padding
    #     y_min = min_v - padding
    #     y_max = max_v + padding
    #
    #     self.ax.set_xlim(x_min, x_max)
    #     self.ax.set_ylim(y_min, y_max)
    #
    #     # Redraw the canvas
    #     self.canvas.draw_idle()
    #
class CurveSimManualFit:
    def __init__(self, p, bodies, time_s0, time_d, measured_tt):
        root = tk.Tk()
        app = FittingGUI(root, p, bodies, time_s0, time_d, measured_tt)
        root.focus_set()  # Give focus to the main window so key presses are immediately captured
        root.mainloop()

    def save_lmfit_results(self, p):
        pass

# from curvesimulator.cs_flux_data import csv2df
# from curvesimulator.curvesim import CurveSimulator  circular import!


def stopwatch():
    def decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            tic = time.perf_counter()
            result = func(self, *args, **kwargs)
            toc = time.perf_counter()
            print(f" {func.__name__}: {toc - tic:6.1f} seconds")
            return result

        return wrapper

    return decorator

def append_line_locked(filename, line, wait=0.1):
    """
    Append a line to `filename` ensuring all processes succeed:
    - try a non-blocking lock
    - on contention sleep `wait` seconds and retry
    - write, flush, fsync, release lock
    Works on Windows (msvcrt) and Unix (fcntl).
    """
    dirpath = os.path.dirname(filename)
    if dirpath:
        os.makedirs(dirpath, exist_ok=True)

    while True:
        if os.name == "nt":
            try:
                with open(filename, "a", encoding="utf8") as f:
                    fd = f.fileno()
                    try:
                        msvcrt.locking(fd, msvcrt.LK_NBLCK, 1)  # try non-blocking lock of 1 byte
                    except OSError:
                        time.sleep(wait)
                        continue
                    try:
                        f.write(line + "\n")
                        f.flush()
                        os.fsync(fd)
                    finally:
                        try:
                            msvcrt.locking(fd, msvcrt.LK_UNLCK, 1)
                        except Exception:
                            pass
                    return
            except Exception:
                raise  # propagate unexpected errors (e.g., permission issues)
        else:
            try:
                with open(filename, "a", encoding="utf8") as f:
                    try:
                        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)  # try non-blocking exclusive lock
                    except (BlockingIOError, OSError):
                        time.sleep(wait)
                        continue
                    try:
                        f.write(line + "\n")
                        f.flush()
                        os.fsync(f.fileno())
                    finally:
                        try:
                            fcntl.flock(f, fcntl.LOCK_UN)
                        except Exception:
                            pass
                    return
            except Exception:
                raise

class CurveSimMCMC:

    def __init__(self, p, bodies, time_s0, time_d, measured_flux_array, flux_err, measured_tt, dummy_object=False):
        self.results_directory = p.results_directory
        if dummy_object:
            return
        os.environ["OMP_NUM_THREADS"] = "1"  # Some builds of NumPy automatically parallelize some operations. This can cause problems when multi processing inside emcee is enabled. Turn that off by setting the environment variable OMP_NUM_THREADS=1.
        if not (p.flux_file or p.tt_file or p.rv_file):
            print(f"{Fore.RED}\nERROR: No measurements for fitting have been provided.{Style.RESET_ALL}")
            sys.exit(1)
        if os.path.exists("residual.tmp"):
            os.remove("residual.tmp")
        if os.path.exists("iteration.tmp"):
            os.remove("iteration.tmp")
        self.fitting_parameters = p.fitting_parameters
        self.moves = p.moves
        self.walkers = p.walkers
        self.thin_samples = p.thin_samples
        self.thin_samples_plot = max(p.thin_samples, 10)  # avoid unnecessary memory usage for some plots
        self.burn_in = p.burn_in
        self.chunk_size = p.chunk_size
        self.bins = p.bins
        self.unit = p.unit
        self.scale = p.scale
        self.steps = p.steps
        self.credible_mass = 0.68
        self.param_references = [(fp.body_index, fp.parameter_name) for fp in self.fitting_parameters]  # list of names of fitting parameters. Needed so these parameters can be updated inside log_likelihood().
        self.body_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in self.fitting_parameters]
        self.long_body_parameter_names = [fpn + " [" + self.unit[fpn.split(".")[-1]] + "]" for fpn in self.body_parameter_names]
        for fp, fpn, fpnu in zip(p.fitting_parameters, self.body_parameter_names, self.long_body_parameter_names):
            fp.body_parameter_name = fpn
            fp.long_body_parameter_name = fpnu
        self.param_bounds = [(fp.lower, fp.upper) for fp in self.fitting_parameters]
        self.ndim = len(self.param_references)
        self.theta0 = self.random_initial_values()
        self.args = (self.param_bounds, self.param_references, bodies, time_s0, time_d, measured_flux_array, flux_err, measured_tt, p)
        self.moves = eval(self.moves)
        self.acceptance_fractions = []
        self.integrated_autocorrelation_time = []
        self.start_real_time = time.strftime("%d.%m.%y %H:%M:%S")
        self.start_timestamp = time.perf_counter()
        self.max_likelihood_avg_residual_in_std = []
        self.mean_avg_residual_in_std = []
        self.median_avg_residual_in_std = []
        self.trace_plot_ok = True
        self.corner_plot_ok = True
        self.autocorrelation_function_plot_ok = True
        self.acceptance_plot_ok = True


        if p.backend:
            self.backend = emcee.backends.HDFBackend(p.backend)
            if p.load_backend:
                print(f"Loading backend from {p.backend}. Contains {self.backend.iteration} iterations.")
                steps_done = self.backend.iteration - self.burn_in
                self.loaded_steps = steps_done
                if steps_done < 0:
                    print(f"{Fore.RED}\nERROR: Backend contains less iterations than burn-in. Uncomment load_backend in the config file to start from scratch.{Style.RESET_ALL}")
                    sys.exit(1)
            else:
                print("Ignoring and resetting backend.")
                self.backend.reset(p.walkers, self.ndim)  # clear/reset the backend in case the file already exists
                steps_done, self.loaded_steps = 0, 0
        else:
            print("Running without backend.")
            self.backend = None
            steps_done, self.loaded_steps = 0, 0

        if p.mcmc_multi_processing:
            with Pool() as pool:  # enable multi processing
                self.mcmc_fit(p, bodies, time_s0, time_d, measured_flux_array, flux_err, measured_tt, steps_done, pool)
        else:
            self.mcmc_fit(p, bodies, time_s0, time_d, measured_flux_array, flux_err, measured_tt, steps_done)

    def mcmc_fit(self, p, bodies, time_s0, time_d, measured_flux_array, flux_err, measured_tt, steps_done, pool=None):
        self.sampler = emcee.EnsembleSampler(p.walkers, self.ndim, CurveSimMCMC.log_probability, pool=pool, moves=self.moves, args=self.args, backend=self.backend)
        # self.sampler = emcee.EnsembleSampler(p.walkers, self.ndim, CurveSimMCMC.log_probability, pool=pool, moves=self.moves, args=self.args)
        if not p.load_backend:
            self.theta = self.sampler.run_mcmc(self.theta0, self.burn_in, progress=True)
        else:
            self.theta = self.theta0.copy()
        for chunk in range(1, self.steps // self.chunk_size):
            self.theta = self.sampler.run_mcmc(self.theta, self.chunk_size, progress=True)
            steps_done += self.chunk_size
            self.mcmc_results(p, bodies, steps_done, time_s0, time_d, measured_tt, measured_flux_array, flux_err, chunk)

    def __repr__(self):
        return f"CurveSimMCMC with {self.walkers} walkers."

    @staticmethod
    def single_run(p, bodies=None, time_s0=None, time_d=None):
        if time_s0 is None and time_d is None:
            time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
        sim_rv, sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # Calculate all body positions and the resulting lightcurve
        results = bodies.find_transits(rebound_sim, p, sim_flux, time_s0, time_d)
        results["chi_squared_tt"], results["chi_squared_rv"], results["chi_squared_flux"], results["chi_squared_total"] = 0, 0, 0, 0
        results["measurements_tt"], results["measurements_rv"], results["measurements_flux"], results["measurements_total"] = 0, 0, 0, 0
        if p.video_file:
            CurveSimAnimation(p, bodies, sim_rv, sim_flux, time_s0)  # Create a video
        if p.tt_file:
            measured_tt = CurveSimResults.get_measured_tt(p)
            residuals_tt_sum_squared, measured_tt = CurveSimMCMC.match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0)
            measured_tt = results.calc_tt_chi_squared(measured_tt, p.free_parameters)  # store chi squared and p-value in results
            CurveSimMCMC.tt_delta_plot(p, 0, "tt_o_vs_c.png", measured_tt)  # compare observed vs. computed TT
        else:
            measured_tt = None
        if p.rv_file:
            measured_rv = CurveSimResults.get_measured_rv(p)
            measured_rv = CurveSimResults.calc_rv_residuals(measured_rv, p.rv_body, rebound_sim)  # compare observed vs. computed RV
            measured_rv = results.calc_rv_chi_squared(measured_rv, p.free_parameters)  # store chi squared and p-value in results
            CurveSimResults.sim_rv_plot(p, sim_rv, time_d, "rv_computed")  # plot computed RV
            CurveSimResults.rv_observed_computed_plot(p, sim_rv, time_d, "rv_o_vs_c", measured_rv)  # plot computed and observed RV
            CurveSimResults.rv_residuals_plot(p, "rv_residuals", measured_rv)  # plot RV residuals
        if p.flux_file:
            time_s0, _, _, _, measured_flux = CurveSimResults.get_measured_flux(p)
            for body in bodies:  # HACK because length of body.positions is initialized with the correct value for simulation, NOT measurements
                body.positions = np.ndarray((len(time_s0), 3), dtype=float)
            _, sim_flux, _ = bodies.calc_physics(p, time_s0)  # run simulation
            measured_flux = CurveSimResults.calc_flux_residuals(measured_flux, sim_flux)  # compare observed vs. computed flux
            results.calc_flux_chi_squared(measured_flux, p.free_parameters)  # store chi squared and p-value in results
            CurveSimResults.flux_observed_computed_plots_time(p, "flux_o_vs_c_x=time", measured_flux, measured_tt)  # plot computed and observed flux
            CurveSimResults.flux_observed_computed_plot_data(p, "flux_o_vs_c_x=data", measured_flux)  # plot computed and observed flux
            CurveSimResults.flux_chi_squared_plot_data(p, "flux_chi2_x=data", measured_flux)  # plot flux chi squared per datapoint
            CurveSimResults.flux_residuals_plots_time(p, "flux_residuals_x=time", measured_flux, measured_tt)  # plot Flux residuals
            CurveSimResults.flux_residuals_plot_data(p, "flux_residuals_x=data", measured_flux)  # plot Flux residuals
            # plot something
        if p.tt_file or p.rv_file or p.flux_file:
            results.calc_total_chi_squared(p.free_parameters)
        if p.result_file:
            results.save_results(p)
        if p.sim_flux_file:
            sim_flux.save_sim_flux(p, time_d)
        # self.sim_flux = sim_flux
        # self.results = results
        vitkova_debug = False
        if vitkova_debug:
            p.eclipsers = ["TOI4504d"]
            p.eclipsees = ["TOI4504"]
            # results.plot_parameter("TOI4504c", "TOI4504", "T14", time_d[0], time_d[-1],
            #                         filename=f"TOI4504c_i={bodies[2].i*p.rad2deg:.2f}_T14.png")
            results.plot_parameter("TOI4504d", "TOI4504", "T14", time_d[0], time_d[-1],
                                   filename=f"TOI4504d_i={bodies[1].i * p.rad2deg:.2f}_T14.png")
            results.plot_parameter("TOI4504d", "TOI4504", "depth", time_d[0], time_d[-1],
                                   filename=f"TOI4504d_i={bodies[1].i * p.rad2deg:.2f}_depth.png")

            # measured_tt = CurveSimMCMC.get_measured_tt(p)
            # p.bodynames2bodies(bodies)
            # _, measured_tt = CurveSimMCMC.match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0)
            # dummy_mcmc = CurveSimMCMC(None, None, None, None, None, None, None, dummy_object=True)
            # dummy_mcmc.tt_delta_plot(1, "Vitkova_MaxL_tt_delta.png", measured_tt)
        return bodies, sim_flux, results

    @staticmethod
    def match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0):
        sim_tt = CurveSimBodies.find_tts(rebound_sim, p, sim_flux, time_s0, time_d)  # sim_tt is a list of tuples (eclipser, eclipsee, tt)
        nearest_sim_tt = []
        for idx, row in measured_tt.iterrows():
            eclipser = row["eclipser"]
            measured_tt_val = row["tt"]
            sim_tt_filtered = [tt for tt in sim_tt if tt[0] == eclipser]  # Filter sim_tt for matching eclipser
            if sim_tt_filtered:
                closest_tt = min(sim_tt_filtered, key=lambda x: abs(x[2] - measured_tt_val))  # Find sim_tt with minimal |measured_tt - sim_tt|
                nearest_sim_tt.append(closest_tt[2])
            else:
                nearest_sim_tt.append(0)  # No match found
        measured_tt["nearest_sim"] = nearest_sim_tt  # add 2 columns to data frame
        measured_tt["delta"] = measured_tt["nearest_sim"] - measured_tt["tt"]
        residuals_tt = measured_tt["delta"] / measured_tt["tt_err"]  # residuals are weighted with uncertainty!
        residuals_tt_sum_squared = np.sum(residuals_tt ** 2)
        return residuals_tt_sum_squared, measured_tt

    @staticmethod
    def log_probability(theta, param_bounds, param_references, bodies, time_s0, time_d, measured_flux_array, flux_err, measured_tt, p):
        lp = CurveSimMCMC.log_prior(theta, param_bounds)
        if not np.isfinite(lp):
            return -np.inf
        return lp + CurveSimMCMC.log_likelihood(theta, param_references, bodies, time_s0, time_d, measured_flux_array, flux_err, measured_tt, p)

    @staticmethod
    def log_prior(theta, param_bounds):
        """# If any parameter is outside resonable bounds: return -np.inf"""
        for val, (lower, upper) in zip(theta, param_bounds):
            if not (lower < val < upper):
                return -np.inf
        return 0

    @staticmethod
    def log_likelihood(theta, param_references, bodies, time_s0, time_d, measured_flux_array, flux_err, measured_tt, p):
        # def log_likelihood(theta, param_references, bodies, time_s0, measured_flux_array, flux_err, measured_tt, tt_err, measured_rv, rv_err, p):
        """
        theta:
            List containing the current numerical values of the `param_references` (see below).
            It is automatically modified by the MCMC process.
            Before the simulated lightcurve is recalculated in `log_likelihood()`,
            the parameters are updated using the values from `theta`.

        param_references:
            List containing the names of the parameters to be fitted.
            For example: ["Tmin_pri", "P_days", "incl_deg", "R1a", "R2R1"]
        """
        residuals_sum_squared = 0
        if p.flux_file:
            residuals_sum_squared += p.flux_weight * CurveSimMCMC.residuals_flux_sum_squared(theta, param_references, bodies, time_s0, measured_flux_array, flux_err, p)
        if p.tt_file:
            residuals_sum_squared += p.tt_weight * CurveSimMCMC.residuals_tt_sum_squared(theta, param_references, bodies, time_s0, time_d, measured_tt, p)
            # residuals_sum_squared += p.tt_weight * CurveSimMCMC.residuals_tt_sum_squared_simple(theta, param_references, bodies, time_s0, p)
        # if p.rv_file:
        #     residuals_sum_squared += p.rv_weight * CurveSimMCMC.residuals_rv_sum_squared(theta, param_references, bodies, time_s0, time_d, measured_flux_array, flux_err, p)
        return -0.5 * residuals_sum_squared

    @staticmethod
    def residuals_flux_sum_squared(theta, param_references, bodies, time_s0, measured_flux_array, flux_err, p):
        i = 0
        for body_index, parameter_name in param_references:
            bodies[body_index].__dict__[parameter_name] = theta[i]  # update all parameters from theta
            i += 1
        sim_rv, sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # run simulation
        residuals_flux = (measured_flux_array - sim_flux) / flux_err  # residuals are weighted with uncertainty!
        residuals_flux_sum_squared = np.sum(residuals_flux ** 2)
        return residuals_flux_sum_squared

    @staticmethod
    def bodies_from_fitting_params(bodies, fitting_parameters, param_type=None):
        if param_type == "startvalue":
            for fp in fitting_parameters:
                bodies[fp.body_index].__dict__[fp.parameter_name] = fp.startvalue
        elif param_type == "max_likelihood":
            for fp in fitting_parameters:
                bodies[fp.body_index].__dict__[fp.parameter_name] = fp.max_likelihood / fp.scale
        elif param_type == "mean":
            for fp in fitting_parameters:
                bodies[fp.body_index].__dict__[fp.parameter_name] = fp.mean / fp.scale
        elif param_type == "median":
            for fp in fitting_parameters:
                bodies[fp.body_index].__dict__[fp.parameter_name] = fp.median / fp.scale
        return bodies

    @staticmethod
    def residuals_tt_sum_squared(theta, param_references, bodies, time_s0, time_d, measured_tt, p):
        # measured_tt: pandas DataFrame with columns eclipser, tt, tt_err
        i = 0
        for body_index, parameter_name in param_references:
            bodies[body_index].__dict__[parameter_name] = theta[i]  # update all parameters from theta
            i += 1
        sim_rv, sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # run simulation
        residuals_tt_sum_squared, measured_tt = CurveSimMCMC.match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0)
        return residuals_tt_sum_squared

    @staticmethod
    def residuals_tt_sum_squared_simple(theta, param_references, bodies, time_s0, p):
        """Useful when the config file values of starts and ends are
        chosen so that all simulated flux should be inside transits.
         In that case it is sufficient to merely compare all simulated
         flux to the target flux, which is a single value of about the flux at TT """
        i = 0
        for body_index, parameter_name in param_references:
            bodies[body_index].__dict__[parameter_name] = theta[i]  # update all parameters from theta
            i += 1
        sim_rv, sim_flux, _ = bodies.calc_physics(p, time_s0)  # run simulation

        residuals_tt = sim_flux - p.target_flux
        residuals_tt_sum_squared = np.sum(residuals_tt ** 2)
        return residuals_tt_sum_squared

    @staticmethod
    def hdi_std_mean(data, credible_mass=0.68):
        # Calculate HDI (1-sigma interval with highest density).
        # Data contains samples from the flattened and thinned mcmc-chains and gets sorted.
        # The interval's length is calculated from the number of data and the credible mass percentage.
        # 68% means the interval should contain data from mean minus the standard_deviation til mean plus the standard_deviation.
        # The interval with the smallest difference between its highest (last) and lowest (first) item is chosen
        # Calculate also standard deviation and mean of the sample.
        sorted_data = np.sort(data)
        n = len(sorted_data)
        interval_idx_inc = int(np.floor(credible_mass * n))
        intervals = sorted_data[interval_idx_inc:] - sorted_data[:n - interval_idx_inc]
        min_idx = np.argmin(intervals)
        hdi_min = sorted_data[min_idx]
        hdi_max = sorted_data[min_idx + interval_idx_inc]
        std = np.std(data, ddof=1)
        mean = np.mean(data)
        median = np.median(data)
        return hdi_min, hdi_max, std, mean, median

    def random_initial_values(self):
        """return randomized initial values of the fitting parameters"""
        rng = np.random.default_rng()  # init random number generator
        initial_values = [fp.initial_values(rng, self.walkers) for fp in self.fitting_parameters]
        theta0 = np.array(initial_values)
        return theta0.T

    def scale_samples(self, flat_thin_samples):
        self.scaled_samples = np.copy(flat_thin_samples)
        self.scales = []
        for fpn, ss in zip(self.body_parameter_names, self.scaled_samples.T):
            param = fpn.split(".")[-1]
            ss *= self.scale[param]
            self.scales.append(self.scale[param])

    # @stopwatch()
    def trace_plots(self, steps_done, plot_filename):
        if self.trace_plot_ok:
            try:
                plot_filename = self.results_directory + plot_filename
                fig, axes = plt.subplots(self.ndim, figsize=(10, self.ndim * 2), sharex=True)
                fig.text(0.1, 0.99, f"Traces after {steps_done} steps", ha="left", va="top", fontsize=14, transform=fig.transFigure)
                plt.subplots_adjust(top=0.975)
                if self.ndim == 1:
                    axes = [axes]
                chains = np.moveaxis(self.sampler.get_chain(flat=False, thin=self.thin_samples_plot), -1, 0)
                for i, (chain, ax, name, scale) in enumerate(zip(chains, axes, self.long_body_parameter_names, self.scales)):
                    nsteps = chain.shape[0]
                    x = np.arange(1, nsteps + 1) * self.thin_samples_plot  # 1*thin, 2*thin, ...
                    ax.plot(x, chain * scale, color="xkcd:black", alpha=0.05)
                    ax.set_ylabel(name)
                    ax.axvline(self.burn_in, color="xkcd:tomato", linestyle="solid", label="burn-in")
                    ax.tick_params(labelbottom=True)  # Show x-tick labels for all
                    if i == len(axes) - 1:
                        ax.set_xlabel("Steps including burn-in (red line)")  # Only last subplot
                try:
                    plt.savefig(plot_filename)
                except:
                    print(f"{Fore.RED}\nERROR: Saving Trace Plot failed.{Style.RESET_ALL}")
                plt.close(fig)
            except:
                print(f"{Fore.RED}\nERROR: Trace Plot failed.{Style.RESET_ALL}")
                self.trace_plot_ok = False


    def max_likelihood_parameters(self, flat_thin_samples):
        log_prob_samples = self.sampler.get_log_prob(flat=True, discard=self.burn_in, thin=self.thin_samples)
        if len(log_prob_samples):
            max_likelihood_idx = np.argmax(log_prob_samples)
            self.max_likelihood_params_scaled = self.scaled_samples[max_likelihood_idx]
            self.max_likelihood_params = flat_thin_samples[max_likelihood_idx]
            self.max_log_prob = log_prob_samples[max_likelihood_idx]
        else:
            self.max_likelihood_params_scaled = None
            self.max_likelihood_params = None
            self.max_log_prob = None

    def save_max_likelihood_bodies(self):
        pass

    def high_density_intervals(self):
        # Calculate HDI and other mcmc results.
        self.mean_params = []
        self.median_params = []
        for i, fp in enumerate(self.fitting_parameters):
            hdi_min, hdi_max, std, mean, median = CurveSimMCMC.hdi_std_mean(self.scaled_samples[:, i], self.credible_mass)
            fp.hdi_min = hdi_min
            fp.hdi_max = hdi_max
            fp.std = std
            fp.mean = mean
            fp.median = median
            fp.max_likelihood = self.max_likelihood_params_scaled[i]
            self.mean_params.append(mean)
            self.median_params.append(median)

    def get_max_likelihood_bodies(self, bodies):
        i = 0
        for body_index, parameter_name in self.param_references:
            bodies[body_index].__dict__[parameter_name] = self.max_likelihood_params[i]  # update all parameters from theta
            i += 1
        return bodies

    def max_likelihood_tt(self, max_likelihood_bodies, p, time_s0, time_d, measured_tt):
        sim_rv, sim_flux, rebound_sim = max_likelihood_bodies.calc_physics(p, time_s0)  # run simulation
        _, measured_tt = CurveSimMCMC.match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0)
        return measured_tt

    @staticmethod
    def add_new_best_delta(measured_tt, steps_done):
        first_time = measured_tt.columns[-1] == "delta"
        new_max_likelihood = not measured_tt["delta"].equals(measured_tt.iloc[:, -1])
        if first_time or new_max_likelihood:
            measured_tt[f"step_{steps_done}"] = measured_tt["delta"]
        return measured_tt

    # @stopwatch()
    def mcmc_histograms(self, steps_done, bins, plot_filename):
        plot_filename = self.results_directory + plot_filename
        fig, axes = plt.subplots(self.ndim, figsize=(10, self.ndim * 2))
        fig.text(0.02, 0.99, f"Histograms, {steps_done} steps after burn-in.", ha="left", va="top", fontsize=14, transform=fig.transFigure)
        startvalues = [fp.startvalue * fp.scale for fp in self.fitting_parameters]
        if self.ndim == 1:
            axes = [axes]
        for i, (sample, ax, fp, startvalue) in enumerate(zip(self.scaled_samples.T, axes, self.fitting_parameters, startvalues)):
            densities, bin_edges, _ = ax.hist(sample, bins=bins, density=True, alpha=0.7, color="xkcd:light blue", edgecolor="xkcd:black")
            ax.axvline(fp.hdi_min, color="xkcd:tree green", linestyle="dashed", label="HDI Lower Bound")
            ax.axvline(fp.mean - fp.std, color="xkcd:warm grey", linestyle="dotted", label="Mean - Std")
            ax.axvline(fp.max_likelihood, color="xkcd:tomato", linestyle="solid", label="Max Likelihood")
            ax.axvline(fp.mean, color="xkcd:black", linestyle="solid", label="Mean")
            ax.axvline(fp.hdi_max, color="xkcd:tree green", linestyle="dashed", label="HDI Upper Bound")
            ax.axvline(fp.mean + fp.std, color="xkcd:warm grey", linestyle="dotted", label="Mean + Std")
            ax.axvline(fp.median, color="xkcd:nice blue", linestyle="solid", label="Median")
            ax.axvline(startvalue, color="xkcd:mango", linestyle="solid", label="Startvalue")
            ax.set_xlabel(fp.long_body_parameter_name)
            ax.set_ylabel("Density")
            ax.ticklabel_format(useOffset=False, style="plain", axis="x")  # show x-labels as they are
            if i == 0:
                ax.legend(loc="lower left", bbox_to_anchor=(0.5, 1.02), ncol=3, borderaxespad=0.)
                # ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3, borderaxespad=0.)
        plt.tight_layout()
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}\nERROR: Saving histogram plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    # @stopwatch()
    def mcmc_corner_plot(self, steps_done, plot_filename):
        if self.corner_plot_ok:
            try:
                plot_filename = self.results_directory + plot_filename
                if self.ndim > 1:
                    fig = corner.corner(
                        self.scaled_samples,
                        labels=self.long_body_parameter_names,
                        truths=self.max_likelihood_params_scaled,
                        title_fmt=".4f",
                        quiet=True
                    )
                    fig.suptitle(f"Corner plot. {steps_done} steps after burn-in.", fontsize=16)
                    try:
                        plt.savefig(plot_filename)
                    except:
                        print(f"{Fore.RED}\nERROR: Saving corner plot failed.{Style.RESET_ALL}")
                    plt.close(fig)
            except:
                    print(f"{Fore.RED}\nERROR: Corner plot failed.{Style.RESET_ALL}")
                    self.corner_plot_ok = False

    # @stopwatch()
    def autocorrelation_function_plot(self, steps_done, plot_filename):
        if self.autocorrelation_function_plot_ok:
            try:
                plot_filename = self.results_directory + plot_filename
                samples = self.sampler.get_chain(discard=0, flat=False, thin=self.thin_samples_plot)  # shape: (steps, walkers, ndim)
                nsteps = samples.shape[0]
                nwalkers = samples.shape[1]
                fig, axes = plt.subplots(self.ndim, figsize=(10, self.ndim * 2), sharex=True)
                fig.text(0.1, 0.99, f"Autocorrelation after {steps_done} steps", ha="left", va="top", fontsize=14, transform=fig.transFigure)
                plt.subplots_adjust(top=0.975)
                if self.ndim == 1:
                    axes = [axes]
                x = np.arange(1, nsteps + 1) * self.thin_samples_plot
                for dim, param_name in enumerate(self.long_body_parameter_names):
                    ax = axes[dim]
                    for walker in range(nwalkers):
                        chain_1d = samples[:, walker, dim]
                        ac = np.asarray(emcee.autocorr.function_1d(chain_1d))  # ensure numpy array to avoid `array.pyi` type issue
                        ax.plot(x, ac, alpha=0.2, color="xkcd:royal blue")
                    ax.set_ylabel(param_name)
                    ax.axvline(self.burn_in, color="xkcd:tomato", linestyle="solid", label="burn-in")
                    ax.tick_params(labelbottom=True)  # Show x-tick labels for all
                    if dim == self.ndim - 1:
                        ax.set_xlabel("Steps including burn-in (red line)")  # Only last subplot
                try:
                    plt.savefig(plot_filename)
                except:
                    print(f"{Fore.RED}\nERROR: Saving autocorrelation plot failed.{Style.RESET_ALL}")
                plt.close(fig)
            except:
                print(f"{Fore.RED}\nERROR: Autocorrelation plot failed.{Style.RESET_ALL}")
                self.autocorrelation_function_plot_ok = False

    # @stopwatch()
    def integrated_autocorrelation_time_plot(self, steps_done, plot_filename1, plot_filename2):
        plot_filename1 = self.results_directory + plot_filename1
        plot_filename2 = self.results_directory + plot_filename2
        integrated_autocorrelation_time = np.array(self.integrated_autocorrelation_time).T
        steps = [step for step in range(self.chunk_size + self.loaded_steps, steps_done + 1, self.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = ["xkcd:royal blue", "xkcd:red", "xkcd:black", "xkcd:frog green", "xkcd:piss yellow", "xkcd:purply blue", "xkcd:sepia", "xkcd:wine", "xkcd:ocean", "xkcd:rust", "xkcd:forest", "xkcd:pale violet", "xkcd:robin's egg", "xkcd:pinkish purple", "xkcd:azure", "xkcd:hot pink", "xkcd:mango", "xkcd:baby pink", "xkcd:fluorescent green", "xkcd:medium grey"]
        linestyles = ["solid", "dashed", "dashdot", "dotted"]
        for i, (autocorr_times, fpn) in enumerate(zip(integrated_autocorrelation_time, self.long_body_parameter_names)):
            color = colors[i % len(colors)]
            linestyle = linestyles[(i // len(colors)) % len(linestyles)]
            ax.plot(steps, autocorr_times, label=fpn, color=color, linestyle=linestyle)
        ax.set_xlabel("Steps after burn-in")
        ax.set_title(f"Integrated Autocorrelation Time per Dimension after {steps_done} steps")
        ax.legend(loc="upper left")
        plt.tight_layout()
        try:
            plt.savefig(plot_filename1)
        except:
            print(f"{Fore.RED}\nERROR: Saving Integrated Autocorrelation Time plot failed.{Style.RESET_ALL}")
        plt.close(fig)

        steps_done_div_integrated_autocorrelation_time = steps / integrated_autocorrelation_time
        fig, ax = plt.subplots(figsize=(10, 6))
        for i, (autocorr_times, fpn) in enumerate(zip(steps_done_div_integrated_autocorrelation_time, self.long_body_parameter_names)):
            color = colors[i % len(colors)]
            linestyle = linestyles[(i // len(colors)) % len(linestyles)]
            ax.plot(steps, autocorr_times, label=fpn, color=color, linestyle=linestyle)
        ax.set_xlabel("Steps after burn-in")
        ax.set_title(f"Steps divided by Integrated Autocorrelation Time per Dimension after {steps_done} steps")
        ax.legend(loc="upper left")
        plt.tight_layout()
        try:
            plt.savefig(plot_filename2)
        except:
            print(f"{Fore.RED}\nERROR: Saving Steps divided by Integrated Autocorrelation Time plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    # @stopwatch()
    def acceptance_fraction_plot(self, steps_done, plot_filename):
        if self.acceptance_plot_ok:
            try:
                plot_filename = self.results_directory + plot_filename
                acceptance_fractions_array = np.stack(self.acceptance_fractions, axis=0).T  # shape: (num_lines, 32)
                steps = [step for step in range(self.chunk_size + self.loaded_steps, steps_done + 1, self.chunk_size)]
                fig, ax = plt.subplots(figsize=(10, 6))
                for i in range(acceptance_fractions_array.shape[0]):
                    ax.plot(steps, acceptance_fractions_array[i], label=f"Line {i + 1}", color="xkcd:tree green", alpha=0.15)
                ax.set_xlabel("Steps after burn-in")
                ax.set_ylabel("Acceptance Fraction")
                ax.set_title(f"Acceptance Fraction per Walker after {steps_done} steps")
                plt.tight_layout()
                try:
                    plt.savefig(plot_filename)
                except:
                    print(f"{Fore.RED}\nERROR: Saving acceptance plot failed.{Style.RESET_ALL}")
                plt.close(fig)
            except:
                print(f"{Fore.RED}\nERROR: Acceptance plot failed.{Style.RESET_ALL}")
                self.acceptance_plot_ok = False

    # @stopwatch()
    def average_residual_in_std_plot(self, p, steps_done, plot_filename):
        plot_filename = self.results_directory + plot_filename
        steps = [step for step in range(self.chunk_size + self.loaded_steps, steps_done + 1, self.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(steps, self.max_likelihood_avg_residual_in_std, label="Max Likelihood Parameters", marker="o", markersize=2, color="xkcd:tomato")
        if p.flux_file:
            ax.plot(steps, self.median_avg_residual_in_std, label="Median Parameters", marker="o", markersize=2, color="xkcd:nice blue")
            ax.plot(steps, self.mean_avg_residual_in_std, label="Mean Parameters", marker="o", markersize=2, color="xkcd:black")
        ax.set_xlabel("Steps after burn-in")
        ax.ticklabel_format(useOffset=False, style="plain", axis="y")  # show y-labels as they are
        ax.set_title(f"Average Residual [Standard Deviations] after {steps_done} steps")
        ax.legend(loc="upper left")
        plt.tight_layout()
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}\nERROR: Saving Average Residual plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    def calc_maxlikelihood_avg_residual_in_std(self, p):
        flux, rv, tt = 0, 0, 0
        if p.flux_file:
            flux = getattr(p, "total_iterations", 0)
        if p.rv_file:
            rv = getattr(p, "rv_datasize", 0)
        if p.tt_file:
            tt = getattr(p, "tt_datasize", 0)
        maxlikelihood_avg_residual_in_std = math.sqrt(-2 * self.max_log_prob / (flux + rv + tt))
        self.max_likelihood_avg_residual_in_std.append(maxlikelihood_avg_residual_in_std)

    # @stopwatch()
    @staticmethod
    def tt_delta_plot(p, steps_done, plot_filename, measured_tt):
        plot_filename = p.results_directory + plot_filename
        unique_eclipsers = measured_tt["eclipser"].unique()
        n_eclipsers = len(unique_eclipsers)
        fig, axes = plt.subplots(n_eclipsers, figsize=(10, 3.5 * n_eclipsers), sharex=True)
        if n_eclipsers == 1:
            axes = [axes]
        abs_min = abs(min(measured_tt["delta"].min(), 0))
        abs_max = abs(max(measured_tt["delta"].max(), 0))
        ylim = (-1.3 * max(abs_min, abs_max), 1.3 * max(abs_min, abs_max))
        dx = 0.005 * (measured_tt["tt"].max() - measured_tt["tt"].min())
        for ax, eclipser in zip(axes, unique_eclipsers):
            df = measured_tt[measured_tt["eclipser"] == eclipser]
            ax.plot(df["tt"], df["delta"], marker="o", markersize=5, linestyle="-", color="xkcd:nice blue")
            for x, tt_err in zip(df["tt"], df["tt_err"]):  # uncertainties
                ax.hlines(tt_err, x - dx, x + dx, colors="xkcd:tomato", linewidth=1)
                ax.hlines(-tt_err, x - dx, x + dx, colors="xkcd:tomato", linewidth=1)
                ax.vlines(x, -tt_err, tt_err, colors="xkcd:tomato", linewidth=1)
            ax.axhline(0, color="xkcd:warm grey", linestyle="dashed", linewidth=1)
            ax.set_ylabel(f"Transit Time Residuals [days]")
            ax.set_title(f"Eclipser: {eclipser}")
            ax.tick_params(labelbottom=True)
            ax.set_ylim(ylim)
            ax.ticklabel_format(useOffset=False, style="plain", axis="x")
        axes[-1].set_xlabel("Transit Time [BJD]")
        if steps_done > 0:
            fig.suptitle(f"TT Delta. {steps_done} steps after burn-in.", fontsize=14)
        else:
            fig.suptitle(f"TT Delta", fontsize=14)
        plt.tight_layout(rect=(0, 0, 1, 0.97))
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}\nERROR: Saving TT delta plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    # @stopwatch()
    def tt_multi_delta_plot(self, steps_done, plot_filename, measured_tt):
        plot_filename = self.results_directory + plot_filename
        # plot_filename = self.results_directory + str(steps_done) + plot_filename
        unique_eclipsers = measured_tt["eclipser"].unique()
        n_eclipsers = len(unique_eclipsers)
        fig, axes = plt.subplots(n_eclipsers, figsize=(10, 4 * n_eclipsers), sharex=True)
        if n_eclipsers == 1:
            axes = [axes]

        delta_columns = [col for col in measured_tt.columns if col.startswith("step_")]  # Find up to 10 last columns with names starting with "step_"
        delta_columns = delta_columns[-10:]  # last up to 10

        abs_min = abs(min(measured_tt["delta"].min(), 0))
        abs_max = abs(max(measured_tt["delta"].max(), 0))
        ylim = (-1.3 * max(abs_min, abs_max), 1.3 * max(abs_min, abs_max))

        num_lines = len(delta_columns)
        if num_lines > 2:  # Create the gray scale colors from light to dark for all but the last two lines
            gray_shades = plt.cm.gray(np.linspace(0.85, 0.15, num_lines - 2))
        else:
            gray_shades = []

        for ax, eclipser in zip(axes, unique_eclipsers):
            df = measured_tt[measured_tt["eclipser"] == eclipser]
            for i, col in enumerate(delta_columns):
                if i == num_lines - 1:
                    color = "xkcd:nice blue"  # last line blue
                elif i == num_lines - 2:
                    color = "xkcd:black"  # second to last black
                else:
                    color = gray_shades[i]  # shades of gray
                ax.plot(df["tt"], df[col], marker="o", linestyle="-", alpha=0.7, label=col, color=color)
            ax.axhline(0, color="xkcd:warm grey", linestyle="dashed", linewidth=1)
            ax.set_ylabel(f"TT Delta [days]")
            ax.set_title(f"Eclipser: {eclipser}")
            ax.tick_params(labelbottom=True)
            ax.set_ylim(ylim)

        axes[-1].set_xlabel("Transit Time [BJD]")
        fig.suptitle(f"TT Delta. {steps_done} steps after burn-in.", fontsize=14)
        plt.tight_layout(rect=(0, 0.12, 1, 0.97))  # leave 12% at bottom
        fig.canvas.draw()  # Draw the figure first to avoid clipping bugs
        handles, labels = axes[0].get_legend_handles_labels()  # Since all subplots have identical handles/labels, get from only first subplot
        fig.legend(handles, labels,
                   bbox_to_anchor=(0.175, 0.01, 0.65, 0.1),  # x, y, width, height
                   ncol=5, mode="expand", fontsize="small", frameon=False)
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}\nERROR: Saving TT delta plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    @staticmethod
    def seconds2readable(seconds):
        days = int(seconds // 86400)
        hours = int((seconds % 86400) // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{days:02d}:{hours:02d}:{minutes:02d}:{secs:02.0f} [dd:hh:mm:ss]"

    def save_mcmc_results(self, p, bodies, steps_done, measured_tt):
        results = {}
        results["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        results["Simulation Parameters"] = {}
        results["Simulation Parameters"]["comment"] = getattr(p, "comment", None)

        results["Simulation Parameters"]["start_realtime"] = self.start_real_time + " [DD.MM.YY hh:mm:ss]"
        results["Simulation Parameters"]["end_realtime"] = time.strftime("%d.%m.%y %H:%M:%S") + " [DD.MM.YY hh:mm:ss]"
        runtime = time.perf_counter() - self.start_timestamp
        results["Simulation Parameters"]["run_time"] = CurveSimMCMC.seconds2readable(runtime)
        results["Simulation Parameters"]["run_time_per_iteration"] = f"{runtime / (self.burn_in + steps_done):.3f} [s]"
        results["Simulation Parameters"]["simulations_per_second"] = f"{(self.burn_in + steps_done) * self.walkers / runtime:.0f} [iterations*walkers/runtime]"

        results["Simulation Parameters"]["results_directory"] = self.results_directory
        results["Simulation Parameters"]["start_date"] = p.start_date
        results["Simulation Parameters"]["default_dt"] = p.dt
        results["Simulation Parameters"]["flux_data_points"] = getattr(p, "total_iterations", None)
        results["Simulation Parameters"]["walkers"] = self.walkers
        results["Simulation Parameters"]["burn_in_steps"] = self.burn_in
        results["Simulation Parameters"]["steps_after_burn_in"] = int(steps_done)
        results["Simulation Parameters"]["moves"] = p.moves
        results["Simulation Parameters"]["thin_samples"] = self.thin_samples

        if p.flux_file:
            results["Simulation Parameters"]["flux_file"] = p.flux_file
            results["Simulation Parameters"]["mean_avg_residual_in_std"] = self.mean_avg_residual_in_std[-1]
            results["Simulation Parameters"]["median_avg_residual_in_std"] = self.median_avg_residual_in_std[-1]
        if p.tt_file:
            results["Simulation Parameters"]["tt_file"] = p.tt_file
            results["Simulation Parameters"]["tt_data_points"] = p.tt_datasize
        # if p.rv_file:
            # results["Simulation Parameters"]["rv_file"] = p.rv_file

        results["Simulation Parameters"]["max_log_prob"] = self.max_log_prob
        results["Simulation Parameters"]["max_likelihood_avg_residual_in_std"] = self.max_likelihood_avg_residual_in_std[-1]

        if p.tt_file:
            results["Simulation Parameters"]["mean_delta"] = float(np.mean(np.abs(measured_tt["delta"])))
            results["Simulation Parameters"]["max_delta"] = float(np.max(np.abs(measured_tt["delta"])))
            results["Simulation Parameters"]["param_json"] = bodies.bodies2param_json(measured_tt, p)
            results["measured_tt_list"] = measured_tt.to_dict(orient="list")  # Convert measured_tt DataFrame to a serializable format

        results["Bodies"] = {}
        params = (["body_type", "primary", "mass", "radius", "luminosity"]
                  + ["limb_darkening_u1", "limb_darkening_u2", "mean_intensity", "intensity"]
                  + ["e", "i", "P", "a", "Omega", "Omega_deg", "omega", "omega_deg", "pomega", "pomega_deg"]
                  + ["L", "L_deg", "ma", "ma_deg", "ea", "ea_deg", "nu", "nu_deg", "T", "t"])
        fitting_param_tuples = [(fp.body_index, fp.parameter_name) for fp in self.fitting_parameters]
        for i, body in enumerate(bodies):
            results["Bodies"][body.name] = {}
            for key in params:
                if (i, key) not in fitting_param_tuples and (i, key.split("_deg")[0]) not in fitting_param_tuples:
                    attr = getattr(body, key)
                    if attr is not None:
                        results["Bodies"][body.name][key] = attr

        fitting_parameters = copy.deepcopy(p.fitting_parameters)
        for fp in fitting_parameters:
            fp.startvalue *= fp.scale
            fp.lower *= fp.scale
            fp.upper *= fp.scale
            fp.sigma *= fp.scale
        results["Fitting Parameters"] = {fp.body_parameter_name: fp.__dict__ for fp in fitting_parameters}

        p_copy = copy.deepcopy(p)
        to_remove = [
            "fitting_parameters", "standard_sections", "eclipsers", "eclipsees",
            "tt_file", "total_iterations", "walkers", "moves", "burn_in",
            "thin_samples", "comment", "start_date", "results_directory",
            "fitting_parameter_dic", "tt_datasize",
        ]
        for name in to_remove:
            if hasattr(p_copy, name):
                delattr(p_copy, name)

        for name in ("starts_s0", "starts_d", "ends_s0", "ends_d", "dts"):
            if hasattr(p_copy, name):
                orig = getattr(p_copy, name)
                p_copy.__dict__[name] = [float(i) for i in orig]

        results["ProgramParameters"] = p_copy.__dict__
        self.mcmc_results2json(results, p)

    def mcmc_results2json(self, results, p):
        """Converts results to JSON and saves it."""
        filename = self.results_directory + "mcmc_results.json"
        try:
            with open(filename, "w", encoding="utf8") as file:
                json.dump(results, file, indent=4, ensure_ascii=False)
            if p.verbose:
                print(f" Saved MCMC results to {filename}")
        except:
            print(f"{Fore.RED}\nERROR: Saving MCMC Results JSON failed.{Style.RESET_ALL}")
            print(results)
            print(f"{Fore.YELLOW}Printed Results to console because saving failed.{Style.RESET_ALL}")

    @stopwatch()
    def mcmc_results(self, p, bodies, steps_done, time_s0, time_d, measured_tt, measured_flux_array, flux_err, chunk):
        flat_thin_samples = self.sampler.get_chain(discard=self.burn_in, thin=self.thin_samples, flat=True)
        # discard the initial self.burn_in steps from each chain to ensure only samples that represent the equilibrium distribution are analyzed.
        # thin=10: keep only every 10th sample from the chain to reduce autocorrelation in the chains and the size of the resulting arrays.
        # flat=True: return all chains in a single, two-dimensional array (shape: (n_samples, n_parameters))
        print(f"{steps_done} steps done.  ")

        self.acceptance_fractions.append(self.sampler.acceptance_fraction)
        if chunk % 5 == 0:
            self.acceptance_fraction_plot(steps_done, "acceptance.png")
        self.scale_samples(flat_thin_samples)
        self.max_likelihood_parameters(flat_thin_samples)
        # self.save_max_likelihood_bodies()

        if p.tt_file:
            max_likelihood_bodies = self.get_max_likelihood_bodies(bodies)
            measured_tt = self.max_likelihood_tt(max_likelihood_bodies, p, time_s0, time_d, measured_tt)
            measured_tt = CurveSimMCMC.add_new_best_delta(measured_tt, steps_done)
            CurveSimMCMC.tt_delta_plot(p, steps_done, "tt_delta.png", measured_tt)
            self.tt_multi_delta_plot(steps_done, "tt_multi_delta.png", measured_tt)
        self.calc_maxlikelihood_avg_residual_in_std(p)
        self.high_density_intervals()

        if p.flux_file:
            median_residuals_flux_sum_squared = CurveSimMCMC.residuals_flux_sum_squared(self.median_params, self.param_references, bodies, time_s0, measured_flux_array, flux_err, p)
            mean_residuals_flux_sum_squared = CurveSimMCMC.residuals_flux_sum_squared(self.mean_params, self.param_references, bodies, time_s0, measured_flux_array, flux_err, p)
            flux_data_points = getattr(p, "total_iterations", 0)
            self.mean_avg_residual_in_std.append(math.sqrt(mean_residuals_flux_sum_squared / flux_data_points))
            self.median_avg_residual_in_std.append(math.sqrt(median_residuals_flux_sum_squared / flux_data_points))
        self.average_residual_in_std_plot(p, steps_done, "avg_residual.png")

        bodies = CurveSimMCMC.bodies_from_fitting_params(bodies, self.fitting_parameters, param_type="max_likelihood")
        bodies.save(prefix=p.comment, suffix="_maxL")
        CurveSimParameters.save_fitting_parameters(self.fitting_parameters, prefix=p.comment, suffix="_maxL")
        CurveSimMCMC.single_run(p, bodies, time_s0, time_d)

        self.integrated_autocorrelation_time.append(list(self.sampler.get_autocorr_time(tol=0)))
        self.integrated_autocorrelation_time_plot(steps_done, "int_autocorr_time.png", "steps_per_i_ac_time.png")
        if chunk % 10 == 0:
            self.autocorrelation_function_plot(steps_done, "autocorrelation.png")

        for bins in self.bins:
            self.mcmc_histograms(steps_done, bins, f"histograms_{bins}.png")

        self.save_mcmc_results(p, bodies, steps_done, measured_tt)
        if chunk % 5 == 0:
            self.trace_plots(steps_done, "traces.png")
        if chunk % 10 == 0:
            flat_thin_samples = self.sampler.get_chain(discard=self.burn_in, thin=self.thin_samples_plot, flat=True)
            self.scale_samples(flat_thin_samples)
            self.mcmc_corner_plot(steps_done, "corner.png")


class CurveSimLMfit:

    def __init__(self, p, bodies, time_s0, time_d, measured_tt):
        if not (p.flux_file or p.tt_file or p.rv_file):
            print(f"{Fore.RED}\nERROR: No measurements for fitting hve been provided.{Style.RESET_ALL}")
            sys.exit(1)
        if os.path.exists("residual.tmp"):
            os.remove("residual.tmp")
        if os.path.exists("iteration.tmp"):
            os.remove("iteration.tmp")
        self.results_directory = p.results_directory
        self.fitting_parameters = p.fitting_parameters
        self.unit = p.unit
        self.scale = p.scale
        self.param_references = [(fp.body_index, fp.parameter_name) for fp in self.fitting_parameters]  # list of names of fitting parameters. Needed so these parameters can be updated inside log_likelihood().
        self.body_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in self.fitting_parameters]
        p.index_from_bodyparamname = {bpn: fp.index for bpn, fp in zip(self.body_parameter_names, self.fitting_parameters)}
        self.long_body_parameter_names = [fpn + " [" + self.unit[fpn.split(".")[-1]] + "]" for fpn in self.body_parameter_names]
        for fp, fpn, fpnu in zip(p.fitting_parameters, self.body_parameter_names, self.long_body_parameter_names):
            fp.body_parameter_name = fpn
            fp.long_body_parameter_name = fpnu
        self.param_bounds = [(fp.lower, fp.upper) for fp in self.fitting_parameters]
        self.args = (self.param_references, bodies, time_s0, time_d, measured_tt, p)
        self.start_real_time = time.strftime("%d.%m.%y %H:%M:%S")
        self.start_timestamp = time.perf_counter()
        for (body_index, parameter_name), fp in zip(self.param_references, p.fitting_parameters):  # update bodies from fitting parameters (in case of changed fitting parameter start values)
            setattr(bodies[body_index], parameter_name, fp.startvalue)
        self.params = lmfit.Parameters()
        for (body_index, parameter_name), (lower, upper) in zip(self.param_references, self.param_bounds):
            self.params.add(bodies[body_index].name + "_" + parameter_name, value=bodies[body_index].__dict__[parameter_name], min=lower, max=upper)

        # self.result = lmfit.minimize(CurveSimLMfit.lmfit_residual_tt, self.params, method="brute", args=(self.param_references, bodies, time_s0, time_d, measured_tt, p))
        self.result = lmfit.minimize(CurveSimLMfit.lmfit_residual_tt, self.params, method=p.lmfit_method, args=(self.param_references, bodies, time_s0, time_d, measured_tt, p))
        # self.result = lmfit.minimize(CurveSimLMfit.lmfit_residual_tt, self.params, method="nelder", args=(self.param_references, bodies, time_s0, time_d, measured_tt, p))
        # self.result = lmfit.minimize(CurveSimLMfit.lmfit_residual_tt, self.params, method="powell", args=(self.param_references, bodies, time_s0, time_d, measured_tt, p))
        # ***** METHODS ******
        # best?                                     powell: Powell’s method
        # fine                                      nelder: Nelder-Mead simplex
        # fine                                      differential_evolution: Differential Evolution (global optimization)
        # n params, n<8 because 7.5*n*20**n bytes memory, slow  brute: Brute force grid search

        # needs Jacobian                            newton: Newton-CG
        # needs Jacobian                            dogleg	Dogleg
        # needs Jacobian                            trust-exact	Exacttrust-region
        # needs Jacobian                            trust-krylov	NewtonGLTRtrust-region
        # needs Jacobian                            trust-ncg	NewtonCGtrust-region
        # does not even find minimum for 3 params   least_squares: SciPy’s least_squares (Trust Region Reflective, Dogbox, Levenberg-Marquardt)
        # does not even find minimum for 3 params   lbfgsb: L-BFGS-B (bounded minimization)
        # does not even find minimum for 3 params   ampgo: Adaptive Memory Programming for Global Optimization
        # does not even find minimum for 3 params   cg: Conjugate Gradient
        # does not even find minimum for 3 params   cobyla: COBYLA
        # does not even find minimum for 3 params   bfgs: BFGS
        # does not even find minimum for 3 params   tnc: Truncated Newton
        # does not even find minimum for 3 params   trust-constr: Trust Region Constrained
        # does not even find minimum for 3 params   basinhopping	Basinhopping
        # does not even find minimum for 3 params   dual_annealing	DualAnnealing
        # does not even find minimum for 3 params   shgo	SimplicialHomologyGlobalOptimization
        # does not even find minimum for 3 params   slsqp	SequentialLinearSquaresProgramming
        # does not find minimum + needs more residual than params  leastsq: Levenberg-Marquardt (default, for least-squares problems)

    @staticmethod
    def lmfit_residual_tt(params, param_references, bodies, time_s0, time_d, measured_tt, p):
        # measured_tt: pandas DataFrame with columns eclipser, tt, tt_err
        for body_index, parameter_name in param_references:
            bodies[body_index].__dict__[parameter_name] = params[bodies[body_index].name + "_" + parameter_name].value  # update all parameters from params
        sim_rv, sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # run simulation
        residuals_tt_sum_squared, measured_tt = CurveSimMCMC.match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0)
        return residuals_tt_sum_squared

    def save_lmfit_results(self, p):
        results = {}
        results["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        results["Simulation Parameters"] = {}
        results["Simulation Parameters"]["comment"] = getattr(p, "comment", None)

        results["Simulation Parameters"]["start_realtime"] = self.start_real_time + " [DD.MM.YY hh:mm:ss]"
        results["Simulation Parameters"]["end_realtime"] = time.strftime("%d.%m.%y %H:%M:%S") + " [DD.MM.YY hh:mm:ss]"
        runtime = time.perf_counter() - self.start_timestamp
        results["Simulation Parameters"]["run_time"] = CurveSimMCMC.seconds2readable(runtime)

        results["Simulation Parameters"]["results_directory"] = self.results_directory

        if p.flux_file:
            results["Simulation Parameters"]["flux_file"] = p.flux_file
            results["Simulation Parameters"]["mean_avg_residual_in_std"] = self.mean_avg_residual_in_std[-1]
            results["Simulation Parameters"]["median_avg_residual_in_std"] = self.median_avg_residual_in_std[-1]
        if p.tt_file:
            results["Simulation Parameters"]["tt_file"] = p.tt_file
            results["Simulation Parameters"]["tt_data_points"] = p.tt_datasize
            # results["Simulation Parameters"]["rv_file"] = p.rv_file

        result_copy = copy.deepcopy(self.result)
        result_copy.last_internal_values = list(result_copy.last_internal_values)
        result_copy.residual = list(result_copy.residual)
        result_copy.x = list(result_copy.x)
        result_copy.params = json.loads(result_copy.params.dumps())

        # results["LMfitParameters"] = find_ndarrays(result_copy.__dict__)

        self.lmfit_results2json(results, p)

    @staticmethod
    def check_for_fit_improvement(residual):
        try:
            with open("residual.tmp", "r", encoding="utf8") as file:
                best_residual = float(file.read().strip())
        except (FileNotFoundError, ValueError):
            best_residual = float("inf")
        improvement = residual < best_residual
        if improvement:
            with open("residual.tmp", "w", encoding="utf8") as file:
                file.write(str(residual))
        return improvement

    @staticmethod
    def get_iteration_from_file():
        try:
            with open("iteration.tmp", "r", encoding="utf8") as file:
                iteration = int(file.read().strip())
        except (FileNotFoundError, ValueError):
            iteration = 0
        with open("iteration.tmp", "w", encoding="utf8") as file:
            file.write(str(iteration + 1))
        return iteration

    @staticmethod
    def save_intermediate_lmfit_results(p, bodies, measured_tt):
        results = {}
        results["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        results["Simulation Parameters"] = {}
        results["Simulation Parameters"]["comment"] = getattr(p, "comment", None)
        results["Simulation Parameters"]["end_realtime"] = time.strftime("%d.%m.%y %H:%M:%S") + " [DD.MM.YY hh:mm:ss]"

        if p.flux_file:
            results["Simulation Parameters"]["flux_file"] = p.flux_file
        if p.tt_file:
            results["Simulation Parameters"]["tt_file"] = p.tt_file
            results["Simulation Parameters"]["tt_data_points"] = p.tt_datasize
            # results["Simulation Parameters"]["rv_file"] = p.rv_file
        # if p.tt_file:
        #     results["Simulation Parameters"]["tt_measured"] = list(p.best_tt_df["tt"])
        #     results["Simulation Parameters"]["tt_best_sim"] = list(p.best_tt_df["nearest_sim"])

        results["Bodies"] = {}
        params = (["body_type", "primary", "mass", "radius", "luminosity"]
                  + ["limb_darkening_u1", "limb_darkening_u2", "mean_intensity", "intensity"]
                  + ["e", "i", "P", "a", "Omega", "omega", "pomega"]
                  + ["L", "ma", "ea", "ea_deg", "nu", "T", "t"])

        for i, body in enumerate(bodies):
            results["Bodies"][body.name] = {}
            for key in params:
                attr = getattr(body, key)
                if attr is not None:
                    if key in p.scale:
                        scale = p.scale[key]
                    else:
                        scale = 1
                    results["Bodies"][body.name][key] = attr * scale

        fitting_parameters = copy.deepcopy(p.fitting_parameters)
        for fp in fitting_parameters:
            fp.startvalue *= fp.scale
            fp.lower *= fp.scale
            fp.upper *= fp.scale
            fp.last_value = bodies[fp.body_index].__dict__[fp.parameter_name]
            fp.last_value *= fp.scale
            # width = 25 - len(fp.body_parameter_name)
            # print(f"{fp.body_parameter_name}:{fp.last_value:{width}.5f}")

        results["Fitting Parameters"] = {fp.body_parameter_name: fp.__dict__ for fp in fitting_parameters}

        results["measured_tt_list"] = measured_tt.to_dict(orient="list")  # Convert measured_tt DataFrame to a serializable format
        # results["measured_tt_records"] = measured_tt.to_dict(orient="records")  # Convert measured_tt DataFrame to a serializable format

        p_copy = copy.deepcopy(p)
        del p_copy.fitting_parameters
        del p_copy.standard_sections
        del p_copy.eclipsers
        del p_copy.eclipsees
        del p_copy.tt_file
        del p_copy.total_iterations
        del p_copy.walkers
        del p_copy.moves
        del p_copy.burn_in
        del p_copy.thin_samples
        del p_copy.tt_datasize
        del p_copy.comment
        del p_copy.start_date
        del p_copy.results_directory
        del p_copy.starts_s0
        del p_copy.starts_d
        del p_copy.ends_s0
        del p_copy.ends_d
        del p_copy.dts
        results["ProgramParameters"] = p_copy.__dict__

        filename = p.results_directory + f"/lmfit_results.tmp.json"
        with open(filename, "w", encoding="utf8") as file:
            json.dump(results, file, indent=4, ensure_ascii=False)
        if p.verbose:
            print(f" Saved intermediate LMfit results to {filename}")

    def lmfit_results2json(self, results, p):
        """Converts results to JSON and saves it."""
        filename = self.results_directory + f"/lmfit_results.json"
        with open(filename, "w", encoding="utf8") as file:
            json.dump(results, file, indent=4, ensure_ascii=False)
        if p.verbose:
            print(f" Saved LMfit results to {filename}")

    def save_best_fit(self, p, bodies, measured_tt):
        max_delta = float(np.max(np.abs(measured_tt["delta"])))
        mean_delta = float(np.mean(np.abs(measured_tt["delta"])))
        color = Fore.WHITE
        if mean_delta < 1.0:
            color = Fore.RED
            if mean_delta < 0.1:
                color = Fore.YELLOW
            if mean_delta < 0.02:
                color = Fore.GREEN
            if mean_delta < 0.004:
                color = Fore.CYAN
            if mean_delta < 0.0024 or max_delta < 0.0048:
                color = Fore.MAGENTA
            line = bodies.bodies2param_json(measured_tt, p)
            filename = p.results_directory + f"/lmfit_best_fits.txt"
            try:
                append_line_locked(filename, line, wait=0.1)
            except OSError as e:
                # non-fatal: print error but continue
                print(f"{Fore.RED}\nERROR: Could not write best fit to `lmfit_best_fits.txt`: {e}{Style.RESET_ALL}")
        runtime = CurveSimMCMC.seconds2readable(time.perf_counter() - self.start_timestamp)
        print(f"{color}Runtime: {runtime}   max_delta: {max_delta:7.4f} days  mean_delta: {mean_delta:7.4f} days{Style.RESET_ALL}    ", end="")


def find_ndarrays(obj, path="root"):
    if isinstance(obj, np.ndarray):
        print(f"{path}: numpy.ndarray, shape={obj.shape}, dtype={obj.dtype}")
        return obj.tolist()
    elif isinstance(obj, dict):
        for k, v in obj.items():
            obj[k] = find_ndarrays(v, f"{path}[{repr(k)}]")
        return obj
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            obj[i] = find_ndarrays(v, f"{path}[{i}]")
        return obj
    elif hasattr(obj, "__dict__"):
        obj.__dict__ = find_ndarrays(obj.__dict__, f"{path}.__dict__")
        return obj
    else:
        return obj


class CurveSimParameters:

    def __init__(self, config_file):
        """Read program parameters and properties of the physical bodies from config file."""
        self.body_parameter_names = None
        self.long_body_parameter_names = None
        self.PARAMS = (["body_type", "primary", "mass", "radius", "luminosity"]
                       + ["limb_darkening_u1", "limb_darkening_u2", "mean_intensity", "intensity"]
                       + ["e", "i", "P", "a", "Omega", "omega", "pomega"]
                       + ["L", "ma", "ea", "nu", "T", "t"])
        self.standard_sections = ["Astronomical Constants", "Results", "Simulation", "Fitting", "Video", "Plot", "Scale", "Debug"]  # These sections must be present in the config file.
        config = configparser.ConfigParser(inline_comment_prefixes="#")  # Inline comments in the config file start with "#".
        config.optionxform = str  # Preserve case of the keys.
        CurveSimParameters.find_and_check_config_file(config_file)  #, standard_sections=self.standard_sections)
        config.read(config_file, encoding="utf-8")
        self.config_file = config_file

        # [Astronomical Constants]
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
        g = eval(config.get("Astronomical Constants", "g", fallback="None"))
        au = eval(config.get("Astronomical Constants", "au", fallback="None"))
        r_sun = eval(config.get("Astronomical Constants", "r_sun", fallback="None"))
        m_sun = eval(config.get("Astronomical Constants", "m_sun", fallback="None"))
        l_sun = eval(config.get("Astronomical Constants", "l_sun", fallback="None"))
        r_jup = eval(config.get("Astronomical Constants", "r_jup", fallback="None"))
        m_jup = eval(config.get("Astronomical Constants", "m_jup", fallback="None"))
        r_earth = eval(config.get("Astronomical Constants", "r_earth", fallback="None"))
        m_earth = eval(config.get("Astronomical Constants", "m_earth", fallback="None"))
        v_earth = eval(config.get("Astronomical Constants", "v_earth", fallback="None"))
        hour = eval(config.get("Astronomical Constants", "hour", fallback="None"))
        day = eval(config.get("Astronomical Constants", "day", fallback="None"))
        year = eval(config.get("Astronomical Constants", "year", fallback="None"))
        rad2deg = eval(config.get("Astronomical Constants", "rad2deg", fallback="None"))
        self.g, self.au, self.r_sun, self.m_sun, self.l_sun = g, au, r_sun, m_sun, l_sun,
        self.r_jup, self.m_jup, self.r_earth, self.m_earth, self.v_earth = r_jup, m_jup, r_earth, m_earth, v_earth
        self.hour, self.day, self.year, self.rad2deg = hour, day, year, rad2deg

        # [Results]
        self.comment = config.get("Results", "comment", fallback="No comment")
        self.verbose = eval(config.get("Results", "verbose", fallback="False"))
        self.transit_precision = eval(config.get("Results", "transit_precision", fallback="1"))

        self.flux_data_directory = config.get("Results", "flux_data_directory", fallback="")

        # [Simulation]
        self.action = config.get("Simulation", "action", fallback="results_only")
        # self.single_run = eval(config.get("Simulation", "single_run", fallback="False"))
        # self.results_only = eval(config.get("Simulation", "results_only", fallback="False"))
        self.dt = eval(config.get("Simulation", "dt"))
        self.start_date = eval(config.get("Simulation", "start_date", fallback="0.0"))

        # [Video]
        self.video_file = config.get("Video", "video_file", fallback=None)

        # [Fitting]
        self.mcmc_multi_processing = eval(config.get("Fitting", "mcmc_multi_processing", fallback="True"))
        self.free_parameters = eval(config.get("Fitting", "free_parameters", fallback="None"))
        self.flux_file = config.get("Fitting", "flux_file", fallback=None)
        self.tt_file = config.get("Fitting", "tt_file", fallback=None)
        self.rv_file = config.get("Fitting", "rv_file", fallback=None)
        self.rv_body = config.get("Fitting", "rv_body", fallback=None)
        self.eclipsers_names = list([x.strip() for x in config.get("Fitting", "eclipsers_names", fallback="None").split("#")[0].split(",")])
        self.eclipsees_names = list([x for x in config.get("Fitting", "eclipsees_names", fallback="None").split("#")[0].split(",")])

        # [Results]
        self.result_file = config.get("Results", "result_file", fallback="None")
        if self.result_file == "None":
            self.result_file = None
        self.result_dt = eval(config.get("Results", "result_dt", fallback="100"))
        self.max_interval_extensions = eval(config.get("Results", "max_interval_extensions", fallback="10"))

        # [Simulation]
        self.results_directory = config.get("Simulation", "results_directory", fallback="None")
        if self.results_directory == "None":
            self.results_directory = None
        if self.results_directory is not None:
            self.find_results_subdirectory()
        self.starts_d = np.array(eval(config.get("Simulation", "starts", fallback="[]")), dtype=float)
        self.ends_d = np.array(eval(config.get("Simulation", "ends", fallback="[]")), dtype=float)
        self.dts = np.array(eval(config.get("Simulation", "dts", fallback="[]")), dtype=float)

        self.sim_flux_file = config.get("Simulation", "sim_flux_file", fallback=None)
        self.start_indices, self.max_iterations, self.total_iterations = self.check_intervals()
        self.tt_padding = eval(config.get("Fitting", "tt_padding", fallback="0.3"))
        self.flux_plots_top = eval(config.get("Fitting", "flux_plots_top", fallback="1.015"))
        self.flux_plots_bottom = eval(config.get("Fitting", "flux_plots_bottom", fallback="0.97"))

        if self.action == "single_run":  # run simulation, generate video and transit results
            if self.sim_flux_file == "None":
                self.sim_flux_file = None
            self.sim_flux_err = eval(config.get("Simulation", "sim_flux_err", fallback="0.0"))

            # [Video]
            self.frames = eval(config.get("Video", "frames"))
            self.fps = eval(config.get("Video", "fps"))
            self.clockwise = eval(config.get("Video", "clockwise", fallback="False"))

            # self.start_indices, self.max_iterations, self.total_iterations = self.check_intervals()
            self.sampling_rate = self.total_iterations / self.frames
            # self.sampling_rate = (self.total_iterations - 1) // self.frames + 1
            if self.sampling_rate < 1:
                print(f"{Fore.YELLOW}\nWARNING: This simulation calculates only {self.total_iterations} iterations for {self.frames} video frames.{Style.RESET_ALL}")
                print(f"{Fore.YELLOW}         Because of this undersampling, the video will be using the same iteration for consecutive frames.{Style.RESET_ALL}")
                print(f"{Fore.YELLOW}         Decrease the number of frames or decrease dt (the real time difference between simulation iterations).{Style.RESET_ALL}")
            # [Scale]
            self.scope_left = eval(config.get("Scale", "scope_left"))
            self.scale_bar_length_left = eval(config.get("Scale", "scale_bar_length_left"))
            self.star_scale_left = eval(config.get("Scale", "star_scale_left"))
            self.planet_scale_left = eval(config.get("Scale", "planet_scale_left"))
            self.scope_right = eval(config.get("Scale", "scope_right"))
            self.scale_bar_length_right = eval(config.get("Scale", "scale_bar_length_right"))
            self.star_scale_right = eval(config.get("Scale", "star_scale_right"))
            self.planet_scale_right = eval(config.get("Scale", "planet_scale_right"))
            self.autoscaling = config.get("Scale", "autoscaling") == "on"
            self.min_radius = eval(config.get("Scale", "min_radius")) / 100.0
            self.max_radius = eval(config.get("Scale", "max_radius")) / 100.0

            # [Plot]
            self.show_left_plot = eval(config.get("Plot", "show_left_plot", fallback="True"))
            self.show_right_plot = eval(config.get("Plot", "show_right_plot", fallback="True"))
            self.show_lc_plot = eval(config.get("Plot", "show_lc_plot", fallback="True"))
            self.show_rv_plot = eval(config.get("Plot", "show_rv_plot", fallback="True"))
            self.main_title = config.get("Plot", "main_title", fallback="")
            self.left_title = config.get("Plot", "left_title", fallback="View from above")
            self.right_title = config.get("Plot", "right_title", fallback="View from Earth")

            self.figure_width = eval(config.get("Plot", "figure_width", fallback="16"))
            self.figure_height = eval(config.get("Plot", "figure_height", fallback="8"))
            self.xlim = eval(config.get("Plot", "xlim", fallback="1.25"))
            self.ylim = eval(config.get("Plot", "ylim", fallback="1.0"))
            self.flux_dot_height = eval(config.get("Plot", "flux_dot_height", fallback="0.077"))
            self.flux_dot_width = eval(config.get("Plot", "flux_dot_width", fallback="0.005"))
            self.rv_dot_height = eval(config.get("Plot", "rv_dot_height", fallback="0.077"))
            self.rv_dot_width = eval(config.get("Plot", "rv_dot_width", fallback="0.005"))
        else:  # run MCMC, fit parameters to flux measurements
            # [Fitting]
            if self.tt_file:
                self.starts_d = np.array(eval(config.get("Simulation", "starts", fallback="[]")), dtype=float)
                self.ends_d = np.array(eval(config.get("Simulation", "ends", fallback="[]")), dtype=float)
                self.dts = np.array(eval(config.get("Simulation", "dts", fallback="[]")), dtype=float)
                # self.start_indices, self.max_iterations, self.total_iterations = self.check_intervals()
                self.best_residuals_tt_sum_squared = 1e99

            self.guifit = eval(config.get("Fitting", "guifit", fallback="False"))
            self.lmfit = eval(config.get("Fitting", "lmfit", fallback="False"))
            self.lmfit_method = config.get("Fitting", "lmfit_method", fallback="powell")
            self.lmfit_max_tt_delta = eval(config.get("Fitting", "lmfit_max_tt_delta", fallback="1/(24*60*60)"))
            self.flux_weight = int(eval(config.get("Fitting", "flux_weight", fallback="1")))
            self.tt_weight = int(eval(config.get("Fitting", "tt_weight", fallback="1")))

            self.backend = config.get("Fitting", "backend", fallback=None)  # e.g. emcee_backend.h5
            self.load_backend = eval(config.get("Fitting", "load_backend", fallback="False"))
            self.walkers = int(eval(config.get("Fitting", "walkers", fallback="32")))
            self.target_flux = eval(config.get("Fitting", "target_flux", fallback="None"))
            self.steps = int(eval(config.get("Fitting", "steps", fallback="10000")))
            self.moves = config.get("Fitting", "moves", fallback="None")
            self.burn_in = int(eval(config.get("Fitting", "burn_in", fallback="500")))
            if self.burn_in < 1:
                self.burn_in = 1
            self.chunk_size = int(eval(config.get("Fitting", "chunk_size", fallback="500")))
            self.bins = tuple([eval(x) for x in config.get("Fitting", "bins", fallback="30").split("#")[0].split(",")])
            self.thin_samples = int(eval(config.get("Fitting", "thin_samples", fallback="10")))

            default_unit = '{"mass": "m_jup", "radius": "r_jup", "e": "1", "i": "deg", "P": "d", "a": "AU", "Omega": "deg", "omega": "deg", "pomega": "deg", "L": "deg", "ma": "deg", "ea": "deg", "nu": "deg", "T": "s", "t": "s"}'
            dict_str = config.get("Fitting", "unit", fallback=default_unit)
            self.unit = eval(dict_str)
            default_scale = '{"mass": 1/m_jup, "radius": 1/r_jup, "e": 1, "i": rad2deg, "P": 1/day, "a": 1/au, "Omega": rad2deg, "omega": rad2deg, "pomega": rad2deg, "L": rad2deg, "ma": rad2deg, "ea": rad2deg, "nu": rad2deg, "T": 1, "t": 1}'
            dict_str = config.get("Fitting", "scale", fallback=default_scale)
            self.scale = eval(dict_str)
            self.fitting_parameters = self.read_fitting_parameters(config)

    def __repr__(self):
        return f"CurveSimParameters from {self.config_file}"

    def check_intervals(self):
        """Checks if the intervals in parameters starts_d, ends_d and dts are well defined.
           Calculates the indices for time_s0, time_d and sim_flux where the intervals start and end.
           Calculates the total number of iterations for which body positions and flux will be simulated and stored.
           Creates alternative parameters starts_s0, ends_s0 in seconds instead of days and starting with 0 at start_date.
         """
        if len(self.starts_d) == 0 or len(self.ends_d) == 0 or len(self.dts) == 0:
            print("At least one of the parameters starts/ends/dts is missing. Default values take effect.")
            self.starts_d = np.array([self.start_date], dtype=float)
            self.dts = np.array([self.dt], dtype=float)
            self.ends_d = np.array([self.start_date + (self.frames * self.fps * self.dt) / self.day], dtype=float)  # default value. Assumes the video shall last "frames" seconds.
        if not (len(self.starts_d) == len(self.ends_d) == len(self.dts)):
            print(f"{Fore.YELLOW}\nWARNING: Parameters starts, ends and dts do not have the same number of items.{Style.RESET_ALL}")
            print(f"{Fore.YELLOW}Only the first {min(len(self.starts_d), len(self.ends_d), len(self.dts))} intervals will be processed.{Style.RESET_ALL}")
        for start, end in zip(self.starts_d, self.ends_d):
            if start > end:
                print(f"{Fore.RED}\nERROR in parameters starts/ends: One interval ends before it begins.{Style.RESET_ALL}")
                sys.exit(1)
        for nextstart, end in zip(self.starts_d[1:], self.ends_d[:-1]):
            if end > nextstart:
                print(f"{Fore.RED}\nERROR in parameters starts/ends: One interval starts before its predecessor ends.{Style.RESET_ALL}")
                sys.exit(1)
        if self.start_date > self.starts_d[0]:
            print(f"{Fore.RED}\nERROR in parameter starts: First interval starts before the simulation's start_date.{Style.RESET_ALL}")
            sys.exit(1)
        self.starts_s0 = (self.starts_d - self.start_date) * self.day  # convert BJD to seconds and start at zero
        self.ends_s0 = (self.ends_d - self.start_date) * self.day  # convert BJD to seconds and start at zero
        max_iterations = [int((end - start) / dt) + 1 for start, end, dt in zip(self.starts_s0, self.ends_s0, self.dts)]  # each interval's number of iterations
        start_indices = [sum(max_iterations[:i]) for i in range(len(max_iterations) + 1)]  # indices of each interval's first iteration
        total_iterations = sum(max_iterations)
        return start_indices, max_iterations, total_iterations

    @staticmethod
    def find_and_check_config_file(config_file):
        """Check if config file can be opened and contains all standard sections."""
        # Check program parameters and extract config file name from them.
        # if len(sys.argv) == 1:
        #     config_file = default
        #     print(f"Using default config file {config_file}. Specify config file name as program parameter if you "
        #           f"want to use another config file.")
        # elif len(sys.argv) == 2:
        #     config_file = sys.argv[1]
        #     print(f"Using {config_file} as config file.")
        # else:
        #     config_file = sys.argv[1]
        #     print(f"Using {config_file} as config file. Further program parameters are ignored.")
        config = configparser.ConfigParser(inline_comment_prefixes="#")
        config.optionxform = str  # Preserve case of the keys.
        if len(config.read(config_file, encoding="utf-8")) < 1:  # does opening the config file fail?
            print(f"{Fore.RED}\nERROR: Config file {config_file} not found.{Style.RESET_ALL}")
            print(f"{Fore.RED}Provide the config file name as the argument of the function curvesim.{Style.RESET_ALL}")
            print(f"{Fore.RED}More information on https://github.com/lichtgestalter/curvesimulator/wiki {Style.RESET_ALL}")
            sys.exit(1)
        if not config_file.endswith(".ini"):
            print(f"{Fore.RED}Please only use config files with the .ini extension. (You tried to use {config_file}.){Style.RESET_ALL}")
            sys.exit(1)

        # for section in standard_sections:  # Does the config file contain all standard sections?
        #     if section not in config.sections() and section != "Debug":
        #         print(f"{Fore.RED}Section {section} missing in config file.{Style.RESET_ALL}")
        #         sys.exit(1)

    @staticmethod
    def init_time_arrays(p):
        time_s0 = np.zeros(p.max_iterations, dtype=float)
        # time_s0 = np.zeros(p.total_iterations, dtype=float)
        i = 0
        for start, dt, max_iteration in zip(p.starts_s0, p.dts, p.max_iterations):
            for j in range(max_iteration):
                time_s0[i] = start + j * dt
                i += 1
        time_d = time_s0 / p.day + p.start_date
        return time_s0, time_d

    def read_param(self, config, section, param, fallback):
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
        g, au, r_sun, m_sun, l_sun = self.g, self.au, self.r_sun, self.m_sun, self.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = self.r_jup, self.m_jup, self.r_earth, self.m_earth, self.v_earth
        hour, day, year = self.hour, self.day, self.year
        line = config.get(section, param, fallback=fallback)
        value = eval(line.split(",")[0])
        # read_param(config, section, "ma", fallback="None")
        if value is not None and param in ["i", "Omega", "omega", "pomega", "ma", "nu", "ea", "L"]:
            value = np.radians(value)
        return value

    def read_param_and_bounds(self, config, section, param):
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
        g, au, r_sun, m_sun, l_sun = self.g, self.au, self.r_sun, self.m_sun, self.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = self.r_jup, self.m_jup, self.r_earth, self.m_earth, self.v_earth
        hour, day, year = self.hour, self.day, self.year
        line = config.get(section, param, fallback=None)
        if line is None:
            return None, None, None, None
        else:
            items = line.split("#")[0].split(",")
        if len(items) == 4:
            value, lower, upper, sigma = items
            return eval(value), eval(lower), eval(upper), eval(sigma)
        else:
            return None, None, None, None

    def read_fitting_parameters(self, config):
        """Search for body parameters in the config file that are meant to
        be used as fitting parameters in MCMC.
        Fitting parameters have 4 values instead of 1, separated by commas:
        Initial Value, Lower Bound, Upper Bound, Standard Deviation of the Initial Values of all chains (with mean = Initial Value)."""
        body_index = 0
        fitting_parameters = []
        if self.verbose:
            print(f"Running MCMC with these fitting parameters:")
        for section in config.sections():
            if section not in self.standard_sections:  # section describes a physical object
                for parameter_name in ["mass", "radius", "e", "i", "a", "P", "Omega", "pomega", "omega", "L", "nu", "ma", "ea", "T"]:
                    value, lower, upper, sigma = self.read_param_and_bounds(config, section, parameter_name)
                    if value is not None:
                        if self.verbose:
                            print(f"body {body_index}: {parameter_name}")
                        if parameter_name in ["i", "Omega", "omega", "pomega", "ma", "nu", "ea", "L"]:
                            value, lower, upper, sigma = np.radians(value), np.radians(lower), np.radians(upper), np.radians(sigma)
                        fitting_parameters.append(FittingParameter(self, section, body_index, parameter_name, value, lower, upper, sigma))
                        fitting_parameters[-1].index = len(fitting_parameters) - 1
                body_index += 1
        # print(f"Fitting {len(fitting_parameters)} parameters.")
        return fitting_parameters

    @staticmethod
    def save_fitting_parameters(fitting_parameters, prefix="", suffix=""):
        """ fitting_parameters is as list of FittingParameter
        Save it in JSON format"""
        data = []
        for fp in fitting_parameters:
            item = {
                "index": fp.index,
                "body_name": fp.body_name,
                # hier (oder in enrich funktion?) lesbare (also mit scale multiplizierte) attribute erzeugen
                "body_index": fp.body_index,
                "parameter_name": fp.parameter_name,
                "unit": fp.unit,
                "long_parameter_name": fp.long_parameter_name,
                "scale": fp.scale,
                "startvalue": fp.startvalue,
                "lower": fp.lower,
                "upper": fp.upper,
                "sigma": fp.sigma,
            }
            data.append(item)

        base_name = f"{prefix}fitting_parameters{suffix}.json"
        filename = base_name
        counter = 1
        while os.path.exists(filename):
            filename = f"../fitting_parameters/{prefix}fitting_parameters{suffix}_{counter}.json"
            counter += 1

        payload = {"created": time.time(), "count": len(data), "fitting_parameters": data}
        with open(filename, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, ensure_ascii=False)

        print(f"CurveSimParameters.save_fitting_parameters: saved {len(data)} entries to {filename}")
        exit(88)

    def find_results_subdirectory(self):
        """Find the name of the non-existing subdirectory with
        the lowest number and create this subdirectory."""
        if not os.path.isdir(self.results_directory):
            print(f"{Fore.RED}\nERROR: Fitting results directory {self.results_directory} does not exist.{Style.RESET_ALL}")
            sys.exit(1)
        # Filter numeric subdirectory names and checks if they are directories.
        existing_subdirectories = [int(subdir) for subdir in os.listdir(self.results_directory)
                                   if subdir.isdigit() and os.path.isdir(os.path.join(self.results_directory, subdir))]
        next_subdirectory = 0
        while next_subdirectory in existing_subdirectories:
            next_subdirectory += 1
        self.results_directory = self.results_directory + f"/{next_subdirectory:04d}/"
        os.makedirs(self.results_directory)

    def bodynames2bodies(self, bodies):
        """ Generates 2 lists of bodies (self.eclipsers, self.eclipsees)
         based on 2 lists of strings (self.eclipsers_names, self.eclipsees_names)"""
        eclipsers, eclipsees = [], []
        for body in bodies:
            if body.name in self.eclipsers_names:
                eclipsers.append(body)
            if body.name in self.eclipsees_names:
                eclipsees.append(body)
        self.eclipsers, self.eclipsees = eclipsers, eclipsees

    def randomize_startvalues_uniform(self):
        for fp in self.fitting_parameters:
            fp.startvalue = fp.lower + random.random() * (fp.upper - fp.lower)

    def TOI4504_startvalue_hack(self):
        d_P = self.get_fitting_parameter(1, "P")
        c_P = self.get_fitting_parameter(2, "P")
        c_P.startvalue = (-2.3493 * d_P.startvalue * d_P.scale + 178.93) / c_P.scale

        d_o = self.get_fitting_parameter(1, "omega")
        d_o.startvalue = (-211.26 * d_P.startvalue * d_P.scale + 9104) / d_o.scale

        d_ma = self.get_fitting_parameter(1, "ma")
        d_ma.startvalue = (-2.1489 * d_o.startvalue * d_o.scale + 572.04) / d_ma.scale

        c_o = self.get_fitting_parameter(2, "omega")
        c_ma = self.get_fitting_parameter(2, "ma")
        c_ma.startvalue = (-1.0343 * c_o.startvalue * c_o.scale + 79.427) / c_ma.scale

        d_m = self.get_fitting_parameter(1, "mass")
        c_m = self.get_fitting_parameter(2, "mass")
        c_m.startvalue = (0.3517 * d_m.startvalue * d_m.scale + 1.8837) / c_m.scale

        # a_m = self.get_fitting_parameter(0, "mass")
        # d_m = self.get_fitting_parameter(1, "mass")
        # c_m = self.get_fitting_parameter(2, "mass")
        # d_m.startvalue = (-0.002505 * a_m.startvalue * a_m.scale + 0.131838) / d_m.scale
        # c_m.startvalue = (-0.002919 * a_m.startvalue * a_m.scale + 0.007444) / c_m.scale

    def get_fitting_parameter(self, body_index, parameter_name):
        return self.fitting_parameters[self.fitting_parameter_dic[(body_index, parameter_name)]]

    def init_fitting_parameter_dic(self):
        self.fitting_parameter_dic = {(fp.body_index, fp.parameter_name): fp.index for fp in self.fitting_parameters}

    def enrich_fitting_params(self, bodies):
        """Add attributes body_parameter_name and long_body_parameter_name to each FittingParameter"""
        self.body_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in self.fitting_parameters]
        self.long_body_parameter_names = [fpn + " [" + self.unit[fpn.split(".")[-1]] + "]" for fpn in self.body_parameter_names]
        for fp, fpn, fpnu in zip(self.fitting_parameters, self.body_parameter_names, self.long_body_parameter_names):
            fp.body_parameter_name = fpn
            fp.long_body_parameter_name = fpnu


class FittingParameter:
    def __init__(self, p, body_name, body_index, parameter_name, startvalue, lower, upper, sigma):
        self.body_name = body_name
        self.body_index = body_index
        self.parameter_name = parameter_name
        self.unit = p.unit[parameter_name]
        self.long_parameter_name = parameter_name + "[" + p.unit[parameter_name] + "]"
        self.scale = p.scale[parameter_name]
        self.startvalue = startvalue
        self.lower = lower
        self.upper = upper
        self.sigma = sigma

    def initial_values(self, rng, size):
        result = []
        while len(result) < size:
            sample = rng.normal(self.startvalue, self.sigma)
            if self.lower <= sample <= self.upper:
                result.append(sample)
        return np.array(result)

class CurveSimPhysics:

    @staticmethod
    def kepler_equation(ea, e, ma):
        """ea: eccentric anomaly [rad], e: eccentricity, ma: mean anomaly [rad]"""
        if not -2 * math.pi < ea < 2 * math.pi:
            raise ValueError("eccentric anomaly ea must be in radians but is outside of the range ]-2π;2π[")
        if not -2 * math.pi < ma < 2 * math.pi:
            raise ValueError("mean anomaly ma must be in radians but is outside of the range ]-2π;2π[")
        if not 0 <= e < 1:
            raise ValueError("eccentricity e is outside of the range [0;1[")
        return ea - e * math.sin(ea) - ma

    @staticmethod
    def kepler_equation_derivative(ea, e):
        """ea: eccentric anomaly [rad], e: eccentricity"""
        return 1.0 - e * math.cos(ea)

    @staticmethod
    def kepler_equation_root(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
        """Calculate the root of the Kepler Equation with the Newton–Raphson method.
            e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
        for n in range(max_steps):
            delta = CurveSimPhysics.kepler_equation(ea_guess, e, ma) / CurveSimPhysics.kepler_equation_derivative(ea_guess, e)
            if abs(delta) < tolerance:
                return ea_guess - delta
            ea_guess -= delta
        raise RuntimeError("Newton\'s root solver did not converge.")

    @staticmethod
    def gravitational_parameter(bodies, g):
        """Calculate the gravitational parameter of masses orbiting a common barycenter
        https://en.wikipedia.org/wiki/Standard_gravitational_parameter"""
        mass = 0.0
        for body in bodies:
            mass += body.mass
        # print(f"Gravitational parameter {g * mass:.3f}")
        return g * mass

    @staticmethod
    def distance_2d_body(body1, body2, i):
        """Return distance of the centers of 2 physical bodies as seen by a viewer (projection x->0)."""
        # dy = body1.positions[i][1] - body2.positions[i][1]
        # dz = body1.positions[i][2] - body2.positions[i][2]
        # return math.sqrt((dy ** 2 + dz ** 2))
        """Return distance of the centers of 2 physical bodies as seen by a viewer (projection z->0)."""
        dx = body1.positions[i][0] - body2.positions[i][0]
        dy = body1.positions[i][1] - body2.positions[i][1]
        return math.sqrt((dx ** 2 + dy ** 2))

    @staticmethod
    def distance_2d_particle(particle1, particle2):
        """Return distance of the centers of 2 rebound simulation particles (projection z->0)."""
        dx = particle1.x - particle2.x
        dy = particle1.y - particle2.y
        return math.sqrt((dx ** 2 + dy ** 2))

    @staticmethod
    def get_limbdarkening_parameters(parameter_1, parameter_2, parameter_type):
        """converts limb darkening parameters to quadratic law parameters u1,u2 if necessary"""
        if parameter_type == "u":
            return parameter_1, parameter_2
        if parameter_type == "a":
            u1 = parameter_1 + 2 * parameter_2
            u2 = -parameter_2
            return u1, u2
        if parameter_type == "q":
            u1 = 2 * math.sqrt(parameter_1) * parameter_2
            u2 = np.sqrt(parameter_1) * (1 - 2 * parameter_2)
            return u1, u2
        if parameter_type is None:
            return None, None
        print(f"{Fore.RED}\nERROR in config file: limb_darkening_parameter_type must be a or u or q.")
        print(f"                      limb_darkening must be [a0,a1,a2] or [u1,u2] or [q1,q2] correspondingly.")
        print(f"                      But config file contains: limb_darkening_parameter_type = {parameter_type}{Style.RESET_ALL}")
        sys.exit(1)

    @staticmethod
    def intensity(mu, u1, u2):
        """Apply quadratic limb darkening law"""
        return 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2

    @staticmethod
    def calc_mean_intensity(limb_darkening_u1, limb_darkening_u2):
        """Calculates the ratio of the mean intensity to the central intensity of a star based on
        the given quadratic law parameters for limb darkening by integrating the intensity over the stellar disk"""
        if limb_darkening_u1 is None or limb_darkening_u2 is None:
            return None
        mu_values = np.linspace(0, 1, 1000)
        intensities = CurveSimPhysics.intensity(mu_values, limb_darkening_u1, limb_darkening_u2)
        return 2 * np.trapz(intensities * mu_values, mu_values)

    @staticmethod
    def limbdarkening(relative_radius, limb_darkening_u1, limb_darkening_u2):
        """
        Approximates the flux of a star at a point on the star seen from a very large distance.
        The point's apparent distance from the star's center is relative_radius * radius.

        Parameters:
        relative_radius (float): The normalized radial coordinate (0 <= x <= 1).
        limb_darkening_parameters: list of coefficients for the limb darkening model.

        Returns:
        float: intensity relative to the intensity at the midlle of the star at the given relative radius.
        """
        if relative_radius < 0:  # handling rounding errors
            relative_radius = 0.0
        if relative_radius > 1:
            relative_radius = 1.0
        mu = math.sqrt(1 - relative_radius ** 2)
        return CurveSimPhysics.intensity(mu, limb_darkening_u1, limb_darkening_u2)

    @staticmethod
    def distance_3d(point1, point2):
        x1, y1, z1 = point1
        x2, y2, z2 = point2
        return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

    @staticmethod
    def calc_transit_intervals(t1, t2, t3, t4):
        t12 = None if t2 is None or t1 is None else t2 - t1
        t23 = None if t3 is None or t2 is None else t3 - t2
        t34 = None if t4 is None or t3 is None else t4 - t3
        t14 = None if t4 is None or t1 is None else t4 - t1
        return t12, t23, t34, t14


class CurveSimRebound:
    def __init__(self, rebound_sim):
        self.sim = rebound_sim
        self.energy = rebound_sim.energy()
        self.total_momentum = self.calc_total_momentum()
        self.total_angular_momentum = self.calc_total_angular_momentum()
        self.center_of_mass_position = self.calc_center_of_mass_position()

    def sim_check_deltas(self, newer):
        energy_change = (newer.energy - self.energy) / self.energy
        # print(f"Energy at Start: {self.energy=:.1e}   End: {newer.energy=:.1e}")
        # to avoid div/0 i calculate the norm _before_ dividing by self.total_momentum
        # something is wrong with self.calc_total_momentum(). Therefore I currently ignore it.
        # total_momentum_change = (np.linalg.norm(newer.total_momentum) - np.linalg.norm(self.total_momentum)) / np.linalg.norm(self.total_momentum)
        # total_angular_momentum_change = np.linalg.norm((newer.total_angular_momentum - self.total_angular_momentum) / self.total_angular_momentum)
        # center_of_mass_position_change = np.linalg.norm(newer.center_of_mass_position - self.center_of_mass_position)
        # print(f"{energy_change=:.2e}  {total_momentum_change=:.2e}  {total_angular_momentum_change=:.2e}  {center_of_mass_position_change=:.2e}  ")
        if abs(energy_change) < 1e-20:
            return -20
        else:
            return math.log10(abs(energy_change))


    def calc_total_momentum(self):
        """
        Returns the total linear momentum vector of the system.
        """
        momentum_x, momentum_y, momentum_z = 0.0, 0.0, 0.0
        for p in self.sim.particles:
            momentum_x += p.m * p.vx
            momentum_y += p.m * p.vy
            momentum_z += p.m * p.vz
        return np.array([momentum_x, momentum_y, momentum_z], dtype=float)

    def calc_total_angular_momentum(self):
        """
        Returns the total angular momentum vector of the system.
        """
        Lx, Ly, Lz = 0.0, 0.0, 0.0
        for p in self.sim.particles:
            # Position vector
            r = np.array([p.x, p.y, p.z], dtype=float)
            # Momentum vector
            v = np.array([p.vx, p.vy, p.vz], dtype=float)
            p_vec = p.m * v
            # Angular momentum: L = r x p
            L = np.cross(r, p_vec)
            Lx += L[0]
            Ly += L[1]
            Lz += L[2]
        return np.array([Lx, Ly, Lz], dtype=float)

    def calc_center_of_mass_position(self):
        """
        Returns the position vector of the center of mass.
        """
        mass = 0.0
        x, y, z = 0.0, 0.0, 0.0
        for p in self.sim.particles:
            mass += p.m
            x += p.m * p.x
            y += p.m * p.y
            z += p.m * p.z
        if mass == 0.0:
            return 0
        return np.array([x/mass, y/mass, z/mass], dtype=float)
# from colorama import Fore, Style


class Transit(dict):
    def __init__(self, eclipsed_body):
        super().__init__()
        self["Transit_params"] = {}
        transit_params = ["EclipsedBody", "T1", "T2", "TT", "T3", "T4", "T12", "T23", "T34", "T14", "b"]
        for key in transit_params:
            self["Transit_params"][key] = None
        self["Transit_params"]["EclipsedBody"] = eclipsed_body.name


class CurveSimResults(dict):
    def __init__(self, bodies, empty=False):
        super().__init__()
        if empty:
            return
        self["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        self["ProgramParameters"] = {}
        self["Bodies"] = {}
        exclude = ["positions", "velocity", "circle_left", "circle_right", "acceleration", "d", "h", "angle", "eclipsed_area", "patch_radius"]
        for body in bodies:
            self["Bodies"][body.name] = {}
            self["Bodies"][body.name]["BodyParameters"] = {}
            self["Bodies"][body.name]["Transits"]= []
            for key in body.__dict__.keys():
                if key not in exclude:
                    self["Bodies"][body.name]["BodyParameters"][key] = getattr(body, key)
            if body.Omega is not None:
                self["Bodies"][body.name]["BodyParameters"]["Omega_deg"] = body.Omega * (180 / math.pi)
            if body.omega is not None:
                self["Bodies"][body.name]["BodyParameters"]["omega_deg"] = body.omega * (180 / math.pi)
            if body.pomega is not None:
                self["Bodies"][body.name]["BodyParameters"]["pomega_deg"] = body.pomega * (180 / math.pi)
            if body.L is not None:
                self["Bodies"][body.name]["BodyParameters"]["L_deg"] = body.L * (180 / math.pi)
            if body.ma is not None:
                self["Bodies"][body.name]["BodyParameters"]["ma_deg"] = body.ma * (180 / math.pi)
            if body.ea is not None:
                self["Bodies"][body.name]["BodyParameters"]["ea_deg"] = body.ea * (180 / math.pi)
            if body.nu is not None:
                self["Bodies"][body.name]["BodyParameters"]["nu_deg"] = body.nu * (180 / math.pi)
            # for key in list(body.__dict__.keys()):  # uncomment to prevent null-values in result file
            #     if body.__dict__[key] is None:
            #         del body.__dict__[key]


    def __repr__(self):
        string = ""
        for body in self["Bodies"]:
            if len(self["Bodies"][body]["Transits"]) == 1:
                string += f"{body:15} {len(self["Bodies"][body]["Transits"]):3} transit\n"
            elif len(self["Bodies"][body]["Transits"]) > 1:
                string += f"{body:15} {len(self["Bodies"][body]["Transits"]):3} transits\n"
        return string[:-1]

    @staticmethod
    def iteration2time(iteration, p):
        """Calculate the date of an iteration in BJD"""
        return p.start_date + iteration * p.dt / p.day

    def results2json(self, p):
        """Converts self to JSON and saves it."""
        filename = p.results_directory + p.result_file
        with open(filename, "w", encoding="utf8") as file:
            json.dump(self, file, indent=4, ensure_ascii=False)
        if p.verbose:
            print(filename, "saved")

    @staticmethod
    def check_resultfilename(resultfilename):
        """Check if resultfilename already exists and attach a number if it does."""
        if not os.path.exists(resultfilename):
            return resultfilename
        base, ext = os.path.splitext(resultfilename)
        match = re.search(r"\.v(\d+)$", base)
        if match:
            num = int(match.group(1)) + 1
            base = base[:match.start()]
        else:
            num = 1
        new_resultfilename = f"{base}.v{num:04}{ext}"
        while os.path.exists(new_resultfilename):
            num += 1
            new_resultfilename = f"{base}.v{num:04}{ext}"
        return new_resultfilename

    def save_results(self, p):
        p_copy = copy.deepcopy(p)
        to_remove = [
            "fitting_parameters", "standard_sections", "eclipsers", "eclipsees",
            "fitting_parameter_dic",
        ]
        for name in to_remove:
            if hasattr(p_copy, name):
                delattr(p_copy, name)
        for name in ("starts_s0", "starts_d", "ends_s0", "ends_d", "dts"):
            if hasattr(p_copy, name):
                orig = getattr(p_copy, name)
                p_copy.__dict__[name] = [float(i) for i in orig]
        self["ProgramParameters"] = p_copy.__dict__


        # diagnostic helper
        def _find_unserializable(obj, path="self", visited=None, max_depth=1000):
            if visited is None:
                visited = set()
            obj_id = id(obj)
            if obj_id in visited or max_depth <= 0:
                return []
            visited.add(obj_id)
            try:
                json.dumps(obj)
                return []
            except Exception:
                # drill down for containers
                failures = []
                if isinstance(obj, dict):
                    for k, v in obj.items():
                        failures.extend(_find_unserializable(v, f"{path}[{repr(k)}]", visited, max_depth - 1))
                    if not failures:
                        failures.append((path, type(obj).__name__, "dict contains non-serializable contents"))
                elif isinstance(obj, (list, tuple, set)):
                    for i, v in enumerate(obj):
                        failures.extend(_find_unserializable(v, f"{path}[{i}]", visited, max_depth - 1))
                    if not failures:
                        failures.append((path, type(obj).__name__, f"{type(obj).__name__} contains non-serializable contents"))
                else:
                    # leaf non-serializable object
                    failures.append((path, type(obj).__name__, repr(obj)[:200]))
                return failures

        # check serialization of self
        failures = _find_unserializable(self)
        if failures:
            # format a concise message with a few examples
            msg_lines = ["Failed to JSON-serialize `self`. Problematic paths (path, type, sample):"]
            for path, tname, sample in failures[:20]:
                msg_lines.append(f" - {path}: {tname} -> {sample}")
            raise RuntimeError("\n".join(msg_lines))

        self.results2json(p_copy)
        if p.verbose:
            print(self)

    @staticmethod
    def load_results(resultfilename):
        """Load JSON results from `resultfilename` into and return a CurveSimResults object."""
        if not os.path.exists(resultfilename):
            raise FileNotFoundError(resultfilename)
        try:
            with open(resultfilename, "r", encoding="utf8") as f:
                data = json.load(f)
        except Exception as e:
            raise RuntimeError(f"Failed to read or parse {resultfilename}: {e}")
        if not isinstance(data, dict):
            raise ValueError(f"Invalid results file {resultfilename}: top-level JSON must be an object")
        results = CurveSimResults(None, empty=True)
        results.clear()
        results.update(data)
        return results

    @staticmethod
    def calc_rv_residuals(measured_rv, body_name, rebound_sim):

        def rv_at_t(t, rebound_sim, body):
            rebound_sim.integrate(t)
            return -body.vz

        body = rebound_sim.particles[body_name]
        measured_rv["rv_sim"] = [rv_at_t(t, rebound_sim, body) for t in measured_rv["time_s0"]]
        measured_rv["residual"] = measured_rv["rv_rel"] - measured_rv["rv_sim"]
        return measured_rv

    @staticmethod
    def calc_flux_residuals(measured_flux, sim_flux):
        measured_flux["flux_sim"] = sim_flux
        measured_flux["residual"] = measured_flux["flux"] - measured_flux["flux_sim"]
        return measured_flux

    def calc_rv_chi_squared(self, measured_rv, free_parameters):
        measured_rv["chi_squared"] = measured_rv["residual"] / measured_rv["rv_jit"]
        measured_rv["chi_squared"] = measured_rv["chi_squared"] * measured_rv["chi_squared"]
        self["chi_squared_rv"] = measured_rv["chi_squared"].sum()
        self["measurements_rv"] = measured_rv.shape[0]
        self["pvalue_rv"] = CurveSimResults.chi_squared_pvalue(self["chi_squared_rv"], self["measurements_rv"], free_parameters)
        return measured_rv

    def calc_flux_chi_squared(self, measured_flux, free_parameters):
        measured_flux["chi_squared"] = measured_flux["residual"] / measured_flux["flux_err"]
        measured_flux["chi_squared"] = measured_flux["chi_squared"] * measured_flux["chi_squared"]
        self["chi_squared_flux"] = measured_flux["chi_squared"].sum()
        self["measurements_flux"] = measured_flux.shape[0]
        self["pvalue_flux"] = CurveSimResults.chi_squared_pvalue(self["chi_squared_flux"], self["measurements_flux"], free_parameters)
        return measured_flux

    def calc_tt_chi_squared(self, measured_tt, free_parameters):
        measured_tt["chi_squared"] = measured_tt["delta"] / measured_tt["tt_err"]
        measured_tt["chi_squared"] = measured_tt["chi_squared"] * measured_tt["chi_squared"]
        self["chi_squared_tt"] = measured_tt["chi_squared"].sum()
        self["measurements_tt"] = measured_tt.shape[0]
        self["pvalue_tt"] = CurveSimResults.chi_squared_pvalue(self["chi_squared_tt"], self["measurements_tt"], free_parameters)
        return measured_tt

    def calc_total_chi_squared(self, free_parameters):
        self["chi_squared_total"] = self["chi_squared_rv"] + self["chi_squared_flux"] + self["chi_squared_tt"]
        self["measurements_total"] = self["measurements_rv"] + self["measurements_flux"] + self["measurements_tt"]
        self["pvalue_total"] = CurveSimResults.chi_squared_pvalue(self["chi_squared_total"], self["measurements_total"], free_parameters)


    @staticmethod
    def chi_squared_pvalue(chi_squared, n_measurements, n_parameters):
        """
        Calculate the p-value for a chi-squared test.
        This is the probability of observing a chi-squared value >= your observed value.
        chi_square :     The chi-squared test statistic
        n_measurements : Number of measurements/observations
        n_parameters :   Number of free parameters in the model
        """
        df = n_measurements - n_parameters  # degrees of freedom
        p_value = 1 - stats.chi2.cdf(chi_squared, df)  # cumulative distribution function
        # print(f"Chi-squared value:      {chi_squared}")
        # print(f"Number of measurements: {n_measurements}")
        # print(f"Number of parameters:   {n_parameters}")
        # print(f"Degrees of freedom:     {df}")
        # print(f"P-value:                {p_value:.4f}\n")
        # if p_value > 0.05:
        #     print(f"The fit is acceptable (p > 0.05). There is no significant evidence that your model is inconsistent with the data.")
        # else:
        #     print(f"The fit is poor (p < 0.05). Your model may not adequately describe the data.")
        return p_value

    @staticmethod
    def get_measured_flux(p):
        measured_flux = pd.read_csv(p.flux_file)
        measured_flux = measured_flux[measured_flux["time"] >= p.start_date]
        measured_flux["time_s0"] = (measured_flux["time"] - p.start_date) * p.day
        time_s0 = np.array(measured_flux["time_s0"], dtype=float)
        measured_flux_array = np.array(measured_flux["flux"])
        flux_err = np.array(measured_flux["flux_err"], dtype=float)
        p.total_iterations = len(time_s0)
        time_d = time_s0 / p.day + p.start_date
        return time_s0, time_d, measured_flux_array, flux_err, measured_flux

    @staticmethod
    def get_measured_tt(p):
        df = pd.read_csv(p.tt_file)
        df = df[df["tt"] >= p.start_date]
        p.tt_datasize = len(df["tt"])
        return df

    @staticmethod
    def get_measured_rv(p):
        df = pd.read_csv(p.rv_file)
        df = df[df["time"] >= p.start_date]
        p.rv_datasize = len(df["time"])
        df["time_s0"] = (df["time"] - p.start_date) * p.day
        return df

    @staticmethod
    def _to_list(val, default):
        if val is None:
            return [default]
        if isinstance(val, (list, tuple)):
            return list(val)
        return [val]

    @staticmethod
    def _expand_param(lst, name, n):
        if len(lst) == 1:
            return lst * n
        if len(lst) == n:
            return lst
        raise ValueError(f"{name} must have length 1 or {n}")

    @staticmethod
    def plot_this(
            x_lists,                   # positions of data points on x-axis
            y_lists: list,             # each list item is a list or numpy array which will be displayed as a curve
            data_labels: list = None,  # each list item is a string representing the label of a curve
            title: str = None,         # plot title
            x_label: str = None,       # label of x-axis
            y_label: str = None,       # label of y-axis
            markers=None,              # single marker or list/tuple of markers
            markersizes=None,          # single markersize or list/tuple of markersizes
            linestyles=None,           # single linestyle or list/tuple of linestyles
            colors=None,               # single color or list/tuple of colors
            linewidths=None,           # single linewidt or list/tuple of linewidts
            left: float = None,        # cut off x-axis
            right: float = None,       # cut off x-axis
            bottom: float = None,      # cut off y-axis
            top: float = None,         # cut off y-axis
            legend: bool = None,       # display legend?
            grid: bool = None,         # display grid?
            show_plot: bool = False,   # show plot?
            plot_file: str = None,     # file name if the plot shall be saved as .png
    ) -> None:
        if data_labels is None:
            data_labels = [f"data{i}" for i in range(len(y_lists))]
        plt.figure(figsize=(10, 6))
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.ticklabel_format(useOffset=False, style="plain", axis="x")   # show x-labels as they are

        n = len(y_lists)
        markers_list = CurveSimResults._to_list(markers, "o")
        markersizes_list = CurveSimResults._to_list(markersizes, 1)
        linestyles_list = CurveSimResults._to_list(linestyles, "None")
        colors_list = CurveSimResults._to_list(colors, "xkcd:tomato")
        linewidths_list = CurveSimResults._to_list(linewidths, 1.5)

        markers_per_curve = CurveSimResults._expand_param(markers_list, "markers", n)
        markersizes_per_curve = CurveSimResults._expand_param(markersizes_list, "markersizes", n)
        linestyles_per_curve = CurveSimResults._expand_param(linestyles_list, "linestyles", n)
        colors_per_curve = CurveSimResults._expand_param(colors_list, "colors", n)
        linewidths_per_curve = CurveSimResults._expand_param(linewidths_list, "linewidths", n)

        # Determine x arrays for each y-list:
        # If x_lists is a list/tuple with same length as y_lists and length > 1, treat as per-curve x arrays.
        # Otherwise treat x_lists as a single x array and broadcast to all curves.
        if isinstance(x_lists, (list, tuple)) and ((len(x_lists) == len(y_lists) and len(x_lists) > 1) or len(x_lists) == 1):
            x_iter = x_lists
        else:
            # single x array (can be ndarray or list); convert to numpy array for plotting
            single_x = np.array(x_lists)
            x_iter = [single_x] * len(y_lists)

        if left or right:
            plt.xlim(left=left, right=right)
        if bottom or top:
            plt.ylim(bottom=bottom, top=top)

        for x, data, data_label, marker, msize, ls, col, lw in zip(
                x_iter, y_lists, data_labels, markers_per_curve, markersizes_per_curve, linestyles_per_curve, colors_per_curve, linewidths_per_curve):
            plt.plot(x, data, marker=marker, markersize=msize, linestyle=ls, label=data_label, color=col, linewidth=lw)
        # for x, data, data_label in zip(x_iter, y_lists, data_labels):
        #     plt.plot(x, data, marker=marker, markersize=markersize, linestyle=linestyle, label=data_label)
        # for data, data_label in zip(y_lists, data_labels):
        #     plt.plot(x, data, marker=marker, markersize=markersize, linestyle=linestyle, label=data_label)

        if legend:
            plt.legend()
        if grid:
            plt.grid(True)
        if plot_file:
            plt.savefig(plot_file)
        if show_plot:
            plt.show()
        plt.close()


    def plot_parameter(self, eclipser, eclipsee, parameter, start, end, filename="", title=None):
        if not title:
            title = f"{os.path.splitext(os.path.basename(filename))[0]}, {self["ProgramParameters"]["comment"]} (dt={self["ProgramParameters"]["dt"]})"
        transit_times = self.get_transit_data(eclipser, eclipsee, "TT")
        parameter_list = self.get_transit_data(eclipser, eclipsee, parameter)
        CurveSimResults.plot_this(
            x_lists=transit_times,
            y_lists=[parameter_list],
            data_labels=[eclipser],
            title=title,
            x_label="Transit Times [BJD]",
            y_label=parameter,
            linestyles="-",
            markersizes=4,
            grid=True,
            left=start,
            right=end,
            plot_file=filename,
        )

    def get_transit_data(self, eclipser, eclipsee, parameter):
        transits = self["Bodies"][eclipser]["Transits"]
        parameter_list = [transit["Transit_params"][parameter] for transit in transits if transit["Transit_params"]["EclipsedBody"] == eclipsee]
        return parameter_list

    @staticmethod
    def ttv_to_date_plot(p, amplitude, period, x_offset, osc_per):
        measured_tt = CurveSimResults.get_measured_tt(p)  # reads from p.tt_file
        tt_tess = np.array(measured_tt["tt"][:13], dtype=float)
        transit_numbers = np.array(measured_tt["nr"][:13], dtype=float)
        ttv_to_date = [tt - tt_tess[0] - n * osc_per for n, tt in zip(transit_numbers, tt_tess)]

        x4sine = np.linspace(tt_tess[0], tt_tess[-1], 300)
        sine_curve = amplitude * np.sin(2 * np.pi * (x4sine - tt_tess[0] - x_offset) / period)
        sine_curve -= sine_curve[0]

        CurveSimResults.plot_this(
            title=f"TESS TT vs. mean osculating Period of TOI-4504 c \nOsc.per.={osc_per:.2f} d, Amplitude={amplitude:.2f} d, x-Offset={x_offset:.2f} d, Super Period={period:.2f})",
            x_label="Transit Times [BJD]",
            y_label="TTV to date [days]",
            x_lists=    [tt_tess,       x4sine],
            y_lists=    [ttv_to_date,   sine_curve],
            data_labels=["TTV to date", "Sine Curve"],
            linestyles= ["",            "-"],
            markersizes=[4,             0],
            colors=     ["xkcd:tomato",         "xkcd:black"],
            linewidths= [0,             1],
            grid=True,
            legend=True,
            plot_file=f"TTV_to_date_{osc_per:.4f}.png",
        )

    @staticmethod
    def sim_rv_plot(p, sim_rv, time_d, plot_filename):
        CurveSimResults.plot_this(
            title=f"Simulated Radial Velocity",
            x_label="Time [BJD]",
            y_label="RV [m/s]",
            x_lists=    [time_d],
            y_lists=    [sim_rv],
            data_labels=["sim_rv"],
            linestyles= ["-"],
            markersizes=[0],
            colors=     ["xkcd:black"],
            linewidths= [1],
            grid=False,
            legend=False,
            plot_file=p.results_directory + plot_filename,
        )

    @staticmethod
    def rv_observed_computed_plot(p, sim_rv, time_d, plot_filename, measured_rv):
        CurveSimResults.plot_this(
            title=f"Radial Velocity: observed vs. computed",
            x_label="Time [BJD]",
            y_label="RV [m/s]",
            x_lists=    [time_d,   measured_rv["time"]],
            y_lists=    [sim_rv,   measured_rv["rv_rel"]],
            data_labels=["computed", "observed"],
            linestyles= ["-",      ""],
            markersizes=[0,        3],
            colors=     ["xkcd:black",  "xkcd:nice blue"],
            linewidths= [1,        0],
            grid=False,
            legend=True,
            left=np.min(measured_rv["time"]) - 5,  # debug: offset in Parameterfile aufnehmen?
            right=np.max(measured_rv["time"]) + 5,
            plot_file=p.results_directory + plot_filename,
        )

    @staticmethod
    def bin_time_window(time, value, half_window_size):
        """ time, value: pandas DataSeries or numpy array
            half_window_size: float
            returns array with binned values
            Bins for each data point with time=t all values with time between
            t - half_window_size and t + half_window_size"""
        time = time.to_numpy()
        value = value.to_numpy()
        mean = np.empty(len(time), dtype=float)
        for i, t in enumerate(time):
            left = bisect.bisect_right(time, t - half_window_size)
            right = bisect.bisect_left(time, t + half_window_size)
            values_to_bin = value[left:right]
            mean[i] = values_to_bin.mean() if values_to_bin.size > 0 else np.nan
        return mean


    @staticmethod
    def bin_time_window2(time, value, half_window_size):
        """
        Vectorized binning: for each time[i], compute mean of values within
        [time[i] - half_window_size, time[i] + half_window_size].

        Returns a numpy array of length len(time) with NaN where no points fall in window.

        Ensure `half_window_size` is in the same units as `time` (e.g. days).
        """
        time_arr = np.asarray(time)
        val_arr = np.asarray(value, dtype=float)

        if time_arr.size == 0:
            return np.array([], dtype=float)

        # sort by time for fast search
        order = np.argsort(time_arr)
        sorted_time = time_arr[order]
        sorted_val = val_arr[order]

        # mask NaNs in values: contribute 0 to sum and 0 to count
        valid_mask = ~np.isnan(sorted_val)
        vals_for_sum = np.where(valid_mask, sorted_val, 0.0)

        # cumulative sums for sums and counts (prepend zero for easy range subtraction)
        cumsum_vals = np.concatenate(([0.0], np.cumsum(vals_for_sum)))
        cumsum_counts = np.concatenate(([0], np.cumsum(valid_mask.astype(int))))

        # compute left/right indices for each original time (use original times so result keeps input order)
        left = np.searchsorted(sorted_time, time_arr - half_window_size, side="left")
        right = np.searchsorted(sorted_time, time_arr + half_window_size, side="right")

        # window sums and counts
        window_sums = cumsum_vals[right] - cumsum_vals[left]
        window_counts = cumsum_counts[right] - cumsum_counts[left]

        # compute means, set NaN where count == 0
        means = np.full(len(time_arr), np.nan, dtype=float)
        nonzero = window_counts > 0
        means[nonzero] = window_sums[nonzero] / window_counts[nonzero]

        return means


    @staticmethod
    def flux_observed_computed_plots_time(p, plot_filename, measured_flux, measured_tt):
        measured_flux["gt"] = measured_flux["time"].diff().gt(0).astype(int)
        left = np.min(measured_flux["time"])
        right = np.max(measured_flux["time"])
        left -= (right - left) * 0.02
        right += (right - left) * 0.02
        measured_flux["bin_30min"] = CurveSimResults.bin_time_window2(measured_flux["time"], measured_flux["flux"], 30 / (2 * 60 * 24))
        CurveSimResults.plot_this(
            title=f"Flux: observed vs. computed",
            x_label="Time [BJD]",
            y_label="Normalized Flux",
            x_lists=    [measured_flux["time"], measured_flux["time"]],
            y_lists=    [measured_flux["flux"], measured_flux["flux_sim"]],
            data_labels=["observed",            "computed"],
            linestyles= ["",                    ""],
            markersizes=[1,                     1],
            colors=     ["xkcd:nice blue",                 "xkcd:black"],
            # linewidths= [1,                     0],
            grid=False,
            legend=True,
            left=left,
            right=right,
            top=p.flux_plots_top,
            bottom=p.flux_plots_bottom,
            plot_file=p.results_directory + plot_filename,
        )
        if measured_tt is not None:
            directory = p.results_directory + "flux_per_transit/"
            os.makedirs(directory)
            for transit in measured_tt.itertuples(index=False):
                CurveSimResults.plot_this(
                    title=f"{transit.eclipser} Transit nr. {transit.nr}: observed vs. computed flux",
                    x_label="Time [BJD]",
                    y_label="Normalized Flux",
                    x_lists=    [measured_flux["time"], measured_flux["time"],      measured_flux["time"]],
                    y_lists=    [measured_flux["flux"], measured_flux["bin_30min"], measured_flux["flux_sim"]],
                    data_labels=["observed",            "obs. 30 min",              "computed"],
                    linestyles= ["",                    "-",                        ""],
                    markersizes=[1,                     0,                          1],
                    colors=     ["xkcd:nice blue",      "xkcd:nice blue",           "xkcd:black"],
                    # linewidths= [1,                     0],
                    grid=False,
                    legend=True,
                    left=transit.tt - p.tt_padding,  # debug: offset in Parameterfile aufnehmen?
                    right=transit.tt + p.tt_padding,
                    top=p.flux_plots_top,
                    bottom=p.flux_plots_bottom,
                    plot_file=f"{directory}{transit.eclipser}_{transit.nr}_{plot_filename}",
                )

    @staticmethod
    def flux_observed_computed_plot_data(p, plot_filename, measured_flux):
        CurveSimResults.plot_this(
            title=f"Flux: observed vs. computed",
            x_label="Datapoints",
            y_label="Normalized Flux",
            x_lists=    [[x for x in range(measured_flux.shape[0])],     [x for x in range(measured_flux.shape[0])]],
            y_lists=    [measured_flux["flux"],                          measured_flux["flux_sim"]],
            data_labels=["observed",                                     "computed"],
            linestyles= ["",                                             ""],
            markersizes=[1,                                              1],
            colors=     ["xkcd:nice blue",                                          "xkcd:black"],
            # linewidths= [1,                                              0],
            grid=False,
            legend=True,
            plot_file=p.results_directory + plot_filename,
        )

    @staticmethod
    def flux_chi_squared_plot_data(p, plot_filename, measured_flux):
        CurveSimResults.plot_this(
            title=f"Flux: Chi Squared per data point",
            x_label="Datapoints",
            y_label="Chi Squared",
            x_lists=    [[x for x in range(measured_flux.shape[0])]],
            y_lists=    [measured_flux["chi_squared"]],
            data_labels=["chi_squared"],
            linestyles= [""],
            markersizes=[1],
            colors=     ["xkcd:tree green"],
            linewidths= [1],
            grid=False,
            legend=True,
            plot_file=p.results_directory + plot_filename,
        )

    @staticmethod
    def rv_residuals_plot(p, plot_filename, measured_rv):
        title = f"Radial Velocity: Residuals (observed minus computed)"
        x_label = "Time [BJD]"
        y_label = "RV [m/s]"
        x = [measured_rv["time"]]
        y = [measured_rv["residual"]]
        data_labels = ["residual"]
        linestyles = [""]
        markers = ["o"]
        markersizes = [4]
        colors = ["xkcd:nice blue"]
        linewidths = [0]
        xpaddings = [0.01]
        # xpaddings = [0.01 * (np.max(x) - np.min(x))]
        left = np.min(x) - xpaddings[0] * (np.max(x[0]) - np.min(x[0]))
        right = np.max(x) + xpaddings[0] * (np.max(x[0]) - np.min(x[0]))
        # bottom = None
        # top = None
        plot_file = p.results_directory + plot_filename

        plt.figure(figsize=(10, 6))
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        # plt.legend()
        # plt.grid(True)
        plt.ticklabel_format(useOffset=False, style="plain", axis="x")   # show x-labels as they are
        plt.xlim(left=left, right=right)
        # plt.ylim(bottom=bottom, top=top)

        for time, residual, jitter in zip(x, y, measured_rv["rv_jit"]):
            plt.vlines(time, residual - jitter, residual + jitter, colors="xkcd:black", linewidth=1)

        plt.plot(x[0], y[0], marker=markers[0], markersize=markersizes[0], linestyle=linestyles[0], label=data_labels[0], color=colors[0], linewidth=linewidths[0])
        plt.hlines(0, left, right, colors="xkcd:black", linewidth=1)
        plt.hlines(measured_rv["rv_jit"].mean(), left, right, colors="xkcd:warm grey", linewidth=1, linestyles="--")
        plt.hlines(-measured_rv["rv_jit"].mean(), left, right, colors="xkcd:warm grey", linewidth=1, linestyles="--")
        plt.savefig(plot_file)

    @staticmethod
    def flux_residuals_plots_time(p, plot_filename, measured_flux, measured_tt):
        CurveSimResults.flux_residuals_plot_time(p, plot_filename, measured_flux)
        if measured_tt is not None:
            subdirectory = "fluxresiduals_per_transit/"
            os.makedirs(p.results_directory + subdirectory)
            for transit in measured_tt.itertuples(index=False):
                title = f"{transit.eclipser} Transit nr. {transit.nr}: flux residuals"
                left = transit.tt - p.tt_padding
                right = transit.tt + p.tt_padding
                plot_file = f"{transit.eclipser}_{transit.nr}_{plot_filename}"
                CurveSimResults.flux_residuals_plot_time(p, plot_file, measured_flux, left=left, right=right, title=title, subdirectory=subdirectory)

    @staticmethod
    def flux_residuals_plot_time(p, plot_filename, measured_flux, left=None, right=None, title=None, subdirectory=""):
        if title is None:
            title = f"Flux Residuals (observed minus computed)"
        x_label = "Time [BJD]"
        y_label = "Normalized Flux"
        x = [measured_flux["time"]]
        y = [measured_flux["residual"]]
        data_labels = ["residual"]
        linestyles = [""]
        markers = ["o"]
        markersizes = [3]
        colors = ["xkcd:nice blue"]
        linewidths = [0]
        xpaddings = [0.01]
        # xpaddings = [0.01 * (np.max(x) - np.min(x))]
        if left is None:
            left = np.min(x) - xpaddings[0] * (np.max(x[0]) - np.min(x[0]))
        if right is None:
            right = np.max(x) + xpaddings[0] * (np.max(x[0]) - np.min(x[0]))
        # bottom = None
        # top = None
        plot_file = p.results_directory + subdirectory + plot_filename

        plt.figure(figsize=(10, 6))
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        # plt.legend()
        # plt.grid(True)
        plt.ticklabel_format(useOffset=False, style="plain", axis="x")   # show x-labels as they are
        plt.xlim(left=left, right=right)
        # plt.ylim(bottom=bottom, top=top)
        plt.vlines(measured_flux["time"], measured_flux["residual"] - measured_flux["flux_err"], measured_flux["residual"] + measured_flux["flux_err"], colors="xkcd:black", linewidth=1)
        plt.plot(x[0], y[0], marker=markers[0], markersize=markersizes[0], linestyle=linestyles[0], label=data_labels[0], color=colors[0], linewidth=linewidths[0])
        plt.hlines(0, left, right, colors="xkcd:black", linewidth=1)
        plt.savefig(plot_file)
        plt.close()

    @staticmethod
    def flux_residuals_plot_data(p, plot_filename, measured_flux):
        title = f"Flux Residuals (observed minus computed)"
        x_label = "Datapoints"
        y_label = "Normalized Flux"
        x = [[x for x in range(measured_flux.shape[0])]]
        y = [measured_flux["residual"]]
        data_labels = ["residual"]
        linestyles = [""]
        markers = ["o"]
        markersizes = [3]
        colors = ["xkcd:nice blue"]
        linewidths = [0]
        plot_file = p.results_directory + plot_filename

        plt.figure(figsize=(10, 6))
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        # plt.legend()
        # plt.grid(True)
        plt.ticklabel_format(useOffset=False, style="plain", axis="x")   # show x-labels as they are

        # for time, residual, jitter in zip(x, y, measured_flux["flux_err"]):
        #     plt.vlines(time, residual - jitter, residual + jitter, colors="xkcd:black", linewidth=1)
        # for time, residual, jitter in zip(x, y, measured_flux["flux_err"]):
        #     plt.vlines(time, residual - jitter, residual + jitter, colors="xkcd:black", linewidth=1)
        # plt.vlines(measured_flux["time"], measured_flux["residual"] - measured_flux["flux_err"], measured_flux["residual"] + measured_flux["flux_err"], colors="xkcd:black", linewidth=1)
        plt.vlines(range(len(measured_flux["time"])), measured_flux["residual"] - measured_flux["flux_err"], measured_flux["residual"] + measured_flux["flux_err"], colors="xkcd:black", linewidth=1)

        plt.plot(x[0], y[0], marker=markers[0], markersize=markersizes[0], linestyle=linestyles[0], label=data_labels[0], color=colors[0], linewidth=linewidths[0])
        plt.hlines(0, x[0][0], x[0][-1], colors="xkcd:black", linewidth=1)
        plt.savefig(plot_file)
        plt.close()



def try_colors_in_plot():
    plot_filename1 = "color_test.png"
    distance = 10
    n_lines = 100
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ["xkcd:royal blue", "xkcd:red", "xkcd:black", "xkcd:frog green", "xkcd:piss yellow", "xkcd:purply blue", "xkcd:sepia", "xkcd:wine", "xkcd:ocean", "xkcd:shit green", "xkcd:forest", "xkcd:pale violet", "xkcd:robin's egg", "xkcd:pinkish purple", "xkcd:azure", "xkcd:hot pink", "xkcd:mango", "xkcd:baby pink", "xkcd:fluorescent green", "xkcd:medium grey"]
    linestyles = ["solid", "dashed", "dashdot", "dotted"]
    for i in range(n_lines):
        color = colors[i % len(colors)]
        linestyle = linestyles[(i // len(colors)) % len(linestyles)]
        x = i * distance
        ax.plot([x, x], [-1, 1], color=color, linestyle=linestyle, linewidth=1)
    ax.set_xlabel("test")
    ax.set_title(f"test")
    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(plot_filename1)


if __name__ == "__main__":
    try_colors_in_plot()
# When testing, do not run this file directly. Run run_curvesim.py (in the parent directory) instead.


# import numpy as np


def _lmfit_worker_queue(task_queue, result_queue):
    for task in iter(task_queue.get, None):
        config_file, time_s0, time_d, measured_tt, p, run_id = task
        p.randomize_startvalues_uniform()
        p.TOI4504_startvalue_hack()
        bodies_local = CurveSimBodies(p)
        lmfit_run = CurveSimLMfit(p, bodies_local, time_s0, time_d, measured_tt)
        try:
            lmfit_run.save_best_fit(p, bodies_local, measured_tt)
        except Exception:
            pass
        result_queue.put(run_id)
        task_queue.task_done()  # Mark as processed


def run_all_queue(tasks, max_workers):
    task_queue = JoinableQueue()
    result_queue = JoinableQueue()
    total = len(tasks)
    next_index = 0

    # start workers
    workers = [Process(target=_lmfit_worker_queue, args=(task_queue, result_queue))
               for _ in range(max_workers)]
    for w in workers:
        w.start()
        time.sleep(0.2)

    # submit up to max_workers initial tasks
    for _ in range(min(max_workers, total)):
        task_queue.put(tasks[next_index])
        next_index += 1

    completed = 0
    while completed < total:
        run_id = result_queue.get()
        completed += 1
        print(f"LMfit run {run_id} finished ({completed}/{total})")
        # immediately submit the next pending task (if any) so a freed worker starts a new run
        if next_index < total:
            task_queue.put(tasks[next_index])
            next_index += 1

    # wait for workers to mark all tasks done
    task_queue.join()

    # stop workers
    for _ in workers:
        task_queue.put(None)
    for w in workers:
        w.join()


class CurveSimulator:
    def __init__(self, config_file=""):
        warnings.filterwarnings("ignore", module="rebound")
        p = CurveSimParameters(config_file)  # Read program parameters from config file.
        bodies = None
        if p.verbose:
            print(p)
        if p.action in ["lmfit", "guifit", "mcmc"]:
            measured_flux_array, flux_uncertainty, measured_tt, time_s0, time_d, tt_s0, tt_d = (None,) * 7
            if p.flux_file:
                time_s0, time_d, measured_flux_array, flux_uncertainty, measured_flux = CurveSimResults.get_measured_flux(p)
            elif p.tt_file:
                time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
                measured_tt = CurveSimResults.get_measured_tt(p)
            bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
            for body in bodies:  # HACK because length of body.positions is initialized with the correct value for simulation, NOT measurements
                body.positions = np.ndarray((len(time_s0), 3), dtype=float)
            p.init_fitting_parameter_dic()
            print(f"Fitting {len(p.fitting_parameters)} parameters.")
            if p.action == "guifit":
                p.enrich_fitting_params(bodies)
                self.guifit = CurveSimManualFit(p, bodies, time_s0, time_d, measured_tt)
                self.guifit.save_lmfit_results(p)
                sys.exit(0)
            elif p.action == "lmfit":
                num_workers = max(1, os.cpu_count() - 1)  # number of parallel lmfit runs (multiprocessing). Leave one CPU availabe for other programs
                total_runs = 1000  # total number of lmfit runs
                print(f"{num_workers=}, {total_runs=}")
                while True:
                    tasks = [(config_file, time_s0, time_d, measured_tt, p, i) for i in range(total_runs)]
                    run_all_queue(tasks, num_workers)
            if p.action == "mcmc":
                mcmc = CurveSimMCMC(p, bodies, time_s0, time_d, measured_flux_array, flux_uncertainty, measured_tt)
                self.sampler = mcmc.sampler  # mcmc object
                self.theta = mcmc.theta  # current state of mcmc chains. By saving sampler and theta it is possible to continue the mcmc later on.
            else:
                print(f"{Fore.RED}\nERROR: Invalid value for parameter <action> in configuration file {Style.RESET_ALL}")
                sys.exit(1)
        elif p.action == "single_run":
            bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
            # bodies.save(p, "vorne_", "_hinten")
            # new_body = CurveSimBody.load("vorne_TOI4504d_hinten.bdy")
            # new_body.save("abc__")
            # exit(1)
            bodies, sim_flux, results = CurveSimMCMC.single_run(p, bodies)
            self.sim_flux = sim_flux
            self.results = results
        elif p.action == "results_only":
            # time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
            # bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation

            # p.eclipsers = ["TOI4504c"]
            # p.eclipsees = ["TOI4504"]

            # results = CurveSimResults.load_results(p.result_file)
            # tt_sim = results.get_transit_data("TOI4504c", "TOI4504", "TT")
            # print(tt_sim)

            CurveSimResults.ttv_to_date_plot(p, amplitude=2.1, period=965, x_offset=-340, osc_per=82.97213)
            CurveSimResults.ttv_to_date_plot(p, amplitude=2.1, period=965, x_offset=-450, osc_per=82.83)
            CurveSimResults.ttv_to_date_plot(p, amplitude=2.08, period=965, x_offset=-449, osc_per=82.834)
            CurveSimResults.ttv_to_date_plot(p, amplitude=2.0, period=946.5, x_offset=-393, osc_per=82.5438)
            sys.exit(0)
        elif p.action == "get_tess_data":
            CurveSimFluxData.get_tess_flux(p)
        elif p.action == "process_tess_data":
            CurveSimFluxData.process_tess_flux(p)
        else:
            print(f"{Fore.RED}\nERROR: Invalid value for parameter <action> in configuration file {Style.RESET_ALL}")
            sys.exit(1)
        self.parameters = p
        self.bodies = bodies

    def __repr__(self):
        print("CurveSimulator object created with attributes: ", end="")
        for key in self.__dict__:
            print(key, end=", ")
        # print(self.__dict__.keys())
        return ""
# Short version: x
# from curvesimulator import curvesim
# curvesim("../configurations/MyFirstConfigFile.ini")
#
#
# Long version :
# from curvesimulator import curvesim
#
# def main():
#     parameters, bodies, results, sim_flux = curvesim(config_file="../configurations/MyFirstConfigFile.ini")
#     print(parameters)
#     print(bodies)
#     print(results)
#     print(sim_flux)
#
#
# if __name__ == "__main__":
#     main()


# Compatibility shim: some older libraries expect numpy.in1d which may be removed in newer NumPy versions.
# Provide a safe alias to numpy.isin before importing other packages (e.g., astropy/lightkurve) that use it.
# import numpy as np
# if not hasattr(np, 'in1d'):
#     np.in1d = np.isin


# Current version for developers

def main():
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/debug.ini")

    curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_T102_Trifon01.ini")

    # curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V001.ini")
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V002.ini")
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V003.ini")
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V004.ini")
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V005.ini")
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V006.ini")
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V007.ini")

    print(curvesimulation)


if __name__ == "__main__":
    main()
