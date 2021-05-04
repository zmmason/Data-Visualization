# !/usr/bin/env python3
# (04/09/2021)
# University of California, Santa Cruz - BME 163
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason@ucsc.edu)

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches


class CommandLine:
    """ Handle the command line, usage and help requests """
    def __init__(self, in_opts=None):
        """ CommandLine constructor: Implements a parser to interpret the command line argv string using argparse. """
        self.parser = argparse.ArgumentParser(
            description='Homework 1: generate figures, generate panels, understand RGB color, using plot and '
                        'rectangular functions, manipulate axes ticks and labels, use Matlibplot style sheets, '
                        'geometry',
            add_help=True,
            prefix_chars='-',
            usage='%(prog)s --style_sheet (-s) --output_file (-o)')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
        self.parser.add_argument('-s', '--style_sheet', action='store', help='input file name')
        self.parser.add_argument('-o', '--output_file', action='store', help='output file name')
        if in_opts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(in_opts)


class Figures:
    """"""
    def __init__(self, style_sheet, output_file):
        """ Constructor: saves data from input files """
        self.style_sheet = style_sheet
        self.output_file = output_file

    def create_figures(self):
        """ Generate figure """

        plt.style.use(self.style_sheet)
        figure_height, figure_width = 2, 3.42
        plt.figure(figsize=(figure_width, figure_height))

        # panel spec
        panel_height, panel_width = 1, 1

        # figure-panel spec
        relative_panel_width = panel_width / figure_width
        relative_panel_height = panel_height / figure_height

        # panel 1
        # color
        black = (0, 0, 0)
        white = (1, 1, 1)

        panel_1 = plt.axes([0.1, 0.2, relative_panel_width, relative_panel_height])

        # 25 equally spaced 1st quadrant (pi/2) points
        points = np.linspace(0, np.pi/2, 25)

        # plot points with gradient in first panel
        for i in range(len(points)):
            y_val = np.cos(points[i])
            x_val = np.sin(points[i])
            colorVal = abs(x_val)  # color gradient point
            panel_1.plot(x_val, y_val,
                         marker='o',
                         markerfacecolor=(colorVal, colorVal, colorVal),
                         markersize=2,
                         markeredgewidth=0,
                         linewidth=0)

        # remove tick marks
        plt.xticks([])
        plt.yticks([])

        # panel 2
        # color
        blue = (0, 0, 1)
        cyan = (0, 230/255, 1)  # RGB vals from pycharm color picker
        snow = (230/255, 230/255, 1)  # RGB vals from pycharm color picker
        pink = (230/255, 0, 1)  # RGB vals from pycharm color picker

        panel_2 = plt.axes([.55, 0.2, relative_panel_width, relative_panel_height])

        # bottom row gradient (blue -> pink)
        R = np.linspace(blue[0], pink[0], 10)
        G = np.linspace(blue[1], pink[1], 10)
        B = np.linspace(blue[2], pink[2], 10)
        # top row gradient (cyan -> snow)
        R2 = np.linspace(cyan[0], snow[0], 10)
        G2 = np.linspace(cyan[1], snow[1], 10)
        B2 = np.linspace(cyan[2], snow[2], 10)

        # create panels for RGB gradient
        for i in range(0, 10):  # get first and last row of RGB patches
            row_top = mplpatches.Rectangle([i, 9], 1, 1,
                                           facecolor=(R2[i], G2[i], B2[i]),
                                           edgecolor='black',
                                           linewidth=1.0)
            row_bottom = mplpatches.Rectangle([i, 0], 1, 1,
                                              facecolor=(R[i], G[i], B[i]),
                                              edgecolor='black',
                                              linewidth=1.0)
            panel_2.add_patch(row_top)
            panel_2.add_patch(row_bottom)

            # column gradients between top and bottom row for the ith patch index
            R3 = np.linspace(R[i], R2[i], 10)
            G3 = np.linspace(G[i], G2[i], 10)
            B3 = np.linspace(B[i], B2[i], 10)
            for j in range(1, 9):  # Get columns and RGB of patch between
                cols = mplpatches.Rectangle([i, j], 1, 1,
                                            facecolor=(R3[j], G3[j], B3[j]),
                                            edgecolor='black',
                                            linewidth=1.0)
                panel_2.add_patch(cols)

        # removing tick marks
        plt.xticks([])
        plt.yticks([])

        # set panel 2 limits
        panel_2.set_xlim(0, 10)
        panel_2.set_ylim(0, 10)

        # plot
        plt.savefig(self.output_file, dpi=600)


def main(command=None):
    """ Generate specified figures using argparse with input and output file specifications. """
    if command is None:
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main

    style_sheet = my_command.args.style_sheet
    output_file = my_command.args.output_file

    # create plots
    my_figure = Figures(style_sheet, output_file)
    my_figure.create_figures()


if __name__ == "__main__":
    main()
