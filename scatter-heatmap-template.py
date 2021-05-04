# !/usr/bin/env python3
# HW 1 (04/16/2021)
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
            description='Homework 2: generate figures, generate panels, scatter plot, bar plot, optional heatmap',
            add_help=True,
            prefix_chars='-',
            usage='%(prog)s --input_file (-i) --output_file (-o) --style_sheet (-s)')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
        self.parser.add_argument('-i', '--input_file', action='store', help='input file name')
        self.parser.add_argument('-o', '--output_file', action='store', help='output file name')
        self.parser.add_argument('-s', '--style_sheet', action='store', help='style sheet file name')

        if in_opts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(in_opts)


class Figures:
    """"""
    def __init__(self, input_data, style_sheet, output_file):
        """ Constructor: saves data from input files """
        self.input_data = input_data
        self.output_file = output_file
        self.style_sheet = style_sheet
        self.data = []

    def parse_data(self):
        """ Parse input data """
        fn = open(self.data, "r")
        for line in fn:
            data = line.strip().split()
            if data[0].startswith('ENSMUSG'):
                self.data.append((data[1], data[2]))
        fn.close()

    def create_figures(self):
        """ Generate figure """

        plt.style.use(self.style_sheet)
        figure_height, figure_width = 2, 5
        plt.figure(figsize=(figure_width, figure_height))

        # panel 1 cluster specs
        panel_1_left = plt.axes([0.076, 0.15, 1 / (8 * 2.5), 1 / 2])
        panel_1_left.set_yticks(range(0, 16, 5))
        panel_1_left.set_xticks(range(0, 21, 20))
        panel_1_left.invert_xaxis()
        panel_1_main = plt.axes([0.14, 0.15, 1 / 5, 1 / 2])
        panel_1_main.set_yticks([])
        panel_1_main.set_xticks(range(0, 16, 5))
        panel_1_main.set_xlim(0, 15)
        panel_1_main.set_ylim(0, 15)
        panel_1_top = plt.axes([0.14, 0.685, 1 / 5, 1 / 8])
        panel_1_top.set_yticks(range(0, 21, 20))
        panel_1_top.set_xticks([])
        panel_1_top.set_xlim(0, 15)
        panel_1_top.set_ylim(0, 20)

        # panel 2 cluster specs
        panel_2_left = plt.axes([0.476, 0.15, 1 / (8 * 2.5), 1 / 2])
        panel_2_left.set_yticks(range(0, 16, 5))
        panel_2_left.set_xticks(range(0, 21, 20))
        panel_2_left.invert_xaxis()
        panel_2_main = plt.axes([0.54, 0.15, 1 / 5, 1 / 2])
        panel_2_main.set_yticks([])
        panel_2_main.set_xticks(range(0, 16, 5))
        panel_2_main.set_xlim(0, 15)
        panel_2_main.set_ylim(0, 15)
        panel_2_top = plt.axes([0.54, 0.685, 1 / 5, 1 / 8])
        panel_2_top.set_yticks(range(0, 21, 20))
        panel_2_top.set_xticks([])
        panel_2_top.set_xlim(0, 15)
        panel_2_top.set_ylim(0, 20)
        # Heat map
        heat_panel = plt.axes([0.8, 0.15, 1 / 50, 1 / 2])
        heat_panel.set_yticks(range(0, 21, 10))
        labels = ['0', '10', '>20']
        heat_panel.set_yticklabels(labels)
        heat_panel.set_xticks([])
        heat_panel.set_xlim(0, 1)
        heat_panel.set_ylim(0, 20)

        # Parse data
        data_x, data_y = [], []
        fn = open(self.input_data, "r")
        for line in fn:
            data = line.strip().split()
            if data[0].startswith('ENSMUSG'):
                data_x.append(int(data[1]))
                data_y.append(int(data[2]))
        fn.close()

        # reformatting data
        xVal = np.array(data_x)
        yVal = np.array(data_y)
        x_val = np.log2(xVal + 1)
        y_val = np.log2(yVal + 1)

        # Main panel code
        panel_1_main.plot(x_val, y_val,
                          marker='o',
                          markerfacecolor='black',
                          markersize=1.5,
                          markeredgewidth=0,
                          linewidth=0,
                          alpha=0.1)

        panel_2_main.plot(x_val, y_val,
                          marker='s',
                          markerfacecolor='black',
                          markersize=1.5,
                          markeredgewidth=0,
                          linewidth=0,
                          alpha=0.08)

        hist_left, bins_y = np.histogram(y_val, np.arange(0, 15, 0.5))
        hist_top, bins_x = np.histogram(x_val, np.arange(0, 15, 0.5))

        # Left panel code
        for i in range(len(hist_left)):
            rectangle_left1 = mplpatches.Rectangle([0, bins_y[i]], np.log2(hist_left[i] + 1), 0.5,
                                                   facecolor='grey',
                                                   edgecolor='black',
                                                   linewidth=0.1)
            panel_1_left.add_patch(rectangle_left1)
            rectangle_left2 = mplpatches.Rectangle([0, bins_y[i]], np.log2(hist_left[i] + 1), 0.5,
                                                   facecolor='grey',
                                                   edgecolor='black',
                                                   linewidth=0.1)
            panel_2_left.add_patch(rectangle_left2)

        # top panel code
        for i in range(len(hist_top)):
            rectangle_top1 = mplpatches.Rectangle([bins_x[i], 0], 0.5, np.log2(hist_top[i] + 1),
                                                  facecolor='grey',
                                                  edgecolor='black',
                                                  linewidth=0.1)
            panel_1_top.add_patch(rectangle_top1)
            rectangle_top2 = mplpatches.Rectangle([bins_x[i], 0], 0.5, np.log2(hist_top[i] + 1),
                                                  facecolor='grey',
                                                  edgecolor='black',
                                                  linewidth=0.1)
            panel_2_top.add_patch(rectangle_top2)

        # Heat Map code scale
        black = (1, 1, 1)
        white = (0, 0, 0)
        R = np.linspace(black[0], white[0], 20)
        G = np.linspace(black[1], white[1], 20)
        B = np.linspace(black[2], white[2], 20)
        for i in range(0, 20):
            heat = mplpatches.Rectangle([0, i], 1, 1,
                                        facecolor=(R[i], G[i], B[i]),
                                        edgecolor='black',
                                        linewidth=0)
            heat_panel.add_patch(heat)

        # Save figure
        plt.savefig(self.output_file, dpi=600)


def main(command=None):
    """ Generate specified figures using argparse with input and output file specifications. """
    if command is None:
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main

    input_file = my_command.args.input_file
    style_sheet = my_command.args.style_sheet
    output_file = my_command.args.output_file

    # create plots
    my_figure = Figures(input_file, style_sheet, output_file)
    my_figure.create_figures()


if __name__ == "__main__":
    main()


# Python3 Mason_Zachary_BME163_Assignment_Week2.py -i BME163_Input_Data_1.txt -o Mason_Zachary_BME163_Assignment_Week2.png -s BME163