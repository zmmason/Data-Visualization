# !/usr/bin/env python3
# (05/03/2021)
# University of California, Santa Cruz - BME 163
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason@ucsc.edu)

import argparse
import numpy as np
import random
import matplotlib.pyplot as plt


class CommandLine:
    """ Handle the command line, usage and help requests """
    def __init__(self, in_opts=None):
        """ CommandLine constructor: Implements a parser to interpret the command line argv string using argparse. """
        self.parser = argparse.ArgumentParser(
            description='Homework 4',
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
    """create figure"""
    def __init__(self, input_file, output_file, style_sheet):
        """ Constructor: saves data from input files """
        self.input_file = input_file
        self.output_file = output_file
        self.style_sheet = style_sheet
        self.subsample = self.parser()

    def parser(self):
        """Parse input data to get complete subsample lists"""
        fn = open(self.input_file, 'r')
        subsample = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: []}
        for line in fn:
            line = line.strip().split()
            subread = int(line[0].split('_')[3])
            if subread >= 11:
                subread = 11
            percent = float(line[1])
            subsample[subread].append(percent)
        fn.close()
        return self.random_subsample(subsample)

    def random_subsample(self, subsample):
        """Randomly subsample lists"""
        subsample2 = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: []}
        medians = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0}
        iteration = 1000  # random choice subsample
        for i in subsample:
            for j in range(iteration):
                subsample2[i].append(random.choice(subsample[i]))
            median = sum(subsample2[i]) / iteration
            medians[i] = float(median)
        return subsample2

    def swarmplot(self, bin, subset, pt_size):
        """Create swarmplot"""
        used = []
        diam = pt_size/72
        for y_pos in subset:
            x_pos = bin
            proximity = []
            iterations = 0
            placed = False
            for point in used:
                new_pt = point[1]*(2/25)
                new_ypt = y_pos*(2/25)
                if np.abs(new_pt-new_ypt) <= diam:
                    proximity.append(point)
            if len(proximity) == 0:
                used.append((x_pos, y_pos))
                continue
            while not placed:
                proximity_list = []
                iterations += 1
                shift = diam/5
                for point in proximity:
                    # check distance from current point
                    distance = np.sqrt((x_pos*(5/11.5) - point[0]*(5/11.5))**2 + (y_pos*(2/25) - point[1]*(2/25))**2)
                    proximity_list.append(distance)
                min_dist = min(proximity_list)
                if min_dist < diam:
                    # moves point left/right depending on iteration
                    if iterations % 2 == 0:
                        center_shift = shift*(iterations/2)
                    else:
                        center_shift = -(shift * (iterations/2))
                    x_pos = center_shift + x_pos
                elif min_dist > diam:
                    used.append((x_pos, y_pos))
                    placed = True
        return used

    def create_figures(self):
        """ Generate figure """
        plt.style.use(self.style_sheet)
        figure_height, figure_width = 3, 7
        plt.figure(figsize=(figure_width, figure_height))

        # panel spec
        panel_height, panel_width = 2, 5

        # figure-panel spec
        relative_panel_width = panel_width / figure_width
        relative_panel_height = panel_height / figure_height

        # panel 1 ticks/labels
        panel_1 = plt.axes([1/10, 1/5, relative_panel_width, relative_panel_height])
        panel_1.set_xlabel('Subread Coverage')
        panel_1.set_ylabel('Identity %')
        panel_1.set_xticks(range(1, 12, 1))
        labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '>11']
        panel_1.set_xticklabels(labels)
        panel_1.set_yticks(range(75, 101, 5))
        panel_1.set_ylim(75, 100)
        panel_1.set_xlim(0.25, 11.75)

        # line at 95%
        panel_1.plot([-2, 14], [95, 95],
                     lw=0.5,
                     dashes=(1, 2, 2, 2),
                     c='black')

        # call swarm plot function
        for bin in self.subsample:
            subset_points = self.swarmplot(bin, self.subsample[bin], pt_size=0.65)
            y_coord = [i[1] for i in subset_points]
            x_coord = [i[0] for i in subset_points]
            panel_1.plot(x_coord, y_coord,
                         marker='o',
                         markerfacecolor='black',
                         markersize=0.65,
                         markeredgewidth=0,
                         linewidth=0)
            # Median plotting
            panel_1.plot([bin-panel_width/12.5, bin+panel_width/12.5], [np.median(y_coord), np.median(y_coord)],
                         linewidth=1,
                         color='red')

        plt.savefig(self.output_file, dpi=600)


def main(command=None):
    """ Generate specified figures using argparse with input and output file specifications. """
    if command is None:
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main

    input_file = my_command.args.input_file
    output_file = my_command.args.output_file
    style_sheet = my_command.args.style_sheet

    # create plots
    my_figure = Figures(input_file, output_file, style_sheet)
    my_figure.create_figures()


if __name__ == "__main__":
    main()

# python3 Mason_Zachary_BME163_Assignment_Week4.py -i BME163_Input_Data_3.txt -o TEST.png -s BME163
