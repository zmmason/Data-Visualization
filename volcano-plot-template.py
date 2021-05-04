# !/usr/bin/env python3
# HW 3 (04/23/2021)
# University of California, Santa Cruz - BME 163
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason@ucsc.edu)

import argparse
import numpy as np
import matplotlib.pyplot as plt


class CommandLine:
    """ Handle the command line, usage and help requests """
    def __init__(self, in_opts=None):
        """ CommandLine constructor: Implements a parser to interpret the command line argv string using argparse. """
        self.parser = argparse.ArgumentParser(
            description='Homework 3',
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
    def __init__(self, input_file, output_file, style_sheet):
        """ Constructor: saves data from input files """
        self.input_file = input_file
        self.output_file = output_file
        self.style_sheet = style_sheet

    def create_figures(self):
        """ Generate figure """
        plt.style.use(self.style_sheet)
        figure_height, figure_width = 3, 3
        plt.figure(figsize=(figure_width, figure_height))

        # panel spec
        panel_height, panel_width = 2, 2

        # figure-panel spec
        relative_panel_width = panel_width / figure_width
        relative_panel_height = panel_height / figure_height

        # panel 1 ticks/labels
        panel_1 = plt.axes([1/6, 1/6, relative_panel_width, relative_panel_height])
        panel_1.set_xlabel('log$_{2}$(fold change)')
        panel_1.set_ylabel('-log$_{10}$(p-value)')
        panel_1.set_xticks([-10, -5, 0, 5, 10])
        panel_1.set_yticks(range(0, 61, 10))
        panel_1.set_ylim(0, 60)
        panel_1.set_xlim(-12, 12)

        # parser
        file_data = open(self.input_file, 'r')
        all_x = []  # hold all fold changes
        all_y = []  # hold all p-values
        red_gene_x = []  # hold all red x vals
        red_gene_y = []  # hold all red y vals
        red_label_gene = {}  # hold labeled data points
        for line in file_data:
            data = line.strip().split()
            # only get genes with valid data (omit points with 'NA' - BME 163 April 17, 2020 (37:25) )
            if data[1] and data[2] != 'NA':
                gene, fold_change, p_val = data[0], float(data[1]), -np.log10(float(data[2]))
                all_x.append(fold_change)
                all_y.append(p_val)
                # conditionals
                if (2**abs(fold_change)) > 10 and p_val > 8:  # get all red points
                    red_gene_x.append(fold_change)
                    red_gene_y.append(p_val)
                if -(2**fold_change) > -10 and p_val > 30:  # get all labeled points ('mt-Nd1', 'mt-Nd4', 'Gdf11', 'Ccdc24')
                    red_label_gene[gene] = [fold_change, p_val]
        file_data.close()

        # plot all data
        panel_1.plot(all_x, all_y,
                     marker='o',
                     markerfacecolor='black',
                     markersize=2**0.5,
                     markeredgewidth=0,
                     linewidth=0)

        # plot red points over specific points
        panel_1.plot(red_gene_x, red_gene_y,
                     marker='o',
                     markerfacecolor='red',
                     markersize=1.4,
                     markeredgewidth=0,
                     linewidth=0)

        # adding labels to specific significant points
        for point in red_label_gene.items():
            gene = point[0]
            x_coord = point[1][0]
            y_coord = point[1][1]
            panel_1.text(x_coord-0.26, y_coord-0.55,  # adjusting placement of the labels for visual
                         gene,
                         horizontalalignment='right',
                         fontsize=6)

        # plot
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
