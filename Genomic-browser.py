# !/usr/bin/env python3
# Genomic Browser
# University of California, Santa Cruz
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason@ucsc.edu)
# python3 Genomic-browser.py -i1 BME163_Input_Data_5.psl -i2 BME163_Input_Data_6.psl -g gencode.vM12.annotation.gtf -o Mason_Zachary_BME163_Assignment_Final.png

import argparse
import sys

import numpy as np
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches


class CommandLine:
    """ Handle the command line, usage and help requests """

    def __init__(self, in_opts=None):
        """ CommandLine constructor: Implements a parser to interpret the command line argv string using argparse """
        self.parser = argparse.ArgumentParser(
            description='Visually stack reads associated with target chromosome/genomic regions.',
            add_help=True,
            prefix_chars='-',
            usage='python3 %(prog)s --psl1 (-i1) --psl2 (-i2) --gtf (-g) --output_file (-o)')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
        self.parser.add_argument('-i1', '--psl1', action='store', help='/location/of/Input_PSL1.psl')
        self.parser.add_argument('-i2', '--psl2', action='store', help='/location/of/Input_PSL2.psl')
        self.parser.add_argument('-g', '--gtf', action='store', help='/location/of/gencode_GTF.gtf')
        self.parser.add_argument('-o', '--output_file', action='store', help='output file name')

        if in_opts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(in_opts)


def gtf_parser(gtf, target):
    """ Parse .gtf files """
    gtfDict = {}
    open_file = open(gtf, 'r')
    for line in open_file:
        if line[0] != '#':
            read_data = line.strip().split('\t')
            chromosome = read_data[0]
            if chromosome == target[0]:  # verify chromosome 7
                feature = read_data[2]
                if feature in ['exon', 'CDS']:
                    start = int(read_data[3])
                    end = int(read_data[4])
                    if target[1] < start < target[2] or target[1] < end < target[2]:  # verify read in target region
                        transcript = read_data[8].split(' transcript_id "')[1].split('"')[0]
                        if transcript not in gtfDict:
                            gtfDict[transcript] = [[chromosome, start, end, feature]]
                        else:
                            gtfDict[transcript].append([chromosome, start, end, feature])
                    else:
                        continue

    transcript_list = []
    for title, transcripts in gtfDict.items():
        starts = []
        ends = []
        block_starts = []
        block_widths = []
        features = []
        for part in transcripts:
            starts.append(part[1])
            ends.append(part[2])
            block_starts.append(part[1])
            block_widths.append(part[2] - part[1])
            chromosome = part[0]
            features.append(part[3])
        transcript_list.append([chromosome, min(starts), max(ends), block_starts, block_widths, False, features])
    transcript_list.sort(key=itemgetter(2))
    return transcript_list


def psl_parser(psl, target):
    """ Parse .psl files (specific for target gene region: chr7 45232945:45240000)"""
    read_list = []
    open_file = open(psl, 'r')
    for line in open_file:
        read_data = line.strip().split('\t')
        chromosome = read_data[13]
        if chromosome == target[0]:  # verify chromosome 7
            start = int(read_data[15])
            stop = int(read_data[16])
            if target[1] < start < target[2] and target[1] < stop < target[2]:  # verify read in target region
                block_start = np.array(read_data[20].split(',')[:-1], dtype=int)
                block_width = np.array(read_data[18].split(',')[:-1], dtype=int)
                read = [chromosome, start, stop, block_start, block_width, False]
                read_list.append(read)
        else:
            continue

    read_list.sort(key=itemgetter(2))

    return read_list


class Figures:
    """ Create figure components """

    def __init__(self, output_file, psl1_reads, psl2_reads, gtf_data, target):
        """ Constructor: saves data from input files """
        self.output_file = output_file
        self.psl1_reads = psl1_reads
        self.psl2_reads = psl2_reads
        self.start = target[1]
        self.stop = target[2]
        self.gtf_data = gtf_data

    def set_panel(self, figure_width, figure_height, x_panel, y_panel):
        """ set up panel specs"""
        panel_width, panel_height = 10, 1.25
        relative_panelWidth = panel_width / figure_width
        relative_panelHeight = panel_height / figure_height
        panel = plt.axes([x_panel, y_panel, relative_panelWidth, relative_panelHeight])
        panel.set_xticks([])
        panel.tick_params(bottom=False, labelbottom=False,
                          left=False, labelleft=False,
                          right=False, labelright=False,
                          top=False, labeltop=False)
        return panel

    def create_fig(self):
        """ Create figure """
        # figure specs
        plt.style.use('BME163')
        figure_width, figure_height = 10, 5
        plt.figure(figsize=(figure_width, figure_height))

        # top panel
        y_top = 0.65
        top_panel = self.set_panel(figure_width, figure_height, 0, y_top)
        plot_gtf(self.gtf_data, top_panel)
        top_panel.set_xlim(self.start, self.stop)
        top_panel.set_ylim(0.24, 10.0)

        # middle panel
        y_middle = 0.35
        middle_panel = self.set_panel(figure_width, figure_height, 0, y_middle)
        plot_psl(self.psl2_reads, middle_panel)
        middle_panel.set_xlim(self.start, self.stop)
        middle_panel.set_ylim(0.2, 69.5)

        # bottom panel
        y_bottom = 0.05
        bottom_panel = self.set_panel(figure_width, figure_height, 0, y_bottom)
        plot_psl(self.psl1_reads, bottom_panel)
        bottom_panel.set_xlim(self.start, self.stop)
        bottom_panel.set_ylim(0, 436)

        plt.savefig(self.output_file, dpi=1200)


def plot_gtf(data, panel):
    """ Plot data for each panel """
    thin_line, thick_line, gap = 0.1, 0.25, 0.18
    for y_pos in range(1, len(data)):
        plotted_end = 0
        for read in data:
            start, end, block_starts, block_widths, plotted, feature = read[1], read[2], read[3], read[4], read[5], read[6]
            if plotted is False:
                if start > plotted_end:
                    rectangle = mplpatches.Rectangle((start, y_pos + gap), end - start, thin_line,
                                                     facecolor='black',
                                                     edgecolor='black',
                                                     linewidth=0)
                    panel.add_patch(rectangle)
                    for index in np.arange(0, len(block_starts), 1):
                        block_start = block_starts[index]
                        block_width = block_widths[index]
                        the_feature = feature[index]
                        if the_feature == 'exon':
                            rectangle1 = mplpatches.Rectangle((block_start, y_pos + 0.12), block_width, thick_line,
                                                              facecolor='black',
                                                              edgecolor='black',
                                                              linewidth=0)
                            panel.add_patch(rectangle1)
                        if the_feature == 'CDS':
                            thick_line2 = 0.5
                            rectangle2 = mplpatches.Rectangle((block_start, y_pos - 0.01), block_width, thick_line2,
                                                              facecolor='black',
                                                              edgecolor='black',
                                                              linewidth=0)
                            panel.add_patch(rectangle2)

                    read[5] = True
                    plotted_end = end


def plot_psl(data, panel):
    """ Plot data for each panel """
    thin_line, thick_line, gap = 0.1, 0.5, 0.18
    for y_pos in range(1, len(data)):
        plotted_end = 0
        for read in data:
            start, end, block_starts, block_widths, plotted = read[1], read[2], read[3], read[4], read[5]
            if plotted is False:
                if start > plotted_end:
                    rectangle = mplpatches.Rectangle((start, y_pos + gap), end - start, thin_line,
                                                     facecolor='black',
                                                     edgecolor='black',
                                                     linewidth=0)
                    panel.add_patch(rectangle)
                    for index in np.arange(0, len(block_starts), 1):
                        block_start = block_starts[index]
                        block_width = block_widths[index]
                        rectangle1 = mplpatches.Rectangle((block_start, y_pos), block_width, thick_line,
                                                          facecolor='black',
                                                          edgecolor='black',
                                                          linewidth=0)
                        panel.add_patch(rectangle1)
                    read[5] = True
                    plotted_end = end


def main(command=None):
    """ Generate specified figures using argparse with input and output file specifications. """

    if command is None:
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main

    target_chromosome = input("Target Chromosome: ")
    genomic_region_start = input("Target starting location: ")
    genomic_region_end = input("Target ending location: ")
    target = [target_chromosome, int(genomic_region_start), int(genomic_region_end)]
    # target = ['chr7', 45232945, 45240000]  # target region
    print('Generating figure...')
    # get commandline arguments
    i1 = my_command.args.psl1
    i2 = my_command.args.psl2
    gtf = my_command.args.gtf
    output_file = my_command.args.output_file

    # parse data
    psl1_reads = psl_parser(i1, target)
    psl2_reads = psl_parser(i2, target)
    gtf_data = gtf_parser(gtf, target)

    # plot parsed data
    my_figure = Figures(output_file, psl1_reads, psl2_reads, gtf_data, target)
    my_figure.create_fig()
    print('Done')


if __name__ == "__main__":
    main()
