#!/usr/bin/env python
# -*- coding: utf-8 -*-
from SPAdesPipeline.OLCspades.accessoryFunctions import *
__author__ = 'adamkoziol'


class SNPprimers(object):

    def snpprimers(self):
        """
        Run the necessary analyses
        """
        # Extract the coordinates from the file
        self.extractcoords()
        # Load the contig records into a dictionary and extract desired amplicons from the sequences
        self.loadcontigs()
        # Run primer3 to create primers from the desired amplicons
        self.primerthreads()
        # Create reports
        self.reports()
        # Print the metadata to file
        self.printmetadata()

    def extractcoords(self):
        """
        Read the file with contigs and SNP coordinates into memory
        """
        printtime(u'Loading SNP locations', self.start)
        with open(self.coords, 'rb') as coords:
            for line in coords:
                # Split the line on the delimiter
                data = line.split(self.delimiter)
                # Add each SNP to a list of SNPs for each contig
                try:
                    self.contigcoords[data[0]].append(int(data[1].rstrip()))
                except KeyError:
                    self.contigcoords[data[0]] = list()
                    self.contigcoords[data[0]].append(int(data[1].rstrip()))

    def loadcontigs(self):
        """
        Load all the contigs in the assembly .fasta file into a dictionary, and extract the amplicon (upstream and
        downstream sequence) for each SNP
        """
        printtime(u'Extracting SNP context', self.start)
        from Bio import SeqIO
        # Initialise a dictionary to store the SeqIO records of each assembly
        record = dict()
        # Load the records from the assembly into the dictionary
        for rec in SeqIO.parse(self.fastafile, 'fasta'):
            record[rec.id] = str(rec.seq)

        # Iterate through the contigs in the dictionary
        for contig, positions in sorted(self.contigcoords.items()):
            # Iterate through all the SNPs in each contig
            for position in sorted(positions):
                # Initialise the metadata object to store all necessary information
                metadata = MetadataObject()
                metadata.name = '{}_{}'.format(contig, position)
                metadata.amplicon = GenObject()
                # The start and end of the desired amplicon will be the SNP position will be +/- half the desired
                # amplicon size +/- an additional 50 (51 to allow for zero-based indexing) bp to give primer3
                # enough sequence data to work with
                metadata.amplicon.start = position - (self.ampliconsize / 2) - 51
                metadata.amplicon.end = position + (self.ampliconsize / 2) + 50
                # Extract the gene sequence from the contigs
                # The record dictionary has the contig name, and the sequence. Splice out the data using the start and
                # end coordinates specified by ePCR
                metadata.amplicon.sequence = record[contig][metadata.amplicon.start: metadata.amplicon.end]
                # Add the metadata object to a list of metadata objects
                self.metadata.append(metadata)

    def primerthreads(self):
        """
        Multithreaded approach to primer3-mediated primer creation
        """
        printtime(u'Creating primers', self.start)
        from threading import Thread
        # Create the threads for primer creation
        for i in range(len(self.metadata)):
            threads = Thread(target=self.primers, args=())
            threads.setDaemon(True)
            threads.start()
        for sample in self.metadata:
            # Create an object to store the parsed results
            sample.primerzero = GenObject()
            self.primerqueue.put(sample)
        # Join the threads
        self.primerqueue.join()

    def primers(self):
        from subprocess import Popen, PIPE, STDOUT
        while True:
            sample = self.primerqueue.get()
            # Create a Boulder-IO compatible string to be used as input into primer3
            inputdata = 'SEQUENCE_ID={}\n' \
                        'SEQUENCE_TEMPLATE={}\n' \
                        'PRIMER_PRODUCT_SIZE_RANGE={}-{}\n' \
                        '='\
                .format(sample.name, sample.amplicon.sequence, self.ampliconsize - 50, self.ampliconsize + 50)
            # Set the output file name
            outputfile = '{}/{}.txt'.format(self.resultspath, sample.name)
            # Set the system command
            command = ['primer3_core']
            # Set up the system call
            systemcall = Popen(command, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
            # Run the system call with the Boulder-IO input
            data = systemcall.communicate(input=inputdata)
            # Get the data into a dictionary
            primerresults = data[0]
            # Extract only data containing '_0' (e.g. PRIMER_RIGHT_0 or PRIMER_PAIR_0_PRODUCT_SIZE)
            primerzero = {results.split('=')[0]: results.split('=')[1] for results in
                          primerresults.split('\n') if '_0' in results}
            #
            for description, value in primerzero.items():
                setattr(sample.primerzero, description, value)
            # Write the data to file
            with open(outputfile, 'wb') as report:
                report.write(data[0])
            self.primerqueue.task_done()

    def reports(self):
        """Create reports from the abundance estimation"""
        import xlsxwriter
        printtime(u'Creating report', self.start)
        # Create a workbook to store the report. Using xlsxwriter rather than a simple csv format, as I want to be
        # able to have appropriately sized, multi-line cells
        workbook = xlsxwriter.Workbook('{}/primers.xlsx'.format(self.reportpath))
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 8
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 8})
        bold.set_align('center')
        # Format for data cells. Monotype, size 8, top vertically justified
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 8})
        courier.set_align('top')
        # Set the custom width for 1 and 2 to be 15
        worksheet.set_column(1, 2, 20)
        worksheet.set_column(3, 3, 10)
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        col = 0
        # List of the headers to use
        headers = ['SNP', 'LeftPrimer', 'RightPrimer', 'ProductSize', 'LeftTm', 'RightTm']
        # Populate the headers
        for category in headers:
            # Write the data in the specified cell (row, col) using the bold format
            worksheet.write(row, col, category, bold)
            # Move to the next column to write the next category
            col += 1
        # Data starts in row 1
        row = 1
        # Initialise variables to hold the longest name; used in setting the column width
        longestname = 0
        # Return the primer results for each SNP
        for sample in self.metadata:
            # Every record starts at column 0
            col = 0
            # Write the strain name
            worksheet.write(row, col, sample.name, courier)
            col += 1
            # Determine the longest name of all the strains, and use it to set the width of column 0
            if len(sample.name) > longestname:
                longestname = len(sample.name)
                worksheet.set_column(0, 0, longestname)
            # Write the data to the spreadsheet
            for data in [sample.primerzero.PRIMER_LEFT_0_SEQUENCE,
                         sample.primerzero.PRIMER_RIGHT_0_SEQUENCE,
                         sample.primerzero.PRIMER_PAIR_0_PRODUCT_SIZE,
                         sample.primerzero.PRIMER_LEFT_0_TM,
                         sample.primerzero.PRIMER_RIGHT_0_TM]:
                worksheet.write(row, col, data, courier)
                col += 1
            # Increase the row
            row += 1
        # Close the workbook
        workbook.close()

    def printmetadata(self):
        """
        Print the metadata object to file
        """
        printtime(u'Printing metadata to file', self.start)
        import json
        # Iterate through each sample in the analysis
        for sample in self.metadata:
            # Set the name of the json file
            jsonfile = '{}/{}_metadata.json'.format(self.resultspath, sample.name)
            try:
                # Open the metadata file to write
                with open(jsonfile, 'wb') as metadatafile:
                    # Write the json dump of the object dump to the metadata file
                    json.dump(sample.dump(), metadatafile, sort_keys=True, indent=4, separators=(',', ': '))
            except IOError:
                pass

    def __init__(self, args, start):
        import multiprocessing
        from glob import glob
        from Queue import Queue
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.path = os.path.join(args.path, "")
        self.sequencepath = os.path.join(args.sequencepath, '') if args.sequencepath else '{}sequences/' \
            .format(self.path) if os.path.isdir('{}sequences'.format(self.path)) else self.path
        try:
            self.fastafile = glob('{}*.fa*'.format(self.sequencepath))[0]
        except IndexError:
            print u'Could not find a .fasta assembly file'
            quit()
        assert os.path.isdir(self.path), u'Supplied path location is not a valid directory {0!r:s}'.format(self.path)
        self.reportpath = os.path.join(self.path, 'reports')
        self.resultspath = os.path.join(self.path, 'results')
        make_path(self.reportpath)
        make_path(self.resultspath)
        self.start = start
        self.cpus = int(multiprocessing.cpu_count())
        # If an argument for the file is provided, use it
        self.coords = args.coordinateFile
        # If there is no path information present in the argument, use :path + the file name
        if '/' not in self.coords:
            self.coords = '{}{}'.format(self.path, self.coords)
            assert os.path.isfile(self.coords), \
                u'Could not find the contig/coordinate file. Please double check the supplied arguments'
        # Set the delimiter
        self.delimiter = args.delimiter.lower()
        if self.delimiter == 'space':
            self.delimiter = ' '
        elif self.delimiter == 'tab':
            self.delimiter = '\t'
        elif self.delimiter == 'comma' or self.delimiter == ',':
            self.delimiter = ','
        # Set the amplicon size
        self.ampliconsize = args.ampliconSize
        # Create variables for the analysis
        self.metadata = list()
        self.contigcoords = dict()
        self.primerqueue = Queue()
        self.devnull = open(os.devnull, 'wb')
        # Run the analyses
        self.snpprimers()

# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    from argparse import ArgumentParser
    from time import time
    # Parser for arguments
    parser = ArgumentParser(
        description=u'Extracts nucleotide sequence data surrounding SNP coordinates in .fasta assembly files.'
                    u'Create primers to amplify extracted regions using Primer3'
    )
    parser.add_argument('path',
                        help=u'Specify path of folder in which the analyses are to be performed')
    parser.add_argument('-c', '--coordinateFile',
                        required=True,
                        help=u'File with .txt extension containing one contig name and SNP coordinate pair separated '
                             u'by a delimiter (see below) per line. If a path is not specified in the file name, then '
                             u'the program will assume that the file is in the path supplied above.')
    parser.add_argument('-d', '--delimiter',
                        default='space',
                        help=u'The delimiter used to separate contig name and SNP coordinates. Popular options are '
                             u'"space", "tab", and "comma". Default is space. Be aware that a delimiter, '
                             u'such as "-" will break the program if there are hyphens in your contig names')
    parser.add_argument('-a', '--ampliconSize',
                        default=300,
                        help=u'The size of the amplicon to extract from the assembly file. The default is 300; '
                             u'150 bp upstream and 150 bp downstream of the SNP. I suppose the amplicon will actually '
                             u'be 301 bp, since the SNP position needs to be included, but round numbers are easier')
    parser.add_argument('-s', '--sequencepath',
                        help=u'The location of the .fasta assembly file. If not supplied, the program will assume that '
                             u'this file is in path/sequences')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Get the starting time for use in print statements
    starttime = time()
    # Run the function
    SNPprimers(arguments, starttime)
    printtime(u'SNP extraction and primer creation complete', starttime)
    # Print a bold, green exit statement
    print u'\033[92m' + u'\033[1m' + u'\nElapsed Time: %0.2f seconds' % (time() - starttime) + u'\033[0m'
