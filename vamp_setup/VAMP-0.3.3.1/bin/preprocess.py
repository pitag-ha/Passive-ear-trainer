#!/usr/bin/env python
""" Gather statistics from preprocessing
"""

__version__ = "0.2"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2011, Lance Parsons"
__license__ = "BSD 2-Clause"


import glob
import optparse
import os
import sys


def main(argv=None):
    if argv is None:
        argv = sys.argv

    usage = "Usage: %prog [options] preprocessing_directory basename"
    parser = optparse.OptionParser(usage=usage, version='%prog version ' + globals()['__version__'], description=globals()['__doc__'])
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
    try:
        (options, args) = parser.parse_args(argv[1:])  # @UnusedVariable
        if len(args) < 2:
            parser.error('Please specify the preprocessing directory and basename')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    print ("")
    print ("Preprocessing Directory: %s" % args[0])
    preprocessingStats = Preprocess(args[0], args[1])
    preprocessingStats.printStats()
    return 0


class Preprocess():
    "Parse Preprocessing Results"
    def __init__(self, directory, basename):
        self.directory = directory
        self.basename = basename

        self.originalStats = self.__getOriginalStats()
        self.cutadaptStats = self.__getCutadaptStats()
        self.artifactStats = self.__getArtifactStats()
        self.hostContaminantStats = self.__getHostContaminantStats()
        self.tqsStats = self.__getTqsStats()

    def printStats(self):
        print ("Step\tNumber of Fragments Read 1\tRead 2")
        if (len(self.originalStats) != 0):
            r1 = self.originalStats[1]['FastQC'].basicStats['Total Sequences']
            r2 = self.originalStats[2]['FastQC'].basicStats['Total Sequences']
        else:
            r1 = ""
            r2 = ""
        print ("Original\t%s\t%s" % (r1, r2))

        if (len(self.cutadaptStats) != 0):
            r1 = self.cutadaptStats[1]['FastQC'].basicStats['Total Sequences']
            r2 = self.cutadaptStats[2]['FastQC'].basicStats['Total Sequences']
        else:
            r1 = ""
            r2 = ""
        print ("Cutadapt\t%s\t%s" % (r1, r2))

        if (len(self.artifactStats) != 0):
            r1 = self.artifactStats[1]['FastQC'].basicStats['Total Sequences']
            r2 = self.artifactStats[2]['FastQC'].basicStats['Total Sequences']
        else:
            r1 = ""
            r2 = ""
        print ("Artifact Filtered\t%s\t%s" % (r1, r2))

        if (len(self.hostContaminantStats) != 0):
            r1 = self.hostContaminantStats[1]['FastQC'].basicStats['Total Sequences']
            r2 = self.hostContaminantStats[2]['FastQC'].basicStats['Total Sequences']
        else:
            r1 = ""
            r2 = ""
        print ("Host Filtered\t%s\t%s" % (r1, r2))

#        if (len(self.tqsStats) != 0):
#            r1 = self.tqsStats[1]['FastQC'].basicStats['Total Sequences']
#            r2 = self.tqsStats[2]['FastQC'].basicStats['Total Sequences']
#        else:
#            r1 = ""
#            r2 = ""
#        print ("TQS trimmed\t%s\t%s" % (r1, r2))

    def __getOriginalStats(self):
        stats = {}
        for read in xrange(1, 3):
            stats[read] = {}
            fastqcDir = os.path.join(self.directory, '%s_read%s_fastqc' % (self.basename, read))
            print "Original Fastqc: %s" % fastqcDir
            if (os.path.exists(fastqcDir)):
                stats[read]['FastQC'] = FastQC(fastqcDir)
        return stats

    def __getCutadaptStats(self):
        stats = {}
        for read in xrange(1, 3):
            stats[read] = {}
            fastqcDir = os.path.join(self.directory, '%s_cutadapt_read%s_fastqc' % (self.basename, read))
            if (os.path.exists(fastqcDir)):
                stats[read]['FastQC'] = FastQC(fastqcDir)
        return stats

    def __getArtifactStats(self):
        stats = {}
        for read in xrange(1, 3):
            stats[read] = {}
            fastqcDir = os.path.join(self.directory, '%s_cutadapt_artifact_filtered_read%s_fastqc' % (self.basename, read))
            if (os.path.exists(fastqcDir)):
                stats[read]['FastQC'] = FastQC(fastqcDir)
        return stats

    def __getHostContaminantStats(self):
        stats = {}
        for read in xrange(1, 3):
            stats[read] = {}
            fastqcDir = os.path.join(self.directory, '%s_cutadapt_artifact_filtered_host_filtered_%s_fastqc' % (self.basename, read))
            if (os.path.exists(fastqcDir)):
                stats[read]['FastQC'] = FastQC(fastqcDir)
        return stats

    def __getTqsStats(self):
        stats = {}
        for read in xrange(1, 3):
            stats[read] = {}
            fastqcDirs = glob.glob(os.path.join(self.directory, '%s_cutadapt_artifact_filtered_host_filtered_%s_T*C*_fastqc' % (self.basename, read)))
            if len(fastqcDirs) > 1:
                sys.stderr.write('Multiple TQS outputs found, using first one found: %s' % fastqcDirs[1])
            for fastqcDir in fastqcDirs:
                stats[read]['FastQC'] = FastQC(fastqcDir)
        return stats


class FastQC():
    "Parse FastQC results"
    def __init__(self, directory):
        self.directory = directory
        self.readStats()

    def readStats(self):
        # Read data file
        dataFilename = os.path.join(self.directory, "fastqc_data.txt")
        if (not os.path.exists(dataFilename)):
            raise FastqcDataError("Unable to find fastqc_data.txt file")
        dataFile = open(dataFilename, 'rb')

        # Get FastQC Version
        line = dataFile.readline()
        data = line.strip().split('\t')
        if (data[0] != '##FastQC'):
            raise FastqcDataError("Unable to determine version of FastQC")
        self.version = data[1]

        self.__getBasicStats(dataFile)

    def __getBasicStats(self, dataFile):
        line = dataFile.readline()
        data = line.strip().split('\t')
        if (data[0] != '>>Basic Statistics'):
            raise FastqcDataError("Unable to read Basic Statistics Module Results")
        basicStats = {}
        basicStats['Summary'] = data[1]
        for line in dataFile:
            data = line.strip().split('\t')
            if (data[0].startswith('#')):
                continue
            if (data[0] == '>>END_MODULE'):
                break
            basicStats[data[0]] = data[1]
        self.basicStats = basicStats


class FastqcDataError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

if __name__ == '__main__':
    sys.exit(main())
