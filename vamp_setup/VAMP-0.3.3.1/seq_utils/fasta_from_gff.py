#!/usr/bin/env python
""" Extract fasta sequences from regions defined in GFF/BED file and output fasta to stdout
"""

__version__ = "0.1"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2012, Lance Parsons, Lewis-Sigler Institute for Integrative Genomics, Princeton University"
__license__ = "BSD 2-Clause License http://www.opensource.org/licenses/BSD-2-Clause"

import Bio.SeqIO
import optparse
import pybedtools
import sys


def main(argv=None):
    if argv is None:
        argv = sys.argv

    usage = "Usage: %prog [options] sequence_file (fasta/fastq) feature_file (gff/bed)"
    parser = optparse.OptionParser(usage=usage, version='%prog version ' + globals()['__version__'], description=globals()['__doc__'])
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
    parser.add_option('--sequenceformat', default='fasta', help='Sequence file format (default: %default)')
    parser.add_option('--featureformat', default='gff', help='Sequence feature file format (default: %default)')
    try:
        (options, args) = parser.parse_args(argv[1:])
        if len(args) < 1:
            parser.error('argument missing')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    intervals = pybedtools.BedTool(args[1])
    seq_dict = Bio.SeqIO.index(args[0], options.sequenceformat)
    for interval in intervals:
        seq_record = seq_dict[interval.chrom]
        subseq = seq_record.seq[interval.start:interval.end]
        print('>%s\n%s\n' % (interval.name, subseq)),
    return 0


if __name__ == '__main__':
    sys.exit(main())
