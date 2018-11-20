#!/usr/bin/env python
"""
Convert fastq file to fasta
"""

__version__ = "0.1"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2011, Lance Parsons"
__license__ = "BSD 2-Clause License http://www.opensource.org/licenses/BSD-2-Clause"

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import optparse
import sys


def main(argv=None):
    if argv is None:
        argv = sys.argv

    usage = "Usage: %prog [options] fastq_file fasta_file"
    parser = optparse.OptionParser(usage=usage, version='%prog version ' + globals()['__version__'], description=globals()['__doc__'])
    parser.add_option('-w', '--wrap', type='int', default=0, help='Maximum length of lines, 0 means do not wrap (default: %default)')
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
    try:
        (options, args) = parser.parse_args(argv[1:])
        if len(args) < 1:
            parser.error('Please specify a fastq file as input')
        if len(args) < 2:
            parser.error('Please specify a fasta file as output')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    output_fh = open(args[1], 'wb')
    count = 0
    import fileinput
    fi = fileinput.FileInput(args[0], openhook=fileinput.hook_compressed)
    for title, seq, qual in FastqGeneralIterator(fi):  # @UnusedVariable
        count += 1
        output_fh.write(">%s\n" % title)
        if options.wrap > 0:
            # Do line wrapping
            for i in range(0, len(seq), options.wrap):
                output_fh.write(seq[i:i + options.wrap] + "\n")
        else:
            output_fh.write("%s\n" % seq)
    if options.verbose:
        print 'Converted %s records' % count
    return 0

if __name__ == '__main__':
    sys.exit(main())
