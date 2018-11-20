#!/usr/bin/env python
""" Compares genomes using multiple alignment
Input: GFF
Output: GTF File

Very simple and naieve GFF to GTF converter.
Writen to handle GFF output from BioPerls genbank2gff3.pl script for simple DNA viruses
"""

__version__ = "0.1"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2012, Lance Parsons"
__license__ = "BSD 2-Clause License http://www.opensource.org/licenses/BSD-2-Clause"

import optparse
import pybedtools
import re
import sys


def main(argv=None):
    if argv is None:
        argv = sys.argv
    usage = "Usage: %prog [options] gff_file"
    parser = optparse.OptionParser(usage=usage, version='%prog version ' + globals()['__version__'], description=globals()['__doc__'])
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
    try:
        (options, args) = parser.parse_args(argv[1:])
        if len(args) < 1:
            parser.error('You must specify a gff file.')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    gtf_features = ("CDS", "start_codon", "stop_codon", "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS", "exon")

    intervals = pybedtools.BedTool(args[0])
    if intervals.file_type != 'gff':
        parser.error('File specified does not appear to be a valid GFF file.')

    attrib_pattern = re.compile(r'^(.+)\t([a-z_]+) ([^;]+);([a-z_]+) ([^;]+)$')

    # For each interval, filter and generate new record for gtf file
    for feature in iter(intervals):
        if options.verbose:
            print "Feature: %s[%s] (%s:%s-%s)" % (feature.name, feature[2], feature.chrom, feature.start, feature.end)

        if feature[2] not in gtf_features:
            if options.verbose:
                print "Skipping, '%s' is not a valid GTF feature type" % feature[2]
            continue

        # Find pos in aligned reference
        if options.verbose:
            print "Feature: %s - %s" % (feature.start, feature.end)
        # Create new interval record (for GFF, create with no attributes an only add specific ones)
        new_interval = pybedtools.create_interval_from_list([feature[0],
                                                             feature[1],
                                                             feature[2],
                                                             feature[3],
                                                             feature[4],
                                                             feature[5],
                                                             feature[6],
                                                             feature[7],
                                                             'gene_id ""; transcript_id "";'])
        attr_list = feature.attrs.keys()
        gene_id_attrs = ['gene', 'gene_id']
        for attr in gene_id_attrs:
            if attr in attr_list:
                new_interval.attrs['gene_id'] = feature.attrs[attr]
        transcript_id_attrs = ['transcript', 'transcript_id']
        for attr in transcript_id_attrs:
            if attr in attr_list:
                new_interval.attrs['transcript_id'] = feature.attrs[attr]
        if feature[2] == 'CDS' and 'Parent' in attr_list:
            new_interval.attrs['transcript_id'] = feature.attrs['Parent']
        if new_interval.attrs['transcript_id'] == '':
            new_interval.attrs['transcript_id'] = '%s.t01' % new_interval.attrs['gene_id']

        gff_string = str(new_interval)
        m = attrib_pattern.search(gff_string.strip())
        # print m.groups()
        gff_string = "%s\t%s \"%s\"; %s \"%s\";" % (m.group(1), m.group(4), m.group(5), m.group(2), m.group(3))
        print gff_string

    return 0


if __name__ == '__main__':
    sys.exit(main())
