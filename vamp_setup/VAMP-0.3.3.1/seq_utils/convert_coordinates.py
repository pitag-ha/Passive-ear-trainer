#!/usr/bin/env python
""" Convert coordinates from GFF or BED file using multi-fasta alignments
"""

__version__ = "0.2"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2011, Lance Parsons"
__license__ = "BSD 2-Clause License http://www.opensource.org/licenses/BSD-2-Clause"

import Bio.AlignIO
import optparse
import pybedtools
import re
import summarize_alignments
import sys


def main(argv=None):
    if argv is None:
        argv = sys.argv

    usage = "Usage: %prog [options] aligned_fasta gff"
    parser = optparse.OptionParser(usage=usage, version='%prog version ' + globals()['__version__'], description=globals()['__doc__'])
    parser.add_option('--align_format', metavar='FORMAT', default='fasta', help='Alignment format (default: %default)')
    parser.add_option('--diff_summary', default=False, action='store_true', help='Output difference summary for each feature')
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
    try:
        (options, args) = parser.parse_args(argv[1:])
        if len(args) < 2:
            parser.error('argument missing')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    intervals = pybedtools.BedTool(args[1])
    alignment = Bio.AlignIO.read(args[0], options.align_format)
    if options.verbose:
        print ("The alignment is %i bases long" % len(alignment[0, :]))

    if options.diff_summary:
        diff_summary_output_files = {}
    feature_output_files = {}

    gap_positions = {}
    gap_re = re.compile(r'-+')
    gap_iter = gap_re.finditer(str(alignment[0].seq))
    gap_positions[alignment[0].id] = []
    for gap in gap_iter:
        gap_positions[alignment[0].id].append(gap)

    for record in alignment[1:]:
        feature_output_files[record.id] = open('%s.%s' % (record.id, intervals.file_type), 'wb')
        if options.diff_summary:
            diff_summary_output_files[record.id] = open('%s_diff_summary.txt' % record.id, 'wb')
            diff_summary_output_files[record.id].write('Feature\tMatches\tMismatches\tChange Events\tChange Event List\tMistmatch List\n')
        gap_positions[record.id] = []
        gap_iter = gap_re.finditer(str(record.seq))
        for gap in gap_iter:
            gap_positions[record.id].append(gap)
            #print gap.span()

    for feature in iter(intervals):
        if options.verbose:
            print "Feature: %s" % feature.name

        # Find pos in aligned reference
        if options.verbose:
            print "Feature: %s - %s" % (feature.start, feature.end)
        align_ref_start = find_aligned_position(gap_positions[alignment[0].id], feature.start)
        assert(alignment[0].seq[align_ref_start - 1] != '-')
        align_ref_end = find_aligned_position(gap_positions[alignment[0].id], feature.end)
        assert(alignment[0].seq[align_ref_end - 1] != '-')
        if options.verbose:
            print "Aligned Feature: %s - %s" % (align_ref_start, align_ref_end)
        if options.diff_summary:
            diff_summary = summarize_alignments.summary_of_alignment(alignment[:, align_ref_start:align_ref_end])
        num_ref_gaps_start = str(alignment[0, :align_ref_start].seq).count('-')
        num_ref_gaps_end = str(alignment[0, :align_ref_end].seq).count('-')
        for record in alignment[1:]:
            num_new_gaps_start = str(record[:align_ref_start].seq).count('-')
            num_new_gaps_end = str(record[:align_ref_end].seq).count('-')
            new_start = feature.start + num_ref_gaps_start - num_new_gaps_start
            new_end = feature.end + num_ref_gaps_end - num_new_gaps_end
            new_interval = feature
            new_interval.chrom = record.id
            new_interval.start = new_start
            new_interval.end = new_end
            feature_output_files[record.id].write(str(new_interval) + '\n')
            if options.diff_summary:
                diff_summary_output_files[record.id].write('%s\t%s\t%s\t%s\t%s\t%s\n' %
                                                           (feature.name,
                                                            diff_summary[record.id]['match_count'],
                                                            diff_summary[record.id]['mismatch_count'],
                                                            diff_summary[record.id]['noncontig_count'],
                                                            diff_summary[record.id]['noncontigs_string'],
                                                            diff_summary[record.id]['mismatches_string'])
                                                           )

    return 0


def find_aligned_position(gap_positions, pos):
    align_pos = pos
    for gap in gap_positions:
        if gap.start() <= align_pos:
            align_pos = align_pos + (gap.end() - gap.start())
        else:
            continue
    return align_pos

if __name__ == '__main__':
    sys.exit(main())
