#!/usr/bin/env python
"""Compares genomes using multiple alignment

Input: Aligned Fasta and GFF.

Output: For each non-reference record: gff, feature sequences, diff summary,
    vcf file.

Current limited to single sequence (chromosome) at a time.
"""
from operator import itemgetter
from seq_utils import convert_coordinates, summarize_alignments
import Bio.AlignIO
import Bio.Seq
import optparse
import os
import pybedtools
import re
import sys
import vamp.utils

__version__ = "0.1"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2012, Lance Parsons"
__license__ = ("BSD 2-Clause License "
               "http://www.opensource.org/licenses/BSD-2-Clause")


# Set VAMP_HOME for use in makefiles
VAMP_HOME = vamp.utils.get_vamp_home()


def main(argv=None):
    if argv is None:
        argv = sys.argv
    usage = "Usage: %prog [options] aligned_fasta gene_gff"
    parser = optparse.OptionParser(usage=usage,
                                   version='%prog version ' +
                                   globals()['__version__'],
                                   description=globals()['__doc__'])
    parser.add_option('-r', '--reference',
                      help='Sequence id of reference sequence in aligned fasta file')
    parser.add_option('--align_format', metavar='FORMAT', default='fasta', help='Alignment format (default: %default)')
    parser.add_option('-o', '--output', metavar='PREFIX', default='compare_genomes_output%s' % os.sep, help='Output prefix (default: %default)')
    # parser.add_option('--snpeff', metavar='PATH', default=None, help='Path to snpEff.jar (required for snpEff output)')
    parser.add_option('--gff_feature_types', default='CDS,exon,gene,mRNA,stem_loop', help='Comma separated list of gff feature types to parse (default: %default)')
    parser.add_option('--gff_attributes', default='ID,Parent,Note,gene,function,product', help='Comma separated list of feature attributes to carry over (default: %default)')
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
    try:
        (options, args) = parser.parse_args(argv[1:])
        if len(args) < 2:
            parser.error('You must specify an aligned fasta file and a gff file')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    if options.verbose:
        print "VAMP_HOME=%s" % VAMP_HOME

    feature_types = options.gff_feature_types.split(',')
    feature_attributes = options.gff_attributes.split(',')

    alignment_filename = args[0]
    alignment = Bio.AlignIO.read(alignment_filename, options.align_format)
    if options.verbose:
        print ("The alignment is %i bases long" % len(alignment[0, :]))

    intervals = pybedtools.BedTool(args[1])

    # Setup output files and get gap positions for each aligned sequence record
    gap_positions = {}
    gap_re = re.compile(r'-+')
    if not os.path.exists(options.output):
        os.makedirs(options.output)
    output_files = {}
    alignment_basename = os.path.basename(alignment_filename)
    output_files['alignment_positions'] = open('%s%s.%s' % (options.output, alignment_basename, intervals.file_type), 'wb')
    output_files['merged_vcf'] = open('%s%s_merged.vcf' % (options.output, alignment_basename), 'wb')
    output_files['merged_vcf'].write('##fileformat=VCFv4.0\n')
    output_files['merged_vcf'].write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    output_files['merged_vcf'].write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
    reference_found = False
    reference_index = None
    for i in xrange(0, len(alignment)):
        record = alignment[i]
        gap_positions[record.id] = []
        gap_iter = gap_re.finditer(str(record.seq))
        for gap in gap_iter:
            gap_positions[record.id].append(gap)
        if record.id != options.reference:
            output_files[record.id] = {}
            output_files[record.id]['dna_feature_gff'] = open('%s%s.%s' % (options.output, record.id, intervals.file_type), 'wb')
            output_files[record.id]['dna_feature_fasta'] = open('%s%s_features_dna.fasta' % (options.output, record.id), 'wb')
            output_files[record.id]['dna_diff_summary'] = open('%s%s_diff_summary.txt' % (options.output, record.id), 'wb')
            output_files[record.id]['dna_diff_summary'].write('Feature\tType\tMatches\tMismatches\tChange Events\tChange Event List\tMistmatch List\n')
            output_files[record.id]['vcf'] = open('%s%s.vcf' % (options.output, record.id), 'wb')
            output_files[record.id]['vcf'].write('##fileformat=VCFv4.0\n')
            output_files[record.id]['vcf'].write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            output_files[record.id]['vcf'].write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % record.id)
            output_files['merged_vcf'].write('\t%s' % record.id)
        else:
            reference_found = True
            reference_index = i
    output_files['merged_vcf'].write('\n')
    if not reference_found:
        sys.stderr.write("Reference sequence '%s' not found in fasta file '%s'" % (options.reference, alignment_filename))
        return 3

    # For each interval, convert coords, get diff summary
    for feature in iter(intervals):
        feature_type = None
        if intervals.file_type == 'gff':
            feature_type = feature[2]
        if (intervals.file_type == 'gff' and feature[2] == 'gene') or intervals.file_type == 'bed':
            feature_all_output = open('%s%s_all_dna.fasta' % (options.output, feature.name), 'wb')
        if (feature_type is None) or (feature_type in feature_types):
            if options.verbose:
                print "Feature: %s[%s] (%s:%s-%s)" % (feature.name, feature_type, feature.chrom, feature.start, feature.end)

            # Find pos in aligned reference
            reflength = len(alignment[reference_index].seq)
            numgaps = alignment[reference_index].seq.count('-')
            ungapped_reflength = reflength - numgaps
            align_ref_start = convert_coordinates.find_aligned_position(gap_positions[alignment[reference_index].id], feature.start)
            align_ref_end = convert_coordinates.find_aligned_position(gap_positions[alignment[reference_index].id], feature.end)
            if feature.start > ungapped_reflength or feature.end > ungapped_reflength:
                sys.stderr.write("Feature %s[%s] is truncated in the reference!  " % (feature.name, feature_type))
                sys.stderr.write("Position is %s:%i-%i but sequence is only %i bases long.  " %
                                  (feature.chrom, feature.start, feature.end,
                                  ungapped_reflength))
                align_ref_start = min(reflength, align_ref_start)
                align_ref_end = min(reflength, align_ref_end)
                sys.stderr.write("Using %s:%i-%i.\n" % (feature.chrom, align_ref_start, align_ref_end))
            assert(alignment[reference_index].seq[align_ref_start - 1] != '-')
            assert(alignment[reference_index].seq[align_ref_end - 1] != '-')
            if options.verbose:
                print "Aligned Feature: %s - %s" % (align_ref_start, align_ref_end)
            diff_summary = summarize_alignments.summary_of_alignment(alignment[:, align_ref_start:align_ref_end], options.reference)
            if options.verbose:
                print alignment

            # Create new interval record (for GFF, create with no attributes an only add specific ones)
            new_interval = feature
            if intervals.file_type == 'gff':
                new_interval = pybedtools.create_interval_from_list([feature[0],
                                                                     feature[1],
                                                                     feature[2],
                                                                     feature[3],
                                                                     feature[4],
                                                                     feature[5],
                                                                     feature[6],
                                                                     feature[7],
                                                                     "ID="])
                new_attrs = dict([(i, feature.attrs[i]) for i in feature_attributes if i in feature.attrs])
                for k in new_attrs.keys():
                    # print k
                    # print new_attrs[k]
                    new_interval.attrs[k] = new_attrs[k]
                if options.verbose:
                    print new_interval.attrs

            alignment_position_interval = feature
            alignment_position_interval.start = align_ref_start
            alignment_position_interval.end = align_ref_end
            if intervals.file_type == 'gff':
                alignment_position_interval[1] = 'compare_genomes'
            output_files['alignment_positions'].write(str(alignment_position_interval))

            # For each alignment record, get new attribute locations and output fasta and gff files
            for record in alignment[0:]:
                num_new_gaps_start = str(record[:align_ref_start].seq).count('-')
                num_new_gaps_end = str(record[:align_ref_end].seq).count('-')
                new_start = feature.start - num_new_gaps_start
                new_end = feature.end - num_new_gaps_end
                new_interval.chrom = record.id
                new_interval.start = new_start
                new_interval.end = new_end
                feature_seqrecord_aligned = record[align_ref_start:align_ref_end]
                new_seqrecord = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(str(feature_seqrecord_aligned.seq).replace('-', ''),
                                                                    feature_seqrecord_aligned.seq.alphabet))
                new_seqrecord.id = feature.name
                new_seqrecord.description = '%s %s:%i-%i' % (feature_type, new_interval.chrom, (new_interval.start + 1), new_interval.end)

                if record.id != alignment[reference_index].id:
                    output_files[record.id]['dna_feature_fasta'].write(new_seqrecord.format("fasta"))

                    output_files[record.id]['dna_feature_gff'].write(str(new_interval))
                    output_files[record.id]['dna_diff_summary'].write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                                                                      (feature.name,
                                                                       feature_type,
                                                                       diff_summary[record.id]['match_count'],
                                                                       diff_summary[record.id]['mismatch_count'],
                                                                       diff_summary[record.id]['noncontig_count'],
                                                                       diff_summary[record.id]['noncontigs_string'],
                                                                       diff_summary[record.id]['mismatches_string']))

                if (intervals.file_type == 'gff' and feature[2] == 'gene') or intervals.file_type == 'bed':
                    new_seqrecord.id = '%s.%s' % (record.id, feature.name)
                    # new_seqrecord.description = '%s (%s)' % (record.id, new_seqrecord.description)
                    feature_all_output = open('%s%s_all_dna.fasta' % (options.output, feature.name), 'ab')
                    feature_all_output.write(new_seqrecord.format("fasta"))
                    feature_all_output.close()
        else:
            if options.verbose:
                print 'Skipping %s since feature_type is "%s"' % (feature.name, feature_type)

    diff_summary = summarize_alignments.summary_of_alignment(alignment, options.reference)
    merged_vcf_records = dict()
    for record in alignment[0:]:
        if record.id != alignment[reference_index].id:
            for variant in diff_summary[record.id]['noncontig_list']:
                # print variant
                # raw_input("Press Enter")
                vcf_key = (
                    alignment[reference_index].id,  # CHROM
                    variant['pos'],  # POS
                    '.',  # ID
                    re.sub(r'-+', r'', variant['ref']),  # REF
                    re.sub(r'-+', r'', variant['alt']),  # ALT
                    '.',  # QUAL
                    '.',  # FILTER
                    '.',  # INFO
                    'GT',  # FORMAT
                )
                vcf_record = merged_vcf_records.get(vcf_key, dict())
                vcf_record.setdefault(record.id, '1')
                merged_vcf_records[vcf_key] = vcf_record
                vcf_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
                    alignment[reference_index].id,  # CHROM
                    variant['pos'],  # POS
                    '.',  # ID
                    re.sub(r'-+', r'', variant['ref']),  # REF
                    re.sub(r'-+', r'', variant['alt']),  # ALT
                    '.',  # QUAL
                    '.',  # FILTER
                    '.',  # INFO
                    'GT',  # FORMAT
                    '1',  # GENOTYPE
                )
                output_files[record.id]['vcf'].write('%s\n' % vcf_line)

    sorted_keys = sorted(merged_vcf_records.keys(), key=itemgetter(0, 1))
    # sorted_keys = sorted(sorted_keys, key=lambda vcf_keys: vcf_key[0])
    for vcf_key in sorted_keys:
        vcf_line = '\t'.join(map(str, vcf_key))
        for i in xrange(0, len(alignment)):
            record = alignment[i]
            if record.id != options.reference:
                vcf_line = '%s\t%s' % (vcf_line, merged_vcf_records[vcf_key].get(record.id, 0))
        output_files['merged_vcf'].write('%s\n' % vcf_line)
    return 0


if __name__ == '__main__':
    sys.exit(main())
