#!/usr/bin/env python
""" Determine the best MAF block (determined by score) that cover a specified
genome
"""
import Bio.AlignIO
import galaxy.maf_utilities
import optparse
import os
import pybedtools
import re
import seq_utils
import sys
import vamp.utils


__version__ = "0.2 BETA"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2011, Lance Parsons"
__license__ = ("BSD 2-Clause License http://www.opensource.org/licenses/BSD-2-Clause")


VAMP_HOME = vamp.utils.get_vamp_home()


def main(argv=None):
    if argv is None:
        argv = sys.argv

    usage = "Usage: %prog [options] maf_file"
    parser = optparse.OptionParser(usage=usage, version=('%prog version ' + globals()['__version__']), description=globals()['__doc__'])
    parser.add_option('-r', '--reference', help='Reference species')
    parser.add_option('-s', '--species', help='List of species to include')
    parser.add_option('-c', '--chromosome',
                      help='Sequence ID of the chromosome for '
                      'which to generate the alignment net (e.g. NC_001806)')
    parser.add_option('-o', '--output_dir', default=None,
                      help='Directory to store output file, default is maf file directory')
    parser.add_option('--consensus_sequence', action='store_true',
                      default=False,
                      help='Output "consensus sequence" for each species in '
                      'files named [species].[chromosome].consensus.fasta')
    parser.add_option('--reference_fasta', default=False,
                      help='Check MAF file against this fasta '
                      '(for troubleshooting, debugging)')
    parser.add_option('-v', '--verbose', action='store_true', default=False,
                      help='verbose output')
    try:
        (options, args) = parser.parse_args(argv[1:])
        if len(args) < 1:
            parser.error('Must specify MAF file')
        if options.reference is None:
            parser.error('Must specify reference species (-r option)')
        if options.species is None:
            parser.error('Must specify comma separated '
                         'list of species (-s option)')
        if options.chromosome is None:
            parser.error('Must specify the sequence id of the reference '
                         'species for which to generate the alignment '
                         'net (-c option)')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    maf_filename = args[0]
    (output_dir, output_basename) = os.path.split(maf_filename)
    (output_basename, unused_maf_ext) = os.path.splitext(output_basename)
    if options.output_dir:
        output_dir = options.output_dir
    # TODO Check for existence of directory (create if necessary)
    output_basepath = os.path.join(output_dir, output_basename)
    aligned_fasta_filename = "%s.net.afa" % output_basepath
    if options.reference_fasta:
        vamp.utils.verify_maf_fasta(maf_filename,
                                    options.reference,
                                    options.reference_fasta,
                                    options.verbose)

    if options.verbose:
        print "VAMP_HOME=%s" % VAMP_HOME

    # # Split comma separated string
    species = galaxy.maf_utilities.parse_species_option(options.species)
    species.append(options.reference)
    species = set(species)
    if options.verbose:
        print "Species: %s" % species

    ref_src = "%s.%s" % (options.reference, options.chromosome)

    net_alignment, interval_list = vamp.utils.get_sequence_net_alignment(maf_filename, options.reference, options.chromosome, species, options.verbose)

    if options.verbose:
        print net_alignment
    Bio.AlignIO.write(net_alignment, aligned_fasta_filename, "fasta")

    # Consensus sequence and contig tracking
    if options.consensus_sequence:
        interval_list = interval_list.sort()

        species_interval_lists = alignment_interval_to_species_intervals(net_alignment, interval_list, maf_filename, ref_src, options.verbose)
        for src in species_interval_lists.keys():
            interval_list = species_interval_lists[src]
            contig_composition_filename = '%s.consensus_contig_composition.gff' % (src)
            contig_composition_file = open(contig_composition_filename, 'wb')
            contig_composition_file.write(str(interval_list))
            contig_composition_file.close()

            contig_composition_summary_filename = '%s.consensus_contig_composition_summary.txt' % (src)
            contig_composition_summary_file = open(contig_composition_summary_filename, 'wb')
            composition_summary = vamp.utils.summarize_contig_composition(pybedtools.BedTool(contig_composition_filename), 'src_seq', 'src_seq_start', 'src_seq_end', 'src_strand', 'src_size')
            contig_composition_summary_file.write("seq\tstart\tend\tcontig\tcontig_start\tcontig_end\tcontig_strand\tcontig_size\n")
            for section in composition_summary:
                contig_composition_summary_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                                      (section['seq'],
                                                       section['start'],
                                                       section['end'],
                                                       section['contig'],
                                                       section['contig_start'],
                                                       section['contig_end'],
                                                       section['strand'],
                                                       section['contig_size']))
            contig_composition_summary_file.close()

        for src in net_alignment:
            consensus_filename = '%s.consensus.fasta' % (src.id)
            consensus_file = open(consensus_filename, 'wb')
            consensus_seq = re.sub('-', '', str(src.seq))
            consensus_file.write('>%s\n%s\n' % (src.id, consensus_seq))
            # print("%s" % sequence.id)
            consensus_file.close()

    return 0


def alignment_interval_to_species_intervals(net_alignment, interval_list, maf_filename, reference_src, verbose=False):
    for alignment in net_alignment:
        if alignment.id == reference_src:
            ref_seq = alignment.seq
    species_interval_lists = {}
    for alignment in net_alignment:
        species_interval_list = None
        (species, unused_contig) = alignment.id.split('.')
        for interval in interval_list:
            # Get block
            block = vamp.utils.get_block_by_label(maf_filename, interval.attrs['maf_block'])
            if verbose:
                print block.attributes['label']
            # Get reference component
            ref_components = galaxy.maf_utilities.get_components_by_src(block, reference_src)
            if len(ref_components) > 1:
                raise Exception("More than one reference component [%s] in maf block [%s]!" % (reference_src, block))
            elif len(ref_components) < 1:
                raise Exception("No reference component [%s] found in maf block [%s]!" % (reference_src, block))
            ref_component = ref_components[0]
            # Get component by src
            components = galaxy.maf_utilities.get_components_by_src_start(block, species)
            # Add attr for src, start, end where src = compoenent.src,
            # start = block_start with gaps removed + c.start
            # end = block_end with gaps removed + c.start
            src = None
            start = None
            end = None
            src_size = None
            src_strand = None
            # maf_block = interval['maf_block']
            #if len(components) < 1:
            #    raise Exception("No components found by in maf block [%s] for species [%s]!" % (block.attributes['label'], species))
            #elif len(components) > 1:
            #    raise Exception("Multiple components found by in maf block [%s] for species [%s]!" % (block.attributes['label'], species))
            if len(components) == 1:
                component = components[0]
                src = component.src
                # NOTE if interval is reversed that means that ???
                # if component is reversed, that means that ???
                if verbose:
                    print component.src
                    print "component strand %s" % component.strand
                    print "interval.strand %s" % interval.strand
                    print "component start %s" % component.start
                    print "component end %s" % component.end

                src_strand = '+'
                if ref_component.strand == '-':
                    if verbose:
                        print "Reference is minus strand for this block, reverse the component"
                    component = component.reverse_complement()
                if component.strand == '-' and component.src != reference_src:
                    if verbose:
                        print "Component is minus strand for this block, reverse the component"
                    src_strand = '-'
                    component = component.reverse_complement()
                    start, end = seq_utils.utils.convert_interval_gapped_to_nongapped(component.text, int(interval.attrs['block_start']), int(interval.attrs['block_end']))
                    # The start and end should be relative to the forward strand of the source
                    # Since the component was reversed, we should provide reverse positions
                    start = component.size - start
                    end = component.size - end
                else:
                    start, end = seq_utils.utils.convert_interval_gapped_to_nongapped(component.text, int(interval.attrs['block_start']), int(interval.attrs['block_end']))
                if verbose:
                    print "src start after block trimming: %s - %s" % (start, end)
                start = start + component.start
                end = end + component.start
                src_size = component.src_size
            interval.attrs['src_seq'] = src
            interval.attrs['src_seq_start'] = start
            interval.attrs['src_seq_end'] = end
            interval.attrs['src_strand'] = src_strand
            interval.attrs['src_size'] = src_size

            # include_end_gaps = not(component.src == reference_src)

            # aligned_start, aligned_end = seq_utils.utils.convert_interval_nongapped_to_gapped(str(ref_seq), interval.start, interval.end, include_end_gaps=include_end_gaps)
            aligned_start = interval.start
            aligned_end = interval.end
            if verbose:
                print "Ref Seq: %s" % str(ref_seq)
                print "\tinterval start: %s" % interval.start
                print "\taligned start: %s" % aligned_start
                print "\tinterval end : %s" % interval.end
                print "\taligned end: %s" % aligned_end
            new_start, new_end = seq_utils.utils.convert_interval_gapped_to_nongapped(str(alignment.seq), aligned_start, aligned_end)
            interval.start = new_start
            interval.end = new_end
            interval.chrom = str(alignment.id)
            if verbose:
                print interval
            new_bedtool = pybedtools.BedTool(str(interval), from_string=True)
            if verbose:
                print(str(new_bedtool)),
            if species_interval_list:
                species_interval_list = species_interval_list.cat(new_bedtool, postmerge=False)
            else:
                species_interval_list = new_bedtool
        species_interval_lists[alignment.id] = species_interval_list
    return species_interval_lists


if __name__ == '__main__':
    sys.exit(main())
