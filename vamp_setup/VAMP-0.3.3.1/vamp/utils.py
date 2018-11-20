from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from galaxy import maf_utilities
import Bio.Align
import Bio.SeqIO
import bx.align.maf
import copy
import os
import pybedtools
import re
import seq_utils.utils
import sys


def default_vamp_config():
    default_config_dir = os.path.join(get_vamp_home(), 'makefiles')
    default_config = None
    config = os.path.join(default_config_dir, 'config.mk')
    if os.path.exists(config):
        default_config = config
    else:
        config = os.path.join(default_config_dir, 'config.mk.template')
        if os.path.exists(config):
            default_config = config
    return default_config


def parse_makefile(filename):
    variable_rx = re.compile(r'(?P<variable>[a-zA-Z][a-zA-Z0-9_]+)\s*[:?]{0,1}=\s*(?P<value>.*)')
    variables = {}
    fh = open(filename, 'rb')
    for line in fh:
        line = line.strip()
        m = re.match(variable_rx, line)
        if m:
            variables[m.group('variable')] = m.group('value')
    return variables


def get_vamp_home():
    VAMP_HOME = os.path.abspath(os.path.join(os.path.dirname(sys.modules[__name__].__file__), '..'))
    return VAMP_HOME


def get_sequence_net_alignment(maf_filename, reference_species, sequence_id, species, verbose=False):
    '''
    Return alignment created by stitching MAF blocks along entire sequence (including gaps)
    Also returns a list of intervals relative to the alignment that indicate the MAF block,
    block start, and block end of the source of that piece of the alignment
    '''

    maf_file = open(maf_filename, 'rU')

    # Get reference.sequence_id length
    ref_seq_length = get_sequence_length_from_maf(maf_file, reference_species, sequence_id)
    if ref_seq_length is None:
        raise Exception('Source "%s.%s" not found in maf file: "%s"' % (reference_species, sequence_id, maf_filename))
    if verbose:
        print "Size of %s.%s: %i bp" % (reference_species, sequence_id, ref_seq_length)

    # Create alignment with full sequence from fasta for reference.sequence_id, use N for all
    # other species in alignment
    sequence_list = []
    for spec in species:
        seq = Seq('N' * ref_seq_length, generic_dna)
        # print "SPECIES %s" % spec
        # print "SEQ %s" % seq
        src = "%s.%s" % (spec, sequence_id)
        sequence_list.append(SeqRecord(seq, id=src))
    alignment = Bio.Align.MultipleSeqAlignment(sequence_list)
    if verbose:
        print alignment

    # Get all blocks with this species.sequence_id
    # Sort block by score, lowest first
    # index maf for use here, need to remove index_file when finished
    index, index_filename = maf_utilities.open_or_build_maf_index(
        maf_file.name, '%s.idx' % maf_file.name, species=[reference_species])
    if index is None:
        sys.stderr.write("Your MAF file appears to be malformed.\n")

    # Order blocks overlaping this position by score, lowest first
    blocks = []
    source = '%s.%s' % (reference_species, sequence_id)
    if verbose:
        print "Reordering MAF Blocks"
    for block, idx, offset in index.get_as_iterator_with_index_and_offset(source, 0, ref_seq_length):
        score = float(block.score)
        for i in range(0, len(blocks)):
            if score < blocks[i][0]:
                blocks.insert(i, (score, idx, offset))
                if verbose:
                    print "Adding block %s" % block.attributes['label']
                break
        else:
            blocks.append((score, idx, offset))
            if verbose:
                print "Adding block %s" % block.attributes['label']
    # print "REORDERED MAF BLOCKS"

    # Replace section of alignment with alignment from the MAF block
    interval_list = []
    for block_dict in blocks:
        block = block_dict[1].get_at_offset(block_dict[2])
        if verbose:
            print "Analyzing maf block %s" % block.attributes['label']
        for c in block.components:
            # TODO Refactor to function and only allow single period
            (spec, dummy_seq_id) = c.src.split('.',1)
            if verbose:
                print "    Component src %s" % c.src
            if spec == reference_species:
                # Replace alignment and calculate difference in length
                original_alignment_length = len(alignment[0])
                alignment, new_interval = replace_alignment_with_block(alignment, block, reference_species, sequence_id, verbose)
                interval_adjustment = len(alignment[0]) - original_alignment_length
                if verbose:
                    print "interval_adjustment = %s" % interval_adjustment
                # Generate new interval and add to list
                # new_bedtool = pybedtools.BedTool(str(new_interval), from_string=True)
                # if len(interval_list) > 0:
                #    interval_list = interval_list.cat(new_bedtool, postmerge=False)
                # else:
                #    interval_list = new_bedtool
                # For each interval, subtract the new interval and adjust block start, end if needed
                # adjusted_new_interval = interval_list[len(interval_list) - 1]
                interval_list_string = ''
                for interval in interval_list:  # [0:len(interval_list) - 1]:
                    trimmed_intervals = subtract_intervals(interval, new_interval)
                    for trimmed_interval in trimmed_intervals:
                        block_start_adjustment = trimmed_interval.start - interval.start
                        block_end_adjustment = trimmed_interval.end - interval.end
                        block_updated_interval = pybedtools.create_interval_from_list([
                            trimmed_interval.chrom,
                            "maf_net",
                            "alignment_block",
                            str(trimmed_interval.start + 1),
                            str(trimmed_interval.end),
                            ".",
                            '+',
                            ".",
                            "maf_block=%s;block_start=%i;block_end=%i" %
                            (trimmed_interval.attrs['maf_block'],
                             long(trimmed_interval.attrs['block_start']) + block_start_adjustment,
                             long(trimmed_interval.attrs['block_end']) + block_end_adjustment)
                        ])
                        interval_list_string += str(block_updated_interval) + "\n"
                interval_list_string += str(new_interval) + "\n"
                interval_list = pybedtools.BedTool(str(interval_list_string), from_string=True)
                # For each interval, update start and end by interval adjustment to account for gaps
                interval_list_string = ''
                for interval in interval_list:
                    if interval.start >= new_interval.end:
                        interval.start = interval.start + interval_adjustment
                        # interval.attrs['block_start'] = long(interval.attrs['block_start']) + interval_adjustment
                    if interval.end >= new_interval.end:
                        interval.end = interval.end + interval_adjustment
                        # interval.attrs['block_end'] = long(interval.attrs['block_end']) + interval_adjustment
                    interval_list_string += str(interval) + "\n"
                interval_list = pybedtools.BedTool(str(interval_list_string), from_string=True)

                if verbose:
                    for seq in alignment:
                        print "%s: Number of bases %i" % (seq.id, len(re.sub('-', '', str(seq.seq))))

    # Return alignment
    os.unlink(index_filename)
    maf_file.close()
    return alignment, interval_list


def subtract_intervals(interval1, interval2):
    '''
    Subtract two intervals, return list of resulting intervals
    '''
    results = [copy.deepcopy(interval1)]
    if interval2.end <= interval1.end and interval2.end > interval1.start:
        if interval2.start >= interval1.start and interval2.start < interval1.end:
            results.append(copy.deepcopy(interval1))
            results[0].end = interval2.start
            results[1].start = interval2.end
        else:
            results[0].start = interval2.end
    elif interval2.start >= interval1.start and interval2.start < interval1.end:
        results[0].end = interval2.start
    return results


def replace_alignment_with_block(alignment, block, reference_species, sequence_id, verbose=False):
    '''
    Return updated alignment, interval
    '''
    src = '%s.%s' % (reference_species, sequence_id)
    ref_component = block.get_component_by_src(src)
    if verbose:
        print "    Reference component %s" % ref_component
    # Orient block so reference is forward strand
    # reverse = False
    if ref_component.strand == '-':
        block = block.reverse_complement()
        ref_component = block.get_component_by_src(src)
        # reverse = True
        if verbose:
            print "  Reference in reverse orientation"
    ref_component_length_withgaps = len(ref_component.text)
    ref_align_seq = seq_utils.utils.get_seq_by_id_from_list(alignment, src)
    if verbose:
        print "Ref align seq %s" % ref_align_seq
    if ref_align_seq is None:
        return alignment
    find_start = ref_component.start
    find_end = ref_component.start + ref_component.size
    if verbose:
        print (str(ref_align_seq.seq), find_start, find_end)
    alignment_replace_start, alignment_replace_end = seq_utils.utils.convert_interval_nongapped_to_gapped(str(ref_align_seq.seq), find_start, find_end)
    if verbose:
        print '%s - start: %i - end: %i' % (ref_component.src, find_start, find_end)
        print '%s - length (inc gaps): %i' % (ref_component.src, ref_component_length_withgaps)
        print '%s - length (no gaps): %i' % (ref_component.src, ref_component.size)
        print '%s - Replace Start: %i - End: %i' % (ref_component.src, alignment_replace_start, alignment_replace_end)
        print '%s - Size of text to be replaced: %i' % (ref_component.src, (alignment_replace_end - alignment_replace_start))
        print 'Length of alignment before replacement: %i' % len(alignment[0])

    # Replace Sequences one by and one and construct new alignment
    sequence_list = []
    for sequence in alignment:
        (seq_spec, dummy_seq_id) = sequence.id.split('.',1)
        # print "seq_spec: %s" % seq_spec
        replacement_text = ('N' * (ref_component_length_withgaps))
        for c in block.components:
            (c_spec, dummy_c_id) = c.src.split('.',1)
            # print "c_spec: %s" % c_spec
            if c_spec == seq_spec:
                # print "MATCH"
                # if c.strand == '-':
                    # print "REVERSE"
                #    c = c.reverse_complement()
                replacement_text = c.text
                if verbose:
                    print "%s: Replacement text: %s" % (c.src, replacement_text)
                    print "%s: Length of block replacement text: %i" % (c.src, len(replacement_text))
                    print '%s: Number of bases in replacement text: %i' % (c.src, len(re.sub('-', '', replacement_text)))
        new_sequence_text = str(sequence.seq)
        if verbose:
            print '%s: Number of bases to be replaced: %i' % (sequence.id, len(re.sub('-', '', new_sequence_text[alignment_replace_start:alignment_replace_end])))
        new_sequence = Seq(new_sequence_text[0:alignment_replace_start] + replacement_text + new_sequence_text[alignment_replace_end:], generic_dna)
        sequence_list.append(SeqRecord(new_sequence, id=sequence.id))
    new_alignment = Bio.Align.MultipleSeqAlignment(sequence_list)
    block_start = 0 + (find_start - find_start)
    block_end = long(len(replacement_text)) - (find_end - find_end)
    # orientation = '+'
    # if reverse:
    #    orientation = '-'
    new_interval = pybedtools.create_interval_from_list([
        sequence_id,
        "maf_net",
        "alignment_block",
        str(alignment_replace_start + 1),
        str(alignment_replace_end),
        ".",
        '+',
        ".",
        "maf_block=%s;block_start=%i;block_end=%i" % (block.attributes['label'], block_start, block_end)
    ])
    # print new_interval

    if verbose:
        print 'Length of alignment after replacement: %i' % len(new_alignment[0])
        print block.attributes['label']
        # Bio.AlignIO.write(new_alignment, 'tmp%s.afa' % block.attributes['label'], "fasta")
    return new_alignment, new_interval


def get_sequence_length_from_maf(maf_file, reference_species, sequence_id):
    '''
    Return length of the reference_species.sequence_id
    '''
    for block in bx.align.maf.Reader(maf_file):
        for component in block.components:
            if component.src == '%s.%s' % (reference_species, sequence_id):
                return component.get_src_size()
    return None


def verify_maf_fasta(maf_filename, reference_species, fasta_filename, verbose=False):
    maf_file = open(maf_filename, 'rU')
    fasta_file = open(fasta_filename, 'rU')
    record_dict = Bio.SeqIO.to_dict(Bio.SeqIO.parse(fasta_file, 'fasta'))
    if verbose:
        print record_dict.keys()
    for block in bx.align.maf.Reader(maf_file):
        for component in block.components:
            (species, seq_id) = component.src.split('.',1)
            if species == reference_species:
                if component.strand == '-':
                    block = block.reverse_complement()
                    component = block.get_component_by_src(component.src)
                start = component.start
                end = component.start + component.size
                if verbose:
                    print 'Checking block %s: %s:%i-%i' % (block.attributes['label'], seq_id, start, end)
                fasta_sequence = str(record_dict[seq_id][start:end].seq)
                maf_sequence = re.sub('-', '', component.text)
                if fasta_sequence != maf_sequence:
                    print "No Match!"
                    print ">FASTA.%s\n%s" % (seq_id, fasta_sequence)
                    print ">MAF.%s\n%s" % (seq_id, maf_sequence)
    fasta_file.close()
    maf_file.close()


def get_block_by_label(maf_filename, label):
    selected_block = None
    maf_file = open(maf_filename, 'rU')
    for block in bx.align.maf.Reader(maf_file):
        if block.attributes['label'] == label:
            selected_block = block
            break
    return selected_block


def summarize_contig_composition(interval_list, src_tag, start_tag, end_tag, strand_tag, source_size_tag):
    '''
    Summarize the contig composition in a list of tuples
    (seq, start, end, contig, contig_start, contig_end, strand, contig_size)
    '''
    section_list = []
    for interval in interval_list:
        strand = interval.attrs[strand_tag]
        start = interval.attrs[start_tag]
        end = interval.attrs[end_tag]
        section = {'seq': interval.chrom,
                   'start': interval.start,
                   'end': interval.end,
                   'contig': interval.attrs[src_tag],
                   'contig_start': start,
                   'contig_end': end,
                   'strand': strand,
                   'contig_size': interval.attrs[source_size_tag]}
        section_list.append(section)
    return section_list


def read_contig_composition_summary(filename):
    '''
    Read contig composition summary and return attributes
    '''
    linenum = 0
    fh = open(filename, 'rb')
    for line in fh:
        linenum += 1
        if (line[0] != '#') and (linenum != 1):
            line = line.strip()
            fields = line.split('\t')
            section = ContigComposition({'seq': fields[0],
                                         'start': fields[1],
                                         'end': fields[2],
                                         'contig': fields[3],
                                         'contig_start': fields[4],
                                         'contig_end': fields[5],
                                         'strand': fields[6],
                                         'contig_size': fields[7]})
            yield section


def find_deletions(contig_composition_list, verbose=False):
    '''
    Find contigs with deletions and return tuple containing
    list of idicies of contigs to be replaced along with replacement. e.g.::

        ([2, 3],
         {'seq': 'chr', 'start': 3, 'end': 5, 'contig': 'contig1',
          'contig_start': 5, 'contig_end': 10, 'strand': '+', 'contig_size': 20})

    replaces contig_composition_list[2:3] with the new contig composition specified
    '''
    previous_contigs = {}
    replacements = []
    for idx, section in enumerate(contig_composition_list):
        if verbose:
            section.to_tab()
        if section['contig'] in previous_contigs:
            if section['strand'] == previous_contigs[section['contig']]['strand']:
                # Append to deletions and insertions
                if section['strand'] == '+':
                    start = previous_contigs[section['contig']]['end']
                    end = section['start']
                    contig_start = previous_contigs[section['contig']]['contig_end']
                    contig_end = section['contig_start']
                    contig_size = section['contig_size']
                elif section['strand'] == '-':
                    start = previous_contigs[section['contig']]['end']
                    end = section['start']
                    contig_start = previous_contigs[section['contig']]['contig_end']
                    contig_end = section['contig_start']
                    contig_size = section['contig_size']
                if contig_start != contig_end:
                    replacement = ContigComposition({'seq': section['seq'],
                                   'start': start,
                                   'end': end,
                                   'contig': section['contig'],
                                   'contig_start': contig_start,
                                   'contig_end': contig_end,
                                   'strand': section['strand'],
                                   'contig_size': contig_size})
                    to_be_replaced = range(previous_contigs[section['contig']]['last_index'] + 1, idx)
                    replacements.append((to_be_replaced, replacement))
        previous_contigs[section['contig']] = section
        previous_contigs[section['contig']]['last_index'] = idx
    return (replacements)


def update_sequence_with_replacements(seq, replacements, replacement_seq_dict):
    '''
    Update Seq object with replacements
    Replacements must be non-overlapping and sorted
    '''
    offset = 0
    updated_sequence = seq
    for replacement in replacements:
        (dummy_species, contig) = replacement['contig'].split('.',1)
        if replacement['strand'] == '+':
            replacement_seq = replacement_seq_dict[contig][int(replacement['contig_start']):int(replacement['contig_end'])]
        else:
            rev_start = int(replacement['contig_size'] - int(replacement['contig_end']))
            rev_end = int(replacement['contig_size'] - int(replacement['contig_start']))
            replacement_seq = replacement_seq_dict[contig].reverse_complement()[rev_start:rev_end]
        start = int(replacement['start']) + offset
        end = int(replacement['end']) + offset
        updated_sequence = updated_sequence[0:start] + replacement_seq + updated_sequence[end:len(seq)]
        updated_sequence.id = seq.id
        original_length = int(replacement['end']) - int(replacement['start'])
        replacement_length = int(replacement['contig_end']) - int(replacement['contig_start'])
        offset = offset + (replacement_length - original_length)
    return updated_sequence


def update_contig_composition_summary(contig_composition_summary, replacements):
    '''
    Update list of ContigComposition objects with replacements
    Replacements are a list of tuples containing a list of indicies of contigs to be replaced along with replacements
    The replacements must be non-overlapping and sorted
    '''
    last_cc_idx = 0
    updated_summary = []
    offset = 0
    for (replacement_indicies, replacement) in replacements:
        if len(replacement_indicies) >= 1:
            for idx in xrange(last_cc_idx, replacement_indicies[0]):
                updated_contig_composition = copy.deepcopy(contig_composition_summary[idx])
                updated_contig_composition['start'] += offset
                updated_contig_composition['end'] += offset
                updated_summary.append(updated_contig_composition)
            last_cc_idx = replacement_indicies[-1] + 1
        else:
            for idx in xrange(last_cc_idx, len(contig_composition_summary)):
                if contig_composition_summary[idx]['start'] < replacement['start']:
                    updated_contig_composition = copy.deepcopy(contig_composition_summary[idx])
                    updated_contig_composition['start'] += offset
                    updated_contig_composition['end'] += offset
                    updated_summary.append(updated_contig_composition)
                else:
                    last_cc_idx += 1
                    break
        updated_contig_composition = copy.deepcopy(replacement)
        updated_contig_composition['start'] += offset

        original_length = int(replacement['end']) - int(replacement['start'])
        replacement_length = int(replacement['contig_end']) - int(replacement['contig_start'])
        updated_contig_composition['end'] = updated_contig_composition['start'] + replacement_length
        updated_summary.append(updated_contig_composition)
        offset += replacement_length - original_length

    for idx in xrange(last_cc_idx, len(contig_composition_summary)):
        updated_contig_composition = copy.deepcopy(contig_composition_summary[idx])
        updated_contig_composition['start'] += offset
        updated_contig_composition['end'] += offset
        updated_summary.append(updated_contig_composition)
    return updated_summary


class ContigComposition(dict):
    '''Association of two genomic intervals, used to represent the composition
    of one interval by another'''
    def to_tab(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.get("seq"),
                                                  self.get("start"),
                                                  self.get('end'),
                                                  self.get('contig'),
                                                  self.get('contig_start'),
                                                  self.get('contig_end'),
                                                  self.get('strand'),
                                                  self.get('contig_size'))

    @staticmethod
    def tab_headings():
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("seq",
                                                  "start",
                                                  'end',
                                                  'contig',
                                                  'contig_start',
                                                  'contig_end',
                                                  'strand',
                                                  'contig_size')
