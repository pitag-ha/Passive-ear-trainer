#!/usr/bin/env python
""" Find contigs with deletions in contig composition summary.
"""

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from vamp.utils import ContigComposition
import Bio.AlignIO
import Bio.SeqIO
import optparse
import os
import re
import sys
import vamp.utils


__version__ = "0.1 BETA"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2011, Lance Parsons"
__license__ = ("BSD 2-Clause License http://www.opensource.org/licenses/BSD-2-Clause")


VAMP_HOME = vamp.utils.get_vamp_home()


def main(argv=None):
    if argv is None:
        argv = sys.argv

    usage = "Usage: %prog [options] contig_composition.txt aligned_fasta.afa contigs.fasta"
    parser = optparse.OptionParser(usage=usage, version=('%prog version ' + globals()['__version__']), description=globals()['__doc__'])
    parser.add_option('-o', '--output_dir', default=None,
                      help='Directory to store output files, default is aligned_fasta directory')
    parser.add_option('-q', '--quiet', action='store_true', default=False,
                      help='Quiet, replace all deletions found, no prompts')
    parser.add_option('-v', '--verbose', action='store_true', default=False,
                      help='verbose output')
    try:
        (options, args) = parser.parse_args(argv[1:])
        if len(args) < 3:
            parser.error('Must specify a contig composition file, an aligned fasta file, and a contigs fasta file')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    contig_composition_file = args[0]
    aligned_fasta_input = args[1]
    contigs_fasta = args[2]
    (output_dir, output_basename) = os.path.split(aligned_fasta_input)
    (output_basename, dummy_ext) = os.path.splitext(output_basename)
    if options.output_dir:
        output_dir = options.output_dir

#    output_basepath = os.path.join(output_dir, output_basename)
#    aligned_fasta_output = "%s.deletions_added.afa" % output_basepath
#    fasta_output_filename = "%s.deletions_added.fasta" % fasta_output_basepath

    if options.verbose:
        print "VAMP_HOME=%s" % VAMP_HOME

    # Read contigs file
    fasta_dict = Bio.SeqIO.index(contigs_fasta, 'fasta')
    contigs = []
    sequence_source_id = None
    print "\nContig Composition Summary"
    print ContigComposition.tab_headings()
    for contig in vamp.utils.read_contig_composition_summary(contig_composition_file):
        contigs.append(contig)
        print contig.to_tab()
        sequence_source_id = contig['seq']
    fasta_output_filename = os.path.join(output_dir, "%s_deletions_added.fasta" % sequence_source_id)

    # Read Alignment
    alignment = Bio.AlignIO.read(aligned_fasta_input, 'fasta')
    original_sequence_alignment = None
    original_sequence_alignment_idx = None
    for idx, seq in enumerate(alignment):
        if seq.id == sequence_source_id:
            original_sequence_alignment = seq
            original_sequence_alignment_idx = idx
    assembled_sequence_str = re.sub('-', '', str(original_sequence_alignment.seq))
    if options.verbose:
        print "\nOriginal Sequence of %s:\n%s\n" % (original_sequence_alignment.id, assembled_sequence_str)
    assembled_sequence = SeqRecord(Seq(assembled_sequence_str, generic_dna), id=original_sequence_alignment.id)

    # Find deletions
    replacements = vamp.utils.find_deletions(contigs, True)

    # Allow user to select them
    selected_replacements = []
    for (contigs_to_be_replaced, replacement) in replacements:
        print ("\nContig %s:%s-%s is missing" % (replacement['contig'], replacement['contig_start'], replacement['contig_end']))
        if not contigs_to_be_replaced:
            print ("\tinsert\n\t\t%s\n\tat %s:%s" % (replacement.to_tab(), replacement['seq'], replacement['start']))
        else:
            print ("\treplace")
            for idx in contigs_to_be_replaced:
                contig_to_be_replaced = contigs[idx]
                print ("\t\t%s" % contig_to_be_replaced.to_tab())
            print "\twith\n\t\t%s" % replacement.to_tab()
        if (options.quiet or confirm("\t? ")):
            selected_replacements.append(replacement)

    # Generate updated sequence
    updated_sequence = vamp.utils.update_sequence_with_replacements(assembled_sequence, selected_replacements, fasta_dict)
    if options.verbose:
        print "\nUpdated Sequence of %s:\n%s\n" % (original_sequence_alignment.id, updated_sequence.seq)
    SeqIO.write(updated_sequence, fasta_output_filename, "fasta")

    # TODO Generate updated alignment

    # TODO Generate updated contig composition summary


def confirm(prompt=None, resp=False):
    """prompts for yes or no response from the user. Returns True for yes and
    False for no.

    'resp' should be set to the default value assumed by the caller when
    user simply types ENTER.

    >>> confirm(prompt='Create Directory?', resp=True)
    Create Directory? [y]|n:
    True
    >>> confirm(prompt='Create Directory?', resp=False)
    Create Directory? [n]|y:
    False
    >>> confirm(prompt='Create Directory?', resp=False)
    Create Directory? [n]|y: y
    True

    """

    if prompt is None:
        prompt = 'Confirm'

    if resp:
        prompt = '%s [%s]|%s: ' % (prompt, 'y', 'n')
    else:
        prompt = '%s [%s]|%s: ' % (prompt, 'n', 'y')

    while True:
        ans = raw_input(prompt)
        if not ans:
            return resp
        if ans not in ['y', 'Y', 'n', 'N']:
            print 'please enter y or n.'
            continue
        if ans == 'y' or ans == 'Y':
            return True
        if ans == 'n' or ans == 'N':
            return False

if __name__ == '__main__':
    sys.exit(main())
