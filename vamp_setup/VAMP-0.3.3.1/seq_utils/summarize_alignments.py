#!/usr/bin/env python
"""
SYNOPSIS

    TODO summarize_alignments multi_align_fasta reference_sequence [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    TODO This describes how to use this script. This docstring
    will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    TODO: lparsons <lparsons@princeton.edu>

LICENSE

    BSD 2-Clause License http://www.opensource.org/licenses/BSD-2-Clause

VERSION

    $Id$
"""

from Bio import AlignIO
import optparse
import os
import sys
import time
import traceback
#import copy


def main():

    global options, args
    # Read multi-fasta file
    multi_align_file = args[0]
    reference_sequence_id = args[1]

    alignment = AlignIO.read(open(multi_align_file), "fasta")
    print ("There are %i sequences in this alignment" % len(alignment))
    #for i in range(len(alignment)):
    ref_found = False
    for s in alignment:
        print "\t%i: %s" % (s.id),
        if (s.id == reference_sequence_id):
            print " (Reference)",
            ref_found = True
        print "\n",
    print ("The alignment is %i bases long" % len(alignment[0, :]))
    if not ref_found:
        sys.stderr.write("Reference sequence id '%s' not found in fasta file '%s'" % (reference_sequence_id, multi_align_file))
        return 1

    mismatch_summaries = summary_of_alignment(alignment, reference_sequence_id)
    for seq_id, mismatch_summary in mismatch_summaries.iteritems():
        print "\nSequence: %s" % seq_id
        print "Matches:       %i" % mismatch_summary['match_count']
        print "Mismatches:    %i" % mismatch_summary['mismatch_count']
        print "Change Events: %i" % mismatch_summary['noncontig_count']

        #for m in mismatch_summary['mismatches']:
        #    print "%s(%i)%s" % (m.get("reference_base"), (m.get("reference_position")+1), m.get("new_base"))


def summary_of_alignment(alignment, reference_sequence_id):
    '''
    Summarizes changes in given alignment (pairwise only)
    Input: alignment = Bio.AlignIO object
           reference_index = index of the reference sequence in alignment (default is 1)
    Output: dictionary with key for each non-reference sequence in alignment
        Each key has a dictionary with keys (match_count, mismatch_count, mismatches, contiguous_change_count)
                mismatches is list of mismatches by base: 'RefBase(RefPos)NewBase'
                contiguout_change_count is the number of contiguous change "events"
    '''
    mismatch_summaries = {}
    reference_index = None
    for i in range(len(alignment)):
        if (alignment[i].id != reference_sequence_id):
            mismatch_summaries[alignment[i].id] = dict(match_count=0, mismatch_count=0, mismatches=[], noncontig_count=0)
        else:
            reference_index = i

    for sequence in alignment:
        if sequence.id == reference_sequence_id:
            continue
        event_list = []
        event = dict(pos=0, ref='', alt='', align_pos=0, type=None, snp=False, insertion=False, deletion=False)
        mismatch_event_count = 0
        refposition = 0
        previous_position = -2
        for position in range(len(alignment[reference_index])):
            refbase = alignment[reference_index][position]
            altbase = sequence[position]
            if (refbase != '-'):
                refposition += 1
            if (refbase != altbase):
                mismatch_summaries[sequence.id]['mismatch_count'] += 1
                mismatch_summaries[sequence.id]['mismatches'].append(dict(pos=refposition, ref=refbase, alt=altbase))

                if refposition > previous_position + 1:
                    mismatch_event_count += 1
                    if event['pos'] != 0:
                        # Append previous event to list
                        variation = parse_event(event, alignment[reference_index], sequence)
                        event_list.append(variation)
                    event = dict(pos=0, ref='', alt='', align_pos=0, type=None, snp=False, insertion=False, deletion=False)
                    event['pos'] = refposition
                    event['ref'] = refbase
                    event['alt'] = altbase
                    event['align_pos'] = position
                else:
                    event['ref'] += refbase
                    event['alt'] += altbase
                # Set flags used to determine type
                if refbase == '-':
                    event['insertion'] = True
                elif altbase == '-':
                    event['deletion'] = True
                else:
                    event['snp'] = True
                previous_position = refposition
        else:
            mismatch_summaries[sequence.id]['match_count'] += 1
        # Summarize mismatch info
        mismatch_summaries[sequence.id]['noncontig_count'] = len(event_list)
        mismatch_summaries[sequence.id]['noncontig_list'] = event_list
        mismatch_summaries[sequence.id]['noncontigs_string'] = mismatch_string(event_list)
        mismatch_summaries[sequence.id]['mismatches_string'] = mismatch_string(mismatch_summaries[sequence.id]['mismatches'])
    return(mismatch_summaries)


#def summary_of_alignment_orig(alignment, reference_sequence_id):
#    '''
#    Summarizes changes in given alignment (pairwise only)
#    Input: alignment = Bio.AlignIO object
#           reference_index = index of the reference sequence in alignment (default is 1)
#    Output: dictionary with key for each non-reference sequence in alignment
#        Each key has a dictionary with keys (match_count, mismatch_count, mismatches, contiguous_change_count)
#                mismatches is list of mismatches by base: 'RefBase(RefPos)NewBase'
#                contiguout_change_count is the number of contiguous change "events"
#    '''
#    mismatch_summaries = {}
#    reference_index = None
#    for i in range(len(alignment)):
#        if (alignment[i].id != reference_sequence_id):
#            mismatch_summaries[alignment[i].id] = dict(match_count=0, mismatch_count=0, mismatches=[], noncontig_count=0)
#        else:
#            reference_index = i
#    refposition = 0
#    for position in range(len(alignment[0, :])):
#        column = alignment[:, position]
#        refbase = column[reference_index]
#        if refbase != '-':
#            refposition += 1
#        for seq in range(len(column)):
#            if (seq == reference_index):
#                continue
#            seq_id = alignment[seq].id
#            if (refbase != column[seq]):
#                mismatch_summaries[seq_id]['mismatch_count'] += 1
#                mismatch_summaries[seq_id]['mismatches'].append(dict(reference_position=refposition, reference_base=refbase, new_base=column[seq]))
#            else:
#                mismatch_summaries[seq_id]['match_count'] += 1
#        #print alignment[:,position]
#
#    for seq_id in mismatch_summaries.keys():
#        for i in range(len(alignment[:, 0])):
#            if (alignment[i].id == seq_id):
#                seq_id_index = i
#                break
#        noncontig_mismatches_info = noncontig_mismatches(mismatch_summaries[seq_id]['mismatches'], alignment[reference_index], alignment[seq_id_index])
#        mismatch_summaries[seq_id]['noncontig_count'] = noncontig_mismatches_info[0]
#        mismatch_summaries[seq_id]['noncontig_list'] = noncontig_mismatches_info[1]
#        mismatch_summaries[seq_id]['noncontigs_string'] = mismatch_string(mismatch_summaries[seq_id]['noncontig_list'])
#        mismatch_summaries[seq_id]['mismatches_string'] = mismatch_string(mismatch_summaries[seq_id]['mismatches'])
#    return(mismatch_summaries)


#def noncontig_mismatches(mismatches, reference_sequence, alternate_sequence):
#    event_list = list()
#    event = dict(reference_position=0, reference_base='', new_base='')
#    count = 0
#    previous_position = -2
#    for i in range(len(mismatches)):
#        mismatch = mismatches[i]
#        if mismatch["reference_position"] > (previous_position + 1):
#            count += 1
#            if event['reference_position'] != 0:
#                event = parse_event(event, reference_sequence, alternate_sequence)
#                event_list.append(copy.copy(event))
#            #event_list.append(event)
#            event['reference_position'] = mismatch["reference_position"]
#            event['reference_base'] = mismatch["reference_base"]
#            event['new_base'] = mismatch["new_base"]
#        else:
#            event['reference_base'] += mismatch["reference_base"]
#            event['new_base'] += mismatch["new_base"]
#        previous_position = mismatch["reference_position"]
#    event_list.append(copy.copy(event))
#    return (count, event_list)


def parse_event(event, reference_sequence, alternate_sequence):
    '''
    Parse a simple event with reference_position, reference_base, and new_base
    and determine the type and add padding if necessary (for VCF compatibility)
    '''
    refgaps = event.get('ref').count('-')
    altgaps = event.get('alt').count('-')
    if refgaps > 0 or altgaps > 0:
        if event['pos'] > 1:  # Add padding base to start of variant
            event['pos'] -= 1
            event['ref'] = '%s%s' % (reference_sequence[event['align_pos'] - 1], event['ref'])
            event['alt'] = '%s%s' % (alternate_sequence[event['align_pos'] - 1], event['alt'])
        #else:  # Add padding to end of variant
        #    event['ref'] = '%s%s' % (event['ref'], reference_sequence[aligned_position + len(event['ref'])])
        #    event['alt'] = '%s%s' % (event['alt'], alternate_sequence[aligned_position + len(event['alt'])])

    if event.get('snp'):
        if (len(event.get('ref')) == 1 and len(event.get('alt')) == 1):
            event['type'] = 'SNP'
        else:
            if refgaps == 0 and altgaps == 0:
                event['type'] = 'MNP'
            else:
                event['type'] = 'COMPLEX'
    else:
        if refgaps > 0 or altgaps > 0:
            event['type'] = 'INDEL'
    return event


def mismatch_string(mismatches):
    mismatch_string_list = list()
    for m in mismatches:
        mismatch_string_list.append("%s(%i)%s" % (m.get("ref"), (m.get("pos") + 1), m.get("alt")))
    return ', '.join(mismatch_string_list)


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'], version='$Id$')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        #parser.add_option('-r', '--reference', type='int', metavar='INDEX', dest='reference_index', default=0, help='Index of sequence to use a reference (0-based)')
        (options, args) = parser.parse_args()
        if len(args) < 2:
            parser.error('Please specify an aligned fasta file and a reference sequence id')
        if options.verbose:
            print time.asctime()
        main()
        if options.verbose:
            print time.asctime()
            print 'TOTAL TIME IN MINUTES:',
            print (time.time() - start_time) / 60.0
        sys.exit(0)
    except KeyboardInterrupt, e:  # Ctrl-C
        raise e
    except SystemExit, e:  # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)
