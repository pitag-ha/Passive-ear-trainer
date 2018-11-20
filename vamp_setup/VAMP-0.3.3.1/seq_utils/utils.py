"""
Utility classes and methods for working with sequence data
"""


def convert_interval_nongapped_to_gapped(seq, start, end, include_end_gaps=False):
    '''
    Take position without gaps and return position with gaps
    Uses 0-based positions
    '''
    #import pdb
    #pdb.set_trace()
    num_nongap = 0
    position = 0
    for c in seq:
        if c != '-':
            num_nongap += 1
            if num_nongap == start + 1:
                start_including_gaps = position
            if num_nongap == end:
                end_including_gaps = position + 1
        position += 1
    #pdb.set_trace()
    if include_end_gaps:
        for c in seq[:start_including_gaps][::-1]:
            if c == '-':
                start_including_gaps -= 1
                #end_including_gaps += 1
            else:
                break
        for c in seq[end_including_gaps:len(seq)]:
            if c == '-':
                end_including_gaps += 1
            else:
                break
    return start_including_gaps, end_including_gaps


def convert_interval_gapped_to_nongapped(seq, start, end):
    '''
    Take position with gaps and return position without gaps
    Uses 0-based positions
    '''
    gaps_before_start = seq[0:start].count('-')
    gaps_before_end = seq[0:end].count('-')
    start_without_gaps = start - gaps_before_start
    end_without_gaps = end - gaps_before_end
    return start_without_gaps, end_without_gaps


def get_seq_by_id_from_list(alignment, sequence_id):
    for seq in alignment:
        if seq.id == sequence_id:
            return seq


class GenomicRegion():
    """
    A genomic region specified by chr:start-end, using 1-based cooredinates
    """

    def __init__(self, region_string):
        (sequence_id, position_string) = region_string.split(':', 1)
        (start, end) = position_string.split('-', 1)
        self.sequence_id = sequence_id
        self.start = int(start)
        self.end = int(end)

    def chrom(self):
        return self.sequence_id
