#!/usr/bin/env python

'''
Extracts coding sequences (cds) regions from fasta reference and gff file and
translates them into amino acid sequence, output in FASTA format to STDOUT.

Errors during translation are output to STDERR. Genes with translation errors
are not printed.
'''
import sys
import os

import optparse
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from pybedtools import BedTool
from Bio.Data import CodonTable

__version__ = "0.2"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2012, Lance Parsons"
__license__ = "BSD 2-Clause License " \
    "http://www.opensource.org/licenses/BSD-2-Clause"


def main(argv=None):
    if argv is None:
        argv = sys.argv

    usage = "Usage: %prog [options] <genes gff3> <ref fasta>"
    parser = optparse.OptionParser(usage=usage, version='%prog version ' +
                                   globals()['__version__'],
                                   description=globals()['__doc__'])
    parser.add_option('--notrans', action='store_true', default=False,
                      help='Do not translate to amino acid sequence, '
                      'output DNA')
    parser.add_option('-i', '--idattr', default='gene_id',
                      help='GFF attribute to use as gene ID. Features '
                      'with the same ID will be considered parts of the same '
                      'gene. The default "%default" is suitable for GTF '
                      'files.')
    parser.add_option('-t', '--featuretype', default=None,
                      action='append',
                      help='GFF feature type(s) (3rd column) to be used. '
                      'Specify the option multiple times for multiple '
                      'feature types. The default is "CDS" for GFF files '
                      'and "CDS" and "stop_codon" for GTF files.')
    parser.add_option('--table', default=1, help='NCBI Translation table to '
                      'use when translating DNA (see http://www.ncbi.nlm.nih.'
                      'gov/Taxonomy/Utils/wprintgc.cgi). Default: %default.')
    parser.add_option('-v', '--verbose',
                      action='store_true',
                      default=False, help='verbose output')
    try:
        (options, args) = parser.parse_args(argv[1:])
        if len(args) < 2:
            parser.error('Please specify gff and fasta files')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    # Get Arguments
    gff_filename = args[0]
    fasta_filename = args[1]

    # Read GFF file
    basename, extension = os.path.splitext(gff_filename)
    is_gtf = False
    if extension.lower() == '.gtf':
        is_gtf = True
    if options.featuretype is None:
        featuretypes = ['CDS']
        if is_gtf:
            featuretypes.append('stop_codon')
    else:
        featuretypes = options.featuretype
    gff_file = BedTool(gff_filename)
    if gff_file.file_type != 'gff':
        parser.error('Gene file format not recognized, must be GFF or GTF')
    cds_regions = gff_file.filter(lambda b: b.fields[2] in
                                  featuretypes)
    genes = defaultdict(list)
    for feature in cds_regions:
        gene_id = feature.attrs.get(options.idattr)
        if gene_id is not None:
            genes[gene_id].append(feature)
    # TODO Store gene_name to use in fasta description
    # TODO Store start position (and chromosome order) and use to sort

    # Read Fasta File
    with open(fasta_filename) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

    for gene, cds in genes.iteritems():
        #print gene
        seq = Seq("", IUPAC.unambiguous_dna)
        for c in cds:
            ref_seq = ref_recs.get(c.chrom)
            cds_seq = ref_seq[c.start:c.end]
            if c.strand == '-':
                cds_seq = cds_seq.reverse_complement()
            #print "%s:%i-%i - " % (c.chrom, c.start, c.end)
            #print cds_seq.seq
            seq = seq + cds_seq.seq

        # Attempt to translate as complete CDS, report errors
        protein_seq = None
        try:
            protein_seq = seq.translate(table=options.table, cds=True)
        except CodonTable.TranslationError as e:
            sys.stderr.write("Error found in gene '%s': %s\n" %
                             (gene, str(e)))

        # Attempt to translate without requiring complete CDS
        if protein_seq is None:
            try:
                protein_seq = seq.translate(table=options.table, cds=False)
            except Exception as e:
                sys.stderr.write("Unable to translate gene '%s': %s\n" %
                                 (gene, str(e)))
        # Print output
        output_record = None
        if options.notrans:
            output_record = SeqRecord(seq, id=gene, description='')
        else:
            output_record = SeqRecord(protein_seq, id=gene, description='')
        if output_record is not None:
            print output_record.format("fasta"),

    return 0


if __name__ == '__main__':
    sys.exit(main())
