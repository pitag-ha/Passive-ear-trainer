#!/usr/bin/env python
""" Setup VAMP analysis. Prompts for input files and options, creates makefile to run analysis
"""
from string import Template
import optparse
import os
import sys
import vamp.utils

__version__ = "0.1"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2011, Lance Parsons"
__license__ = "BSD Simplified"


def main(argv=None):
    if argv is None:
        argv = sys.argv

    usage = "Usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage, version='%prog version ' + globals()['__version__'], description=globals()['__doc__'])
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
    try:
        (options, args) = parser.parse_args(argv[1:])  # @UnusedVariable
        #if len(args) < 1:
        #    parser.error('argument missing')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    # Get VAMP HOME
    VAMP_HOME = vamp.utils.get_vamp_home()
    print "VAMP home directory: %s" % VAMP_HOME

    # Get Forward Read File
    while True:
        forward_read = raw_input("\nFilename of forward read: ")
        forward_read = os.path.expanduser(forward_read.strip())
        if os.path.exists(forward_read):
            forward_read = os.path.abspath(forward_read)
            break
        else:
            print "The filename you entered does not exist: '%s'" % forward_read
    print "Forward Read Fastq File: %s" % forward_read

    # Get Reverse Read File
    while True:
        reverse_read = raw_input("\nFilename of reverse read [blank for none]: ")
        reverse_read = os.path.expanduser(reverse_read.strip())
        if os.path.exists(reverse_read):
            reverse_read = os.path.abspath(reverse_read)
            if os.path.samefile(forward_read, reverse_read):
                print "Forward and Reverse reads must be different files"
                continue
            else:
                break
        elif reverse_read == '':
            break
        else:
            print "The filename you entered does not exist: '%s'" % reverse_read

    if reverse_read == '':
        print "No reverse read defined, single read analysis"
    else:
        print "Reverse Read Fastq File: %s" % reverse_read

    # Output location
    while True:
        output_dir = raw_input("\nEnter directory to save output files to: ")
        output_dir = os.path.expanduser(output_dir.strip())
        if os.path.isdir(output_dir):
            break
        else:
            print "The directory you entered does not seem to exist: '%s'" % output_dir
    output_dir = os.path.abspath(output_dir)
    print "Output directory: %s" % output_dir

    output_prefix = raw_input("\nEnter the prefix for output files: ")
    output_prefix = output_prefix.strip()

    # Get host reference

    while True:
        host_reference = raw_input("\nHost reference FASTA file or BOWTIE index basename: ")
        host_reference = os.path.expanduser(host_reference.strip())
        host_reference_bowtie_ebwt = "%s.1.ebwt" % host_reference
        bowtie_found = False
        set_host_reference_bowtie_index = ""
        set_host_reference_fasta = ""
        if os.path.exists(host_reference_bowtie_ebwt):
            host_reference_bowtie = os.path.abspath(host_reference)
            bowtie_found = True
            break
        if os.path.exists(host_reference):
            host_reference_fasta = os.path.abspath(host_reference)
            break
        else:
            print "Cannot find bowtie index or fasta file: '%s'" % host_reference
    if (bowtie_found):
        print "Host reference BOWTIE index: %s" % host_reference_bowtie
        set_host_reference_bowtie_index = "HOST_REFERENCE_BOWTIE_INDEX:=%s" % host_reference_bowtie
    else:
        print "Host reference FASTA file: %s" % host_reference_fasta
        set_host_reference_fasta = "HOST_REFERENCE_FASTA:=%s" % host_reference_fasta

    # Print configuration file
    template_file = open(os.path.join(VAMP_HOME, 'makefiles', 'analysis_defaults.mk'))
    template = Template(template_file.read())
    template_string = template.safe_substitute(
        vamp_home_dir=VAMP_HOME,
        forward_read=forward_read,
        reverse_read=reverse_read,
        output_dir=output_dir,
        output_base=output_prefix,
        set_host_reference_fasta=set_host_reference_fasta,
        set_host_reference_bowtie_index=set_host_reference_bowtie_index)
    makeFilename = os.path.join(output_dir, '%s_preprocess.mk' % output_prefix)
    print "Writing preprocessing makefile to: %s" % makeFilename
    makefile = open(makeFilename, 'wb')
    makefile.write(template_string)
    makefile.close()
    return 1


if __name__ == '__main__':
    sys.exit(main())
