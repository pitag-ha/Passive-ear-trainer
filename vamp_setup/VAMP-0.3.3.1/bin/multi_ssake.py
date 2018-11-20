#!/usr/bin/env python
""" Run various iterations of SSAKE, varying input files and parameters
Collect results into single list of contigs
"""
import ConfigParser
import itertools
import optparse
import os
import subprocess
import sys
import vamp.utils

__version__ = "0.3"
__author__ = "Lance Parsons"
__author_email__ = "lparsons@princeton.edu"
__copyright__ = "Copyright 2012, Lance Parsons"
__license__ = "BSD 2-Clause"


# import TQSfastq


def main(argv=None):
    if argv is None:
        argv = sys.argv

    default_config_dir = os.path.join(vamp.utils.get_vamp_home(), 'makefiles')
    default_multi_ssake_config = os.path.join(default_config_dir, 'multi_ssake.config')
    if not os.path.exists(default_multi_ssake_config):
        default_multi_ssake_config = os.path.join(default_config_dir, 'multi_ssake.config.template')

    usage = "Usage: %prog [options] forward_reads reverse_reads"
    parser = optparse.OptionParser(usage=usage, version='%prog version ' + globals()['__version__'], description=globals()['__doc__'])
    parser.add_option('-i', '--insert_size', default=500, help='Mean insert '
                      'size for paired reads (default: %default)')
    parser.add_option('--config', default=default_multi_ssake_config,
                      help='multi_ssake configuration file that specifies options [default: %default]')
    parser.add_option('--vamp_config', default=vamp.utils.default_vamp_config(),
                      help='vamp configuration file that specifies options [default: %default]')
    parser.add_option('--qsub', action='store_true', default=False, help='Use '
                      'qsub to submit commands to cluster (default: %default)')
    parser.add_option('--untrimmed', action='store_true', default=False, help='Include SSAKE assembly of untrimmed reads (not '
                      'recommended)')
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
    parser.add_option('-d', '--debug', action='store_true', default=False, help='debug (do not execute)')
    global options
    try:
        (options, args) = parser.parse_args(argv[1:])
        if len(args) != 2:
            parser.error('Please specify the at forward and reverse fastq files')
    except SystemExit:  # Prevent exit when calling as function
        return 2

    print ("")
    print ("Forward Read File: %s" % args[0])
    if len(args) > 1:
        print ("Reverse Read File: %s" % args[1])
    if len(args) > 2:
        print ("Unpaired Read File: %s" % args[2])

    global config
    config = ConfigParser.ConfigParser()
    config_dir = os.path.join(vamp.utils.get_vamp_home(), 'makefiles')
    try:
        config.readfp(open(os.path.join(config_dir, 'multi_ssake.config')))
    except IOError:
        config.readfp(open(os.path.join(config_dir, 'multi_ssake.config.template')))
    if options.config:
        print ("Using options from config file %s" % options.config)
        config.read(options.config)

    global vamp_config
    vamp_config = vamp.utils.parse_makefile(options.vamp_config)

    # Convert input fastq to fastq
    fasta_files = runFastqToFasta(args, options.qsub, config.get('DEFAULT', 'qsub_options'), options.debug)

    # TQSFastq.py
    tqs_thresholds = map(int, config.get('TQS', 'min_qual').split())
    tqs_consecs = map(int, config.get('TQS', 'min_consec').split())
    print("TQS Quality Threshold(s): %s" % ','.join(map(str, tqs_thresholds)))
    print("TQS Consecutive Threshold(s): %s" % ','.join(map(str, tqs_consecs)))
    tqs_output_files = runTqsFastq(args, tqs_thresholds, tqs_consecs, options.qsub, config.get('TQS', 'qsub_options'), options.debug)

    # SSAKE
    ssake_input_filesets = tqs_output_files
    if options.untrimmed:
        ssake_input_filesets.insert(0, fasta_files)
    for fileset in ssake_input_filesets:
        print 'Fileset: %s' % fileset
    ssake_min_overlaps = map(int, config.get('SSAKE', 'min_overlap').split())
    print("SSAKE Minimum Overlap Value(s): %s" %
          ','.join(map(str, ssake_min_overlaps)))
    ssake_trims = map(int, config.get('SSAKE', 'trim').split())
    print("SSAKE Trim Value(s): %s" % ','.join(map(str, ssake_trims)))
    ssake_jids = runSsake(ssake_input_filesets, options.insert_size, ssake_min_overlaps, ssake_trims, options.qsub, options.debug)

    # Combine SSAKE output and collapse (fastx_collapser?)
    if options.qsub:
        combineSsake(ssake_jids)
    else:
        combineSsake()
    return 0


def runSsake(input_filesets, insert_size=500, min_overlap=[20], trim=[0], qsub=False, debug=False):
    '''
    Executes SSAKE on the input filesets
    Args:
        input_filesets - list of filesets where each fileset is a list of
        tuples, one for each read (forward, reverse, unpaired)
        min_overlap - optional list of minimum overlap values (def 20)
        trim - optional list of trim values (def 0)
        qsub - boolean, use qsub for execution (default: False)
    Return:
        ssake_jids - list of ssake job_ids (if qsub) or return codes
    '''
    count = 1
    ssake_jids = []
    for fileset in input_filesets:
        tempdir = 'ssake_tempdir_%s' % count
        print "Tempdir: %s" % tempdir
        print "\tFileset: %s" % fileset
        count += 1

        # Make Paired
        cmd = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'makePairedOutput2UNEQUALfiles_vamp.pl')
        files = []
        hold_jids = []
        for fileinfo in fileset:
            files.append(os.path.join('..', fileinfo[0]))
            if (fileinfo[2] is not None) or int(fileinfo[2]) != 0:
                hold_jids.append(str(fileinfo[2]))
        if len(hold_jids) > 0:
            make_paired_qsub_options = '-cwd -terse -hold_jid %s' % ",".join(hold_jids)
        cmd_options = "%s %s" % (" ".join(files), insert_size)
        # options = "-m %s -t %s" % (m, t)
        if (not os.path.exists(tempdir)):
            os.mkdir(tempdir)
        os.chdir(tempdir)
        make_paired_retcode = execute(cmd, cmd_options, qsub, make_paired_qsub_options, debug)

        # For each input set, iterate over SSAKE options and run
        ssake_path = config.get('SSAKE', 'ssake_path')
        for (m, t) in itertools.product(min_overlap, trim):
            ssake_options = config.get('SSAKE', 'ssake_options')
            ssake_qsub_options = config.get('SSAKE', 'qsub_options')
            ssake_options = '-m %s -t %s -p 1 -f paired.fa -g unpaired.fa %s' % (m, t, ssake_options)
            if qsub:
                ssake_qsub_options = '%s -hold_jid %s' % (ssake_qsub_options, make_paired_retcode)
            ssake_retcode = execute(ssake_path, ssake_options, qsub, ssake_qsub_options, debug)
            ssake_jids.append(ssake_retcode)

        # Change directory back "up"
        os.chdir('..')
    return ssake_jids


def combineSsake(hold_jids=None):
    # Combine all contigs and mergedcontigs
    combine_script_filename = 'combine_contigs.sh'
    combine_script = open(combine_script_filename, 'wb')
    combine_script.write('#!/bin/bash\n')
    combine_script.write('find . -name "*.contigs" -exec cat {} \; | awk \'NR%2==0\' | sort -S1G | uniq | awk \'{ print length($0) "\t" $0 }\' | sort -k1 -r -n | awk \'{printf(">%s-%s\\n%s\\n",NR,$1,$2)}\' > ssake_combined.contigs\n')
    combine_script.write('find . -name "*.mergedcontigs" -exec cat {} \; | awk \'NR%2==0\' | sort -S1G | uniq | awk \'{ print length($0) "\t" $0 }\' | sort -k1 -r -n | awk \'{printf(">%s-%s\\n%s\\n",NR,$1,$2)}\' > ssake_combined.mergedcontigs\n')

    combine_script.close()
    combine_qsub_options = config.get('DEFAULT', 'qsub_options')
    if hold_jids:
        combine_qsub_options = '%s -hold_jid %s' % (combine_qsub_options, ",".join(map(str, hold_jids)))
    execute('bash', combine_script_filename, options.qsub, combine_qsub_options, options.debug)


def runFastqToFasta(input_files, qsub=False, qsub_options='-terse -cwd', debug=False):
    '''
    Executes format_convert.py script to convert each input input_file
    from fastq to fasta format
    Returns list of tuples, one per input_file
    Tuples are (output_file, log_file, qsub_jid [or return_code])
    '''
    output_files = []
    for input_file in input_files:
        (dirname, base, extension) = splitGzFilename(input_file)  # @UnusedVariable
        output_file = '%s.fasta' % base
        cmd = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'fastq_to_fasta.py')
        fastq_to_fasta_options = '"%s" "%s"' % (input_file, output_file)
        retcode = execute(cmd, fastq_to_fasta_options, qsub, qsub_options, debug)
        output_files.append((output_file, None, retcode))
    return output_files


def runTqsFastq(input_files, thresh=[10], consec=[20], qsub=False, qsub_options='-terse -cwd', debug=False):
    '''
    Executes the TqsFastq script on the specified input file
    Args:
        input_files - list of filenames of fastq files to trim
        thresh - optional list of quality thresholds (def 10)
        consec - options list of consecutive thresholds (def 20)
    Return:
        output_files - list of output files produced
            each member is a list of tuples, one per input file
            each tuple consists of output file, log file, and qsub jobid
            if thresh and consec are None, list of inputs is returned
    '''
    output_files = []
    output_fileset = []

    for (t, c) in itertools.product(thresh, consec):
        output_fileset = []
        for f in input_files:
            # (dirname, base, extension) = splitGzFilename(f)
            cmd = vamp_config.get('TQSFASTQ',
                                  os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                               'TQSfastq_vamp.py'))
            options = "-t %s -c %s -f %s" % (t, c, f)
            retcode = execute(cmd, options, qsub, '-terse -cwd', debug)
            # TQSfastq.main(['-f', f, '-t', str(t), '-c', str(c), '-q'])
            output_root = '%s_T%sC%s.trim' % (f, t, c)
            output_fileset.append(('%s.fasta' % output_root,
                                   '%s.log' % output_root, retcode))
        output_files.append(output_fileset)
    return output_files


def splitGzFilename(filename):
    '''
    Split filename into root, base, and extention
    Ext is standard extention, but includes .gz
    e.g. dir/file.fastq.gz -> ('dir', 'file', 'fastq.gz')
    '''
    (dirname, basename) = os.path.split(filename)  # @UnusedVariable
    (root, ext) = os.path.splitext(basename)
    if (ext == '.gz'):
        (root, ext) = os.path.splitext(root)
        ext = ext + '.gz'
    return (dir, root, ext)


# class SsakeRun:
#    """Class to hold ssake run options and expected output files"""
#
#    def __init__(self):
#        self.config = []
#        self.output_files = []


def execute(cmd, cmd_options, qsub, qsub_options, debug):
    retcode = 0
    if qsub:
        cmd_options = '%s %s %s' % (qsub_options, cmd, cmd_options)
        cmdline = 'qsub %s' % cmd_options
        print cmdline
        if not options.debug:
            retcode = subprocess.Popen(cmdline, shell=True, stdout=subprocess.PIPE).communicate()[0].strip()
    else:
        cmdline = '%s %s' % (cmd, cmd_options)
        print cmdline
        if not options.debug:
            retcode = subprocess.call(cmdline, shell=True)
    if (retcode != 0) and (not qsub):
        print >> sys.stderr, "Child was terminated by signal", -retcode
        print >> sys.stderr, "Error during execution of: %s"
        raise Exception("Error executing command: %s" % cmdline)
    print 'Return code: %s' % retcode
    return retcode


if __name__ == '__main__':
    sys.exit(main())
