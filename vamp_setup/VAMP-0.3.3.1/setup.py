from distribute_setup import use_setuptools
use_setuptools()

from setuptools import setup
# from distutils.core import setup

setup(
    name='VAMP',
    version='0.3.3.1',
    author='Lance Parsons',
    author_email='lparsons@princeton.edu',
    packages=['vamp', 'seq_utils', 'galaxy'],
    scripts=['bin/compare_genomes.py',
             'bin/fastq_to_fasta.py',
             'bin/translate_cds.py',
             'bin/find_contig_deletions.py',
             'bin/gff2gtf_simple.py',
             'bin/maf_net.py',
             'bin/multi_ssake.py',
             'bin/preprocess.py',
             'bin/setup_preprocessing.py',
             'bin/TQSfastq_vamp.py',
             'bin/clean_multi_ssake.sh',
             'bin/makePairedOutput2EQUALfiles_vamp.pl',
             'bin/makePairedOutput2UNEQUALfiles_vamp.pl'],
    url='https://bitbucket.org/lance_parsons/vamp/',
    license='LICENSE.txt',
    description='Virus AsseMbly Pipeline.',
    long_description=open('README.txt').read(),
    data_files=[('makefiles', ['makefiles/analysis_defaults.mk',
                               'makefiles/config.mk.template',
                               'makefiles/preprocess_single_read.mk',
                               'makefiles/preprocess.mk',
                               'makefiles/test.mk'])],
    install_requires=[
        "paired-sequence-utils >= 0.1",
        "cutadapt >= 1.2.1",
        "pybedtools == 0.6",  # pybedtools 0.6.2 has bug in 'print interval'
        "bx-python >= 0.7.1",
        "BioPython >= 1.57",
        "Cython >= 0.17.4",
    ],
)
