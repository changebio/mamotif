from distutils.core import setup
from distutils.core import Extension
from distutils.command.install import INSTALL_SCHEMES

for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']
import os

os.environ['CC'] = 'g++'
os.environ['LDSHARED'] = 'g++'
VERSION = "1.0.4"
CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Enviroment :: Console',
    'Intended Audience :: Science/Research ',
    'Intended Audience :: Education ',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
	'Programming Language :: Python::2.7',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]

install_requires = [
    'numpy',
    'scipy',
    'pandas',
]


setup(
    name="MotifScan",
    description="A motif discovery package developed by Shao lab in SIBS, CAS",
    version=VERSION,
    author="Jiawei Wang",
    author_email="jerryeah@gmail.com",
    url="http://bioinfo.sibs.ac.cn/shaolab/motifscan",
    package_dir={'MotifScan': 'MotifScan'},
    packages=['MotifScan'],
    data_files=[('MotifScan/gene/hg19',['gene/hg19/refSeq.txt']),
                ('MotifScan/gene/mm9',['gene/mm9/refSeq.txt']),
                ('MotifScan/gene/tair10',['gene/tair10/refSeq.txt'])],
    scripts=['bin/MotifScan','bin/MotifCompile'],
    ext_package='MotifScan',
    ext_modules=[Extension(name="score_c",sources=['MotifScan/score.C'],extra_link_args=["--shared"])],
    classifiers=CLASSIFIERS,
)


