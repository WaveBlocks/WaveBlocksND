from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name = 'WaveBlocksND',
      version = '0.5',
      description = u'Reusable building blocks for simulations with semiclassical wavepackets.',
      long_description = readme(),
      author = 'R. Bourquin',
      author_email = 'raoul.bourquin@sam.math.ethz.ch',
      url = 'https://github.com/WaveBlocks/WaveBlocksND',
      packages = ['WaveBlocksND'],
      scripts = [],
      classifiers = [
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Environment :: Console',
        'Programming Language :: Python :: 2.7',
      ],
      license = 'BSD',
      install_requires = [
        'h5py',
        'sympy',
        'numpy',
        'scipy',
      ])
