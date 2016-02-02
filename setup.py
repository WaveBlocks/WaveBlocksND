from setuptools import setup, find_packages
import os


def readme():
    with open('README.md') as f:
        return f.read()


def find_scripts(*paths):
    scripts = []
    for path in paths:
        for fi in os.listdir(path):
            if fi.endswith(".py"):
                fi = os.path.join(path, fi)
                with open(fi) as f:
                    l = f.readlines()[0]
                    if l.startswith("#!"):
                          scripts.append(fi)

    return scripts


setup(name = 'WaveBlocksND',
      version = '0.5',
      description = u'Reusable building blocks for simulations with semiclassical wavepackets.',
      long_description = readme(),
      author = 'R. Bourquin',
      author_email = 'raoul.bourquin@sam.math.ethz.ch',
      url = 'https://github.com/WaveBlocks/WaveBlocksND',
      packages = find_packages(),
      classifiers = ['Development Status :: 4 - Beta',
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
        'h5py >= 2.4.0',
        'matplotlib >= 1.4.3',
        #'mayavi',
        'numpy >= 1.9.2',
        'scipy >= 0.14.1',
        'sympy >= 0.7.6',
      ],
      include_package_data = True,
      scripts = find_scripts('./scripts/scripts', './scripts/plotters/') + ['./scripts/sh/make_video.sh']
)
