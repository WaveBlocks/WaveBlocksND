from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='WaveBlocksND',
      version='0.5',
      description=u'Reusable building blocks for simulations with semiclassical wavepackets for solving the time-dependent Schrödinger equation.',
      long_description=readme(),
      classifiers=[
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
      ],
      url='https://github.com/WaveBlocks/WaveBlocksND',
      author='R. Bourquin',
      license='BSD',
      packages=['WaveBlocksND'],
      install_requires=[
        'h5py',
        'sympy',
        'numpy',
        'scipy',
      ],
      include_package_data=True,
      zip_safe=False)
