WaveBlocksND
============

Reusable building blocks for simulations with semiclassical wavepackets
for solving the time-dependent Schrödinger equation.

Installation
============

Using Python's `virtualenv` for installation into a directory called `waveblocks`:

```Shell
virtualenv -p /usr/bin/python3.4 --system-site-packages waveblocks
cd waveblocks
source ./bin/activate
pip install git+https://github.com/WaveBlocks/WaveBlocksND.git#egg=WaveBlocksND
```

Note: Python in version >= 3.4 is fully supported.

License
=======

This project is covered by the 3-Clause BSD License.

Citation
========

This repository contains source code for published scientific work.
The usual citation requirements apply. For citation of this project
please use the following bibtex snippet:

```Tex
@misc{waveblocksnd,
  author       = {R. Bourquin and V. Gradinaru},
  howpublished = {\url{https://github.com/WaveBlocks/WaveBlocksND}},
  title        = {{WaveBlocks}: Reusable building blocks for simulations with semiclassical wavepackets},
  url          = {https://github.com/WaveBlocks/WaveBlocksND},
  year         = {2010 - 2016},
}
```
