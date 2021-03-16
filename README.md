# MaRs: Motif-based aligned RNA searcher

[![Build Status](https://github.com/seqan/mars/workflows/Mars%20CI/badge.svg)](https://github.com/seqan/mars/actions?query=branch%3Amaster+workflow%3A%22Mars%20CI%22)
[![License](https://img.shields.io/github/license/seqan/mars)](https://github.com/seqan/mars/blob/master/LICENSE.md)
[![Platforms](https://img.shields.io/badge/platform-Linux%20%7C%20BSD%20%7C%20MacOS-informational.svg)](https://github.com/seqan/mars/blob/master/README.md)

Mars is a tool that reads a structural multiple RNA alignment (e.g. from LaRA) and derives fuzzy stem loop descriptors
from it. These descriptors are then subject to a search in a genomic database and Mars returns the hits where the
RNA structure is found, accompanied with a quality value for each hit.


## Download instructions

Clone the repository and use the *-\-recurse-submodules* option for downloading SeqAn as submodule.

```commandline
git clone --recurse-submodules https://github.com/seqan/mars.git
```

Alternatively, you can download a zip package of the repository via the green button at the top of the github page.
If you do so, please unzip the file into a new subdirectory named *mars* and download the dependencies separately.

## Requirements

* platforms: Linux, MacOS
* compiler: gcc ≥ 7
* cmake ≥ 3.8

Mars is dependent on the following libraries:

* [SeqAn 3.0.2](https://github.com/seqan/seqan3.git)
* [IPknot](http://rtips.dna.bio.keio.ac.jp/ipknot) (shipped in the _lib_ directory) with the following dependencies:
  * [ViennaRNA 2](https://www.tbi.univie.ac.at/RNA)
  * [GLPK (GNU Linear Programming Kit) 4](https://www.gnu.org/software/glpk)
  * [Boost C++ Libraries](https://www.boost.org)

*Note:* Users reported problems with installing ViennaRNA, so we provide some hints here.

1. Install the [GNU MPFR Library](https://www.mpfr.org) first.
2. Exclude unnecessary components of ViennaRNA:
   `./configure --without-swig --without-kinfold --without-forester --without-rnalocmin --without-gsl`
3. If you have linker issues use
   `./configure --disable-lto`
4. If your system supports SSE4.1 instructions then we recommend
   `./configure --enable-sse`

If you have further suggestions, we are happy to add them here.

## Build instructions

Please create a new directory and build the program for your platform.

1. create a build directory and visit it: `mkdir build && cd build`
2. run cmake: `cmake ../mars`
3. build the application: `make`
4. optional: build and run the tests: `make test`
5. optional: build the api documentation: `make doc`

## Usage

After building the application binary, running Mars is as simple as

```commandline
bin/mars msa.aln -g genome.fasta
```

The resulting genome positions are printed to stdout.
If you want to store the result in a file instead, please use the *-o* option or redirect the output.

```commandline
bin/mars msa.aln -g genome.fasta -o result.txt
bin/mars msa.aln -g genome.fasta  > result.txt
```

We recommend you to specify the number of threads with the *-j* option, in order to enable parallel execution.
If you specify *-j 0* the program tries to detect the maximal number of threads available on your machine.

```commandline
bin/mars msa.aln -g genome.fasta -j 0
```

For a list of options, please see the help message:

```commandline
bin/mars --help
```

## Authorship & Copyright

LaRA 2 is being developed by [Jörg Winkler](mailto:j.winkler@fu-berlin.de), but it incorporates a lot of work
from other members of the [SeqAn project](http://www.seqan.de).

You can ask questions and report bugs on the [GitHub tracker](https://github.com/seqan/mars/issues).
Please also [subscribe](https://github.com/seqan/mars/subscription) and/or star us!
You can also follow SeqAn on [Twitter](https://twitter.com/SeqAnLib) to receive updates on MaRs.
