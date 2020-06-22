# MaRs: Motif-based aligned RNA searcher

[![Build Status](https://travis-ci.com/joergi-w/mars.svg?branch=master)](https://travis-ci.com/joergi-w/mars)

Mars is a tool that reads a structural multiple RNA alignment (e.g. from LaRA) and derives fuzzy stem loop descriptors
from it. These descriptors are then subject to a search in a genomic database and Mars returns the hits where the
RNA structure is found, accompanied with a quality value for each hit.

Installation:
1. clone this repository: `git clone --recurse-submodules https://github.com/joergi-w/mars.git`
2. create a build directory and visit it: `mkdir build && cd build`
3. run cmake: `cmake ../mars`
4. build the application: `make`
5. optional: build and run the tests: `make test`
6. optional: build the api documentation: `make doc`
7. execute the app: `./mars`
