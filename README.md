# BlackHat ntuple analysis

## Usage

The compiled executables are created in the `bin` directory. Help message can be displayed by running `./bin/executable -h`, or `./bin/executable --help`, or, in most cases, simply `./bin/executable` without any arguments.

Listed below is the list of executables, their purpose, and usage examples.

### reweigh
* Purpose: Reweighting of BlackHat ntuples
* Output: A root ntuple with only new weights
* Usage example: `./bin/reweigh --bh=born_bh.root -o bort_weights.root`

### hist_H3j
* Purpose: Produce plots for different weights
* Output: A root file with histograms in directories corresponding to selected weights
* Usage examples:
  * Minimal: <br />
    `./bin/hist_H3j --bh=born_bh.root -o bort_hist.root` <br />
    Produces histograms for the weight in the BlackHat ntuple; clustering is done by FastJet, AntiKt4 algorithm is the default.
  * Use new weights: <br />
    `./bin/hist_H3j --bh=born_bh.root --wt=born_weights.root -o bort_hist.root` <br />
    Clustering is done by FastJet
  * Use SpartyJet ntuples: <br />
    `./bin/hist_H3j --bh=born_bh.root --sj=born_sj.root --wt=born_weights.root -o bort_hist.root`

---

## Compilation
* To compile for the first time: `make tools parts all`
* To recompile: `make`

### Requirements
The code is written in C++11 and therefore requires a C++ compiler which supports the -std=c++11 flag.
#### Other requirements:
* ROOT: <https://root.cern.ch/drupal/>
* BOOST: <http://www.boost.org/>
* FastJet: <http://fastjet.fr/>

---

## Making changes

### Adding a new analysis
* Make a copy of `src/hist_H3j`, e.g. <br />
  `cp src/hist_H3j src/hist_photon3j`
* Add the new analysis to the `Makefile` by analogy with `src/hist_H2j` and `src/hist_H3j`
* If you use GitHub, submit a pull request, so that your code can be incorporated into the repository for the benefit of others and maintanance.
