# BlackHat ntuple analysis
* Author: Ivan Pogrebnyak, [MSU](http://www.pa.msu.edu/), ATLAS
* Repository: https://github.com/ivankp/bh_analysis
* Source (zip): https://github.com/ivankp/bh_analysis/archive/master.zip

## Usage

Compiled executables are created in the `bin` directory. Help message can be displayed by running `./bin/executable -h`, or `./bin/executable --help`, or, in most cases, simply `./bin/executable` without any arguments.

Below is the list of executables, their purpose, and usage examples. The list is in the order in which the programs need to be run to perform the analysis and obtain plots.

### reweigh
* Purpose: Reweighting of BlackHat ntuples.
* Output: A root ntuple with only new weights.
* Usage example: `./bin/reweigh --bh=born_bh.root -o bort_weights.root`
* Not yet implemented: arguments for selecting scales

### hist_foo
* Purpose: This this the analysis program. It produces plots for different weights.
* Output: A root file with histograms in directories corresponding to selected weights.
* Usage examples:
  * Minimal: <br />
    `./bin/hist_foo --bh=born_bh.root -o bort_hist.root` <br />
    Produces histograms for the weight in the BlackHat ntuple; clustering is done by FastJet, AntiKt4 algorithm is the default.
  * Use new weights: <br />
    `./bin/hist_foo --bh=born_bh.root --wt=born_weights.root -o bort_hist.root` <br />
    Clustering is done by FastJet
  * Use SpartyJet ntuples: <br />
    `./bin/hist_foo --bh=born_bh.root --sj=born_sj.root --wt=born_weights.root -o bort_hist.root`

Note: Numbers of entries in histograms are not numbers of events, but numbers of ntuple entries. These are not the same for real ntuples.

### merge_parts
* Purpose: Merge together different kinds of ntuples (born, real, integrated-subtraction, virtual).
* Output: Root file in the same format with merged histograms.
* Usage example: `./bin/merge_parts NLO.root B.root RS.root I.root V.root`

Note: Numbers of entries in histograms are also combined.

### plot
* Purpose: Plot histograms with scale variation and PDF uncertainty bands.
* Input: A single root file with histograms in directories corresponding to scale and PDF variations.
* Output: A single pdf file with plots on multiple pages.
* Usage example: `./bin/plot NLO.root` will output NLO.pdf. `-o` flag is supported.

---

## Compilation
* To compile for the first time: `make tools parts all`
* To recompile: `make`

### Requirements
The code is written in C++11 and requires a C++ compiler which supports the -std=c++11 flag; -std=c++0x is insufficient.

#### Other requirements:
* ROOT: <https://root.cern.ch/drupal/> -- either version 5 or 6 is ok
* BOOST: <http://www.boost.org/> -- version 1.42.0 and above
* LHAPDF6: <https://lhapdf.hepforge.org/> -- version 5 is not supported
* FastJet: <http://fastjet.fr/>

---

## Making changes

### Adding a new analysis
* Make a copy of `src/hist_H3j`, e.g. <br />
  `cp src/hist_H3j src/hist_photon3j`
* Add the new analysis to the `Makefile` by analogy with `src/hist_H2j` and `src/hist_H3j`
* If you use GitHub, please submit a pull request, so that your code can be incorporated into the repository for the benefit of others and maintanance.
