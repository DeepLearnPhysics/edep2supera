# edep2supera

C++/Python software to run `Supera` ([repository](https://github.com/DeepLearnPhysics/SuperaAtomic)), a library that can generate true labels for a ML based data reconstruction chain called `lartpc_mlreco3d` ([repository](https://github.com/DeepLearnPhysics/lartpc_mlreco3d)). This repository is dedicated to interface [EDepSim](https://github.com/ClarkMcGrew/edep-sim) and [larcv](https://github.com/DeepLearnPhysics/larcv2) frameworks that input and output information in data files respectively. 

## Python binding

Executable scripts are in Python and C++ classes and functions are bound to Python via [PyROOT](https://root.cern/manual/python/) following a pattern used by EDepSim (as opposed to Supera which uses [pybind11](https://pybind11.readthedocs.io/en/stable/)). This allows both edep2supera and EDepSim methods to access C++ objects within Python in a common manner.

## How to build and install

1. Use [this docker image](https://hub.docker.com/layers/deeplearnphysics/dune-nd-sim/ub20.04-cpubase-edep2supera/images/sha256-a9f01459ec3f67c85e2bea813d14d1069ca41a5d8f4b8f301336b992de238ca2?context=repo) (or any other container image derived from it) to get the _most_ of required softwares. Only additional item you would need is [scikit-build](https://scikit-build.readthedocs.io/en/latest/skbuild.html) which you can install with `pip`. 
2. Clone this repository (below) or fork-and-clone.
```
> git clone https://github.com/DeepLearnPhysics/edep2supera
```
3. Build+Install (the example below installs under your `$HOME/.local` path)
```
> pip install . --user
```

## Software validation (unit test)
Simply try:
```
pytest
```
at the top-level directory. 

## How to contribute

1. Fork this repository to your personal github account.
2. Clone the repository to your local machine. Follow the build/install instruction above and make sure you can set up.
3. Create your branch to contain your own development.  Code code code.
4. When it's ready to be shared, make sure a unit test passes. Then request to merge by sending a pull request.

**Optional but strongly recommended**: implement a unit test for the added component of your code so that we can reduce chance of someone else breaking in future development.

## Communication for development
* Feel free to use github issues! We try to be attentive as much as possible.
* Join our weekly ND Technical Software Meeting (noon PST every Wednesday, here's [indico category](https://indico.slac.stanford.edu/category/23/)).
* Join our mailing list dunend-simreco-technical@listserv.slac.stanford.edu
  - Send an email to listserv@listserv.slac.stanford.edu with an empty title. The body should include this text: `SUBSCRIBE DUNEND-SIMRECO-TECHNICAL FIRST LAST` where you should replace the `FIRST` and `LAST` with your first and last names respectively.
* Join our slack channel (contact [Kazu](mailto:kterao@slac.stanford.edu)).

## Status

For RooTracker-backed EDepSim input, the converter also extracts the incoming
neutrino, target, reaction mode, outgoing lepton, and derived lepton kinematics
from each `TG4PrimaryVertex`. These are written as the LArCV
`neutrino_mc_truth` product. A generator must preserve an `initial-state`
informational vertex in EDepSim for this product to be populated.
Neutrino interaction IDs use the same zero-based, event-local ordering as the
Supera particle labels so downstream readers can associate the two products.

Milestones

- [x] Fill Supera datatypes from EDepSim datatypes (i.e. TG4Event attributes)
- [x] Implement particle hierarchy tracking algorithm from [old Supera](https://github.com/DeepLearnPhysics/Supera)
- [x] Implement a python run script to fill input data and run Supera
- [x] Implement a python run script to drive larcv to create an output file and store information from Supera
