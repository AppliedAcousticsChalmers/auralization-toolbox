# The Chalmers Auralization Toolbox

![graphical_abstract](resources/graphical_abstract.png "graphical_abstract")

## Overview

The Chalmers Auralization Toolbox is a collection of MATLAB scripts that auralize sampled sound fields. The sound field quantity that is sampled can be either the sound pressure or the combination of sound pressure and particle velocity (or equivalently the sound pressure gradient). The primary purpose is auralization of data simulated with methods like FDTD, FEM, BEM, and the like.

**Binaural audio examples** that were created with the toolbox are available [here](http://www.ta.chalmers.se/research/audio-technology-group/audio-examples/auditorium-acoustics-2023/). 

The figure below illustrates the types of sampling grids that can be processed (volumetric, spherical surfaces, and cubical surfaces). 

![grids](resources/grids.png "grids")

The sampled data can be converted to either an ambisonic representation, which can then be rendered both binaurally and using loudspeaker arrays. The conventions that we use (i.e. N3D and ACN) are compatible with software tools like [SPARTA](https://leomccormack.github.io/sparta-site/) and the [IEM Plugin Suite](https://plugins.iem.at/). Or, the sampled data can be rendered binaurally without an intermediate format (*direct rendering*).

You'll need to download the employed HRIRs from [here](https://zenodo.org/record/3928297/files/HRIR_L2702.sofa?download=1) and store them in the subfolder `hrtfs` (The MATLAB script is going to do that automatically for you, both the downloading and creating that folder.) as well as the SOFA MATLAB API from [here](https://sourceforge.net/projects/sofacoustics/).

The toolbox was originally presented in 

> Jens Ahrens, "A Software Tool for Auralization of Simulated Sound Fields," Auditorium Acoustics, Athens, Greece, 2023

The work is on-going, so please revisit this repository regularly for updates. Note that some of the implementations are in an experimental state. Make sure that you get in touch with us at jens.ahrens@chalmers.se if things behave differently than you were expecting.

## Example Data

In the folder `resources`, we provide example data for a reverberant room with a reverb decay time around 1 s ("big hall") as well as a much dryer living-room type room ("living room"). 

## Upcoming features

Currently, all methods use sound pressure as input data Those methods that employ both sound pressure and particle velocity obtain the particle velocity from two adjacent pressure sampling points. We will allow for incorporating particle velocity direct as some simulation frameworks provide it.

The example grids comprise approx. 300 sampling points. This produces good results, but the auralization is not perfectly transparent meaning that there are cases in which the auralization can be distinguished from the ground truth (refer to the audio examples above). Increasing the number of sampling point will eventually overcome the deviations. Yet, we are working on updates on the methods that can potentially yield perceptually transparent auralization with the present or even lower number of sampling points. 

We will provide examples that demonstrate how data from the most common room acoustic simulation softwares can be auralized with our toolbox. 

## References

> B. Støfringsdal and U.P. Svensson, 'Conversion of Discretely Sampled Sound Field Data to Aualization Formats,' J. Audio. Eng. Soc. 54(5), pp. 380-400 (May 2006)

> J. Sheaffer, M. van Walstijn, B. Rafaely, and K. Kowalczyk, 'Binaural Reproduction of Finite Difference Simulations Using Spherical Array Processing', IEEE/ACM TASLP 23(12), pp. 2125-2135 (Dec. 2015)

> M. A. Poletti and U. P. Svensson, 'Beamforming Synthesis of Binaural Responses From Computer Simulations of Acoustic Spaces,' J. Acoust. Soc. Am. 124, pp. 301–315 (2008)

