# wspr_spread

This repository contains a Python program for calculating Doppler spread from WSPR data as well as a patch, spread.patch, for wsprd that does the same.
You can find more information in [this article](https://using.tech/posts/wspr-spread/).
It also contains a patch, nodrift.patch, to keep wsprd from performing a frequency drift search.
Finally a third patch, all.patch, adds both features to wsprd.


## Building wsjtx with this patch

wsjt-x supports applying patches and creating packages in its CMake build system.
See the Dockerfile in this directory which creates a build environment with the all patch.
You can then use it to build a Debian package.
The `make_deb.sh` script automates the whole process. 
