#!/bin/bash

# make sure the image is built and tagged
docker build -t wsjtx .

# run the build command and copy the DEB file to the bind-mounted dir
docker run -v $(pwd):/output wsjtx /bin/bash -c "cmake --build . --target package; cp wsjtx-prefix/src/wsjtx-build/wsjtx_2.6.1_amd64.deb /output/"
