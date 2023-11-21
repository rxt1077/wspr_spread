FROM ubuntu:22.04
RUN mkdir /output
WORKDIR /wsjtx
RUN apt-get -y update
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install curl build-essential gcc \
    g++ gfortran cmake git asciidoc asciidoctor texinfo qtmultimedia5-dev \
    libqt5serialport5-dev qttools5-dev qttools5-dev-tools libboost-all-dev \
    libfftw3-dev libreadline-dev libusb-1.0-0-dev libudev-dev portaudio19-dev \
    dos2unix
RUN curl https://gigenet.dl.sourceforge.net/project/wsjt/wsjtx-2.6.1/wsjtx-2.6.1.tgz | tar zx --strip-components=1
COPY all.patch wsjtx.patch
WORKDIR /build
RUN cmake ../wsjtx
RUN cmake --build .
RUN cmake --build . --target wsjtx-update
RUN dos2unix wsjtx-prefix/src/wsjtx/lib/wsprd/wsprd.c
