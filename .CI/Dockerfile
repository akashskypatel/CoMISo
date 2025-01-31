FROM docker.io/debian:trixie-slim
MAINTAINER Martin Heistermann <martin.heistermann@inf.unibe.ch>

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    binutils \
    ca-certificates \
    coinor-libipopt-dev \
    clang \
    cmake \
    curl \
    git \
    libc-dev \
    libcgal-dev \
    liblapack-dev \
    libeigen3-dev \
    libmetis-dev \
    libmumps-seq-5.7 \
    libomp-dev \
    libopenblas64-serial-dev \
    libsuitesparse-dev \
    libtool \
    locales \
    ninja-build \
    time \
    tzdata \
    wget

COPY . /app/CoMISo

RUN cd /app/CoMISo && \
    cmake -B build -G Ninja \
    -D CMAKE_C_COMPILER=clang \
    -D CMAKE_CXX_COMPILER=clang++ \
    -D CMAKE_BUILD_TYPE=Debug \
    -D CMAKE_CXX_FLAGS="-fdiagnostics-color=always" \
    . && CLICOLOR_FORCE=1 cmake --build build
