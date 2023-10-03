FROM ubuntu:20.04

# build tools
RUN apt-get -y update
RUN apt-get install -y \
    wget \
    libssl-dev \
    git \
    build-essential


# install cmake
RUN wget https://github.com/Kitware/CMake/archive/refs/tags/v3.20.2.tar.gz -O cmake.tar.gz
RUN tar zxvf cmake.tar.gz
RUN cd CMake-3.20.2 && env CC=gcc ./bootstrap && make -j4 && make install && cd ..
RUN hash -r
RUN rm -rf CMake-3.20.2

# install yamp-cpp
RUN git clone https://github.com/jbeder/yaml-cpp.git
RUN cd yaml-cpp && mkdir build && cd build && cmake .. && make install && cd ../..
RUN rm -rf yaml-cpp

# install Eigen
RUN ln -sf /usr/share/zoneinfo/Asia/Tokyo /etc/localtime
RUN apt-get install -y libeigen3-dev=3.3.7-2

# install g2o
RUN git clone https://github.com/RainerKuemmerle/g2o.git
RUN cd g2o && mkdir build && cd build && cmake .. && make install && cd ../..
RUN rm -rf g2o

# instal python3
RUN apt-get install -y \
    python3=3.8.2-0ubuntu2 \
    python3-pip

# install python packages
RUN python3 -m pip install \
    matplotlib==3.5.0 \
    numpy==1.21.4

# install OpenGL
RUN apt-get install -y pkg-config \
    mesa-utils \
    libglu1-mesa-dev \
    freeglut3-dev \
    mesa-common-dev \
    libglew-dev \
    libglfw3-dev \
    libglm-dev \
    libao-dev \
    libmpg123-dev
RUN git clone https://github.com/glfw/glfw.git
RUN apt-get install -y libxinerama-dev \
    libxcursor-dev \
    libxi-dev
RUN cd glfw && mkdir build && cd build && cmake .. && make install && cd ../..
RUN rm -rf glfw

# install ceres
RUN apt-get install -y libgoogle-glog-dev libgflags-dev libatlas-base-dev libsuitesparse-dev
RUN wget http://ceres-solver.org/ceres-solver-2.1.0.tar.gz
RUN tar zxf ceres-solver-2.1.0.tar.gz
RUN mkdir ceres-bin && cd ceres-bin && cmake ../ceres-solver-2.1.0 && make -j3 && make install && cd ..
RUN rm ceres-solver-2.1.0.tar.gz && rm -rf ceres-bin

# install boost
RUN apt-get install -y libboost-dev=1.71.0.0ubuntu2

# install pytorch
RUN python3 -m pip install \
    torch==1.13.1

# clear cache
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf ceres-solvers-2.1.0 && rm cmake.tar.gz

# install python packages
RUN python3 -m pip install \
    pyyaml==5.3.1 \
    tensorboard==2.12.0