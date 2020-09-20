# Mosek-Eigen-Blaze-test
In this small project, it is tried to compare the performance of algebraic-operations such as matrix multiplication and the final result on a simple optimization problem is to be evaluated.


This file is yet an ongoing project :)
It has a few requirements:
```
sudo apt-get install gfortran libblas-dev liblapack-dev
sudo apt install libeigen3-dev
sudo apt-get install libboost-all-dev
git clone https://bitbucket.org/blaze-lib/blaze.git
cd blaze/
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/
sudo make install
```
The initial result is that for a particular function (method) in MPC_Mosek_xN class, the calculations with blast-lib takes almost 3~4 times longer than Eigen-lib.
```
For instance let's measure update_obj_matrices(v_h) method by running it over 100 times which uses Eigen3 Library
Using Eigen library, on average basis it took 0.0234767 seconds
For instance let's measure update_obj_matrices(v_h) method by running it over 100 times which uses blaze library
Using Eigen library, on average basis it took 0.085581 seconds
```
