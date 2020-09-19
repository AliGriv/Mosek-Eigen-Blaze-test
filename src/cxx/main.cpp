#include <iostream>
#include <blaze/Math.h>
#include <Eigen/Core>
#include "mpc_mosek_eigen.h"
#include "mpc_mosek_utility_eigen.h"
#include "mpc_mosek_blaze.h"
#include "mpc_mosek_utility_blaze.h"
#include <chrono>
#include <vector>
//
using blaze::DynamicMatrix;
using blaze::StaticMatrix;
using blaze::rowMajor;
using blaze::columnMajor;
using blaze::StaticVector;
using blaze::DynamicVector;
using blaze::columnVector;
using blaze::rowVector;
using blaze::zero;
int main()
{
    std::cout << "inside main function" << std::endl;
    MPC_Mosek_xN test_mpc(3);
    std::vector <float> eigen_time_vec;
    float avg_eigen_time {0.0};
    std::cout << "For instance let's measure update_obj_matrices(v_h) method by running it over 100 times which uses Eigen3 Library" << std::endl;
    Eigen::VectorXd v_h(3);
//    v_h << 3, 0, 0;
    for (int i {0}; i < 100; ++i)
    {
        v_h = 0.5*Eigen::VectorXd::Random(3);
        auto begin1 = std::chrono::high_resolution_clock::now();
        test_mpc.update_obj_matrices(v_h);
        auto end1 = std::chrono::high_resolution_clock::now();
        auto elapsed1 = std::chrono::duration_cast < std::chrono::duration < float >> (end1 - begin1);
        eigen_time_vec.push_back(elapsed1.count());
    }
    avg_eigen_time = std::accumulate( eigen_time_vec.begin(), eigen_time_vec.end(), 0.0)/eigen_time_vec.size();
    std::cout << "Using Eigen library, on average basis it took " << avg_eigen_time << " seconds" << std::endl;


//    std::cout << "creating mpc_mosek_xN_blaze object" <<std::endl;
    MPC_Mosek_xN_blaze test_mpc_blaze(3);
//    std::cout << "mpc_mosek_xN_blaze object created" <<std::endl;
    std::vector <float> blaze_time_vec;
    float avg_blaze_time {0.0};
    std::cout << "For instance let's measure update_obj_matrices(v_h) method by running it over 100 times which uses blaze library" << std::endl;

//    v_h_blaze = {3.0, 0.0, 0.0};
    for(int i {0}; i < 100; ++i)
    {
        DynamicVector <double> v_h_blaze(3);
        v_h_blaze = {rand()/double(RAND_MAX), rand()/double(RAND_MAX), rand()/double(RAND_MAX)};
        auto begin2 = std::chrono::high_resolution_clock::now();
        test_mpc_blaze.update_obj_matrices(v_h_blaze);
        auto end2 = std::chrono::high_resolution_clock::now();
        auto elapsed2 = std::chrono::duration_cast < std::chrono::duration < float >> (end2 - begin2);
        blaze_time_vec.push_back(elapsed2.count());
    }
    avg_blaze_time = std::accumulate( blaze_time_vec.begin(), blaze_time_vec.end(), 0.0)/blaze_time_vec.size();
    std::cout << "Using Eigen library, on average basis it took " << avg_blaze_time << " seconds" << std::endl;


    std::cout << "Program terminating" << std::endl;
    return 0;
}