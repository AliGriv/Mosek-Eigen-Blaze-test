#include <iostream>
#include <blaze/Math.h>
#include <Eigen/Core>
#include "mpc_mosek_eigen.h"
#include "mpc_mosek_utility_eigen.h"
//
using blaze::StaticVector;
using blaze::DynamicVector;

int main()
{
    std::cout << "inside main function" << std::endl;
    MPC_Mosek_xN test_mpc(3);
    // Instantiation of a static 3D column vector. The vector is directly initialized as
    //    ( 4 -2  5 )
    StaticVector<int,3UL> a{ 4, -2, 5 };

    // Instantiation of a dynamic 3D column vector. Via the subscript operator the values are set to
    //    ( 2  5 -3 )
    DynamicVector<int> b( 3UL );
    b[0] = 2;
    b[1] = 5;
    b[2] = -3;

    // Adding the vectors a and b
    DynamicVector<int> c = a + b;

    // Printing the result of the vector addition
    std::cout << "c =\n" << c << "\n";
//    std::cout << "Program terminating" << std::endl;
//    return 0;
}