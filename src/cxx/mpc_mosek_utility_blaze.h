//
// Created by griv on 2020-09-19.
//

#ifndef MEBTEST_MPC_MOSEK_UTILITY_BLAZE_H
#define MEBTEST_MPC_MOSEK_UTILITY_BLAZE_H
#include <blaze/Math.h>
#include <vector>

using blaze::DynamicMatrix;
using blaze::StaticMatrix;
using blaze::rowMajor;
using blaze::columnMajor;
/* Auxiliary functions and variables */

const int MAX_NUM_CVX_CONS_BLAZE = 30;


static const DynamicMatrix <double, columnMajor> I_p_3d_blaze = [] {
    DynamicMatrix <double, columnMajor> temp(3,6);
    temp = { {1, 0, 0, 0, 0, 0},
             {0, 0, 1, 0, 0, 0},
             {0, 0, 0, 0, 1, 0}};
    return temp;
}();
static const DynamicMatrix <double, columnMajor> I_p_2d_blaze = [] {
    DynamicMatrix <double, columnMajor> temp(2,4);
    temp = { {1, 0, 0, 0},
             {0, 0, 1, 0}};
    return temp;
}();

const DynamicMatrix <double, columnMajor> &I_p_blaze(int dim){
    if (dim == 2)
        return I_p_2d_blaze;
    else if (dim == 3)
        return I_p_3d_blaze;
};
static const DynamicMatrix <double, columnMajor> I_v_3d_blaze = [] {
    DynamicMatrix <double, columnMajor> temp(3,6);
    temp ={{0, 1, 0, 0, 0, 0},
           {0, 0, 0, 1, 0, 0},
           {0, 0, 0, 0, 0, 1}};
    return temp;
}();
static const DynamicMatrix <double, columnMajor> I_v_2d_blaze = [] {
    DynamicMatrix <double, columnMajor> temp(2,4);
    temp = {{ 0, 1, 0, 0},
            {0, 0, 0, 1}};
    return temp;
}();

const DynamicMatrix <double, columnMajor> &I_v_blaze(int dim){
    if (dim == 2)
        return I_v_2d_blaze;
    else if (dim == 3)
        return I_v_3d_blaze;
};
#endif //MEBTEST_MPC_MOSEK_UTILITY_BLAZE_H
