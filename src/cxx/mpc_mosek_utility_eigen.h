//
// Created by AliGriv on 2020-09-18.
//

#ifndef TEST_BLAZE_MPC_MOSEK_UTILITY_EIGEN_H
#define TEST_BLAZE_MPC_MOSEK_UTILITY_EIGEN_H

#include <Eigen/Core>
#include <vector>

/* Auxiliary functions and variables */

const int MAX_NUM_CVX_CONS = 30;


static const Eigen::MatrixXd I_p_3d = [] {
    Eigen::MatrixXd temp(3,6);
    temp << 1, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 1, 0;
    return temp;
}();
static const Eigen::MatrixXd I_p_2d = [] {
    Eigen::MatrixXd temp(2,4);
    temp << 1, 0, 0, 0,
            0, 0, 1, 0;
    return temp;
}();

const Eigen::MatrixXd &I_p(int dim){
    if (dim == 2)
        return I_p_2d;
    else if (dim == 3)
        return I_p_3d;
};
static const Eigen::MatrixXd I_v_3d = [] {
    Eigen::MatrixXd temp(3,6);
    temp << 0, 1, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1;
    return temp;
}();
static const Eigen::MatrixXd I_v_2d = [] {
    Eigen::MatrixXd temp(2,4);
    temp << 0, 1, 0, 0,
            0, 0, 0, 1;
    return temp;
}();

const Eigen::MatrixXd &I_v(int dim){
    if (dim == 2)
        return I_v_2d;
    else if (dim == 3)
        return I_v_3d;
};



#endif //TEST_BLAZE_MPC_MOSEK_UTILITY_EIGEN_H
