//
// Created by AliGriv on 2020-09-18.
//

#ifndef TEST_BLAZE_MPC_MOSEK_EIGEN_H
#define TEST_BLAZE_MPC_MOSEK_EIGEN_H


#include <Eigen/Core>
#include <vector>
#include <mosek.h>
#include "mpc_mosek_utility_eigen.h"


class MPCMosekError : public std::exception {
private:
    std::string message;
public:
    explicit MPCMosekError(MSKrescodee res) {
        /* In case of an error print error code and description. */
        char symname[MSK_MAX_STR_LEN];
        char desc[MSK_MAX_STR_LEN];
        MSK_getcodedesc (res,
                         symname,
                         desc);
        message = std::string(symname) + ": " + std::string(desc);
    }

    const char * what () const throw () {
        return message.c_str();
    }
    ~MPCMosekError() throw() {}
};





class MPC_Mosek_xN {
private:
    int dim;
    std::vector <double> dT; /*containing discrete time intervals*/
    double T_h; /* time horizon */
    int m;
    int num_var;
    /* Dynamic-related matrices */
    std::vector <Eigen::MatrixXd> A_d;
    std::vector <Eigen::MatrixXd> B_d;
    Eigen::MatrixXd A_und;
    Eigen::MatrixXd B_und;
    Eigen::VectorXd u_max;
    Eigen::VectorXd u_min;
    /* Initial State */
    Eigen::VectorXd x0;
    /* blc <= A u <= buc */
    /* matrices that are actually sub-matrices of the whole big A, blc and buc */
    /* Limits on actuation */
    Eigen::MatrixXd A1;
    Eigen::VectorXd l1;
    Eigen::VectorXd u1;
    bool is_u_lim_set; /* says that all the coefficients for u has been made and the values of u_max and u_min are injected! */
    /* Continuity of actuators' commands */
    Eigen::MatrixXd A2;
    Eigen::VectorXd l2;
    Eigen::VectorXd u2;
    bool is_continuity_set;
    /* Convex Region constraint */
    Eigen::MatrixXd Ac_und;
    Eigen::VectorXd bc_und;
    Eigen::MatrixXd A3;
//    Eigen::VectorXd l3; /*this one is all infinity, so should be an array */
//    MSKrealt *l3;
    Eigen::VectorXd u3;
    int num_cvx_const;
    bool is_cvx_const_update;
    /* Objective Matrices: 0.5 x^T Qo x + c^T x + cf */
    Eigen::MatrixXd Q1;
    Eigen::VectorXd c1;

    Eigen::MatrixXd I_v_mask; /* extracts velocity  out of states!*/
    Eigen::VectorXd V_h;

    Eigen::MatrixXd Q2; /*uRu*/

    Eigen::MatrixXd Q3; /* for xN*/
    Eigen::VectorXd c3;

    Eigen::MatrixXd Qo; /* total */
    Eigen::VectorXd c;
    double cf;

    /* mosek parameters */
    MSKint32t numqonz;
    MSKint32t *qosubi;
    MSKint32t *qosubj;
    MSKrealt *qovals;
    MSKint32t num_clist;
    MSKint32t *subj_clist;

    /* boundary on variables */
    MSKint32t num_var_bounds;
    MSKint32t *var_bounds_subi;
    MSKboundkeye *var_bound_keys;
    MSKrealt *blx;
    MSKrealt *bux;


    MSKint32t num_alist_1; /* number of coefficients that should be changed */
    MSKint32t *subi_alist_1; /* row indices */
    MSKint32t *subj_alist_1; /* col indices */
    MSKrealt *a_vals_1;
    MSKint32t num_bound_list_1; /* number of bounds that should be changed */
    MSKint32t *bound_sub_list_1; /* constraint indexes */
    MSKboundkeye *bound_key_list_1; /* keys for constraint */
    MSKrealt *blc_list_1; /* lower bound */
    MSKrealt *buc_list_1; /* upper bound */

    int k_steps; /* The numbers of discrete steps of MPC that will be applied */
    /* by default k_steps = 4 */



public:
    double mass;
    Eigen::VectorXd xN;
    Eigen::MatrixXd mpc_sol_u_k;
    Eigen::MatrixXd mpc_sol_x_k;
    MPC_Mosek_xN(int dim_val, int k_steps_val = 4, double t_h = 1.0, double mass_val = 1):
            dim(dim_val),
            T_h(t_h),
            mass(mass_val),
            xN(Eigen::VectorXd::Zero(dim_val)),
            x0(Eigen::VectorXd::Zero(2*dim_val)),
            k_steps(k_steps_val)
    {
        create_dT(t_h);
//            std::cout << "dT created" << std::endl;
        m = dT.size();
        num_var = dim*m;
        construct_dynamic_related_matrices();
//            std::cout << "dynamic matrices created" << std::endl;
        u_max = Eigen::VectorXd::Zero(dim_val);
        u_min = Eigen::VectorXd::Zero(dim_val);
        set_u_constraints(u_max, u_min);
        is_u_lim_set = false;
        is_continuity_set = false;
        set_u_continuity_constraint();
//            std::cout << "continuity matrices created" << std::endl;
        Ac_und = Eigen::MatrixXd::Zero(m*MAX_NUM_CVX_CONS, 2*num_var);
        bc_und = Eigen::VectorXd::Zero(m*MAX_NUM_CVX_CONS);
        A3 = Eigen::MatrixXd::Zero(m*MAX_NUM_CVX_CONS, num_var);
        u3 = Eigen::VectorXd::Zero(m*MAX_NUM_CVX_CONS);

        Q1 = Eigen::MatrixXd::Zero(num_var, num_var);
        c1 = Eigen::VectorXd::Zero(num_var);

        V_h = Eigen::VectorXd::Zero(m*dim);

        Q2 = 2*Eigen::MatrixXd::Identity(num_var, num_var);

        Q3 = Eigen::MatrixXd::Zero(num_var, num_var);
        c3 = Eigen::VectorXd::Zero(num_var);
        create_utility_matrices(dim_val, m);
//            std::cout << "utility matrices created" << std::endl;
        Qo = Eigen::MatrixXd::Zero(num_var, num_var);
        c = Eigen::VectorXd::Zero(num_var);
        cf = 0.0;
        create_msk_const_indices(dim_val, m);
//            std::cout << "MSK const indices created" << std::endl;
        mpc_sol_u_k = Eigen::MatrixXd::Zero(dim_val, k_steps_val);
        mpc_sol_x_k = Eigen::MatrixXd::Zero(2*dim_val, k_steps_val);

    }
    void create_utility_matrices(const int &dim, const int &m);
    void create_dT(const double t_h);
    std::vector <double> get_dT() const;
    void construct_dynamic_related_matrices();
    void set_u_constraints(Eigen::VectorXd u_max1, Eigen::VectorXd u_min1); /* maximum-minimum values for virtual forces */
    void set_u_continuity_constraint();
    void set_x0(const Eigen::VectorXd &init_pos, const Eigen::VectorXd &init_vel);
    void reset_cvx_constraint();
    void update_cvx_constraints(const Eigen::MatrixXd &A_c, const Eigen::VectorXd &b_c);

    void update_obj_matrices(const Eigen::VectorXd &v_h, double w1 = 0.1, double w2 = 0.5, double w3 = 4);
//    void update_obj_matrices(const Eigen::VectorXd &v_h, double w1 = 0, double w2 = 0, double w3 = 1.0);
    void create_msk_const_indices(const int &dim, const int &m);
    void set_k_steps(int k_step);
    int get_k_steps() const;
    bool solve_mpc(MSKenv_t *existing_en=NULL);

    bool do_mpc_mosek_xN(const Eigen::VectorXd &initial_position,
                         const Eigen::VectorXd &initial_velocity,
                         const Eigen::VectorXd &next_way_pt,
                         const Eigen::VectorXd &v_h,
                         const Eigen::MatrixXd &A_c,
                         const Eigen::VectorXd &b_c);

    ~MPC_Mosek_xN()
    {
        delete [] qosubi;
        delete [] qosubj;
        delete [] subj_clist;
        delete [] qovals;
        /* boundary on variables */
        delete [] var_bounds_subi;
        delete [] var_bound_keys;
        delete [] blx;
        delete [] bux;


        delete [] subi_alist_1; /* row indices */
        delete [] subj_alist_1; /* col indices */
        delete [] a_vals_1;
        delete [] bound_sub_list_1; /* constraint indexes */
        delete [] bound_key_list_1; /* keys for constraint */
        delete [] blc_list_1; /* lower bound */
        delete [] buc_list_1; /* upper bound */

    }
};


#endif //TEST_BLAZE_MPC_MOSEK_EIGEN_H
