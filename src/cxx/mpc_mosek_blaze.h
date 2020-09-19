//
// Created by AliGriv on 2020-09-19.
//

#ifndef MEBTEST_MPC_MOSEK_BLAZE_H
#define MEBTEST_MPC_MOSEK_BLAZE_H


#include <blaze/Math.h>
#include <vector>
#include <mosek.h>
#include "mpc_mosek_utility_blaze.h"

using blaze::DynamicMatrix;
using blaze::StaticMatrix;
using blaze::rowMajor;
using blaze::columnMajor;
using blaze::StaticVector;
using blaze::DynamicVector;
using blaze::columnVector;
using blaze::rowVector;
using blaze::zero;
using blaze::uniform;

class MPCMosekErrorBlaze : public std::exception {
private:
    std::string message;
public:
    explicit MPCMosekErrorBlaze(MSKrescodee res) {
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
    ~MPCMosekErrorBlaze() throw() {}
};





class MPC_Mosek_xN_blaze {
private:
    int dim;
    std::vector <double> dT; /*containing discrete time intervals*/
    double T_h; /* time horizon */
    int m;
    int num_var;
    /* Dynamic-related matrices */
    std::vector <DynamicMatrix <double, columnMajor>> A_d;
    std::vector <DynamicMatrix <double, columnMajor>> B_d;
    DynamicMatrix <double, columnMajor> A_und;
    DynamicMatrix <double, columnMajor> B_und;
    DynamicVector <double, columnVector> u_max;
    DynamicVector <double, columnVector> u_min;
    /* Initial State */
    DynamicVector <double, columnVector> x0;
    /* blc <= A u <= buc */
    /* matrices that are actually sub-matrices of the whole big A, blc and buc */
    /* Limits on actuation */
    DynamicMatrix <double, columnMajor> A1;
    DynamicVector <double, columnVector> l1;
    DynamicVector <double, columnVector> u1;
    bool is_u_lim_set; /* says that all the coefficients for u has been made and the values of u_max and u_min are injected! */
    /* Continuity of actuators' commands */
    DynamicMatrix <double, columnMajor> A2;
    DynamicVector <double, columnVector> l2;
    DynamicVector <double, columnVector> u2;
    bool is_continuity_set;
    /* Convex Region constraint */
    DynamicMatrix <double, columnMajor> Ac_und;
    DynamicVector <double, columnVector> bc_und;
    DynamicMatrix <double, columnMajor> A3;
//    DynamicVector <double, columnVector> l3; /*this one is all infinity, so should be an array */
//    MSKrealt *l3;
    DynamicVector <double, columnVector> u3;
    int num_cvx_const;
    bool is_cvx_const_update;
    /* Objective Matrices: 0.5 x^T Qo x + c^T x + cf */
    DynamicMatrix <double, columnMajor> Q1;
    DynamicVector <double, columnVector> c1;

    DynamicMatrix <double, columnMajor> I_v_mask; /* extracts velocity  out of states!*/
    DynamicVector <double, columnVector> V_h;

    DynamicMatrix <double, columnMajor> Q2; /*uRu*/

    DynamicMatrix <double, columnMajor> Q3; /* for xN*/
    DynamicVector <double, columnVector> c3;

    DynamicMatrix <double, columnMajor> Qo; /* total */
    DynamicVector <double, columnVector> c;
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
    DynamicVector <double, columnVector> xN;
    DynamicMatrix <double, columnMajor> mpc_sol_u_k;
    DynamicMatrix <double, columnMajor> mpc_sol_x_k;
    MPC_Mosek_xN_blaze(int dim_val, int k_steps_val = 4, double t_h = 1.0, double mass_val = 1):
            dim(dim_val),
            T_h(t_h),
            mass(mass_val),
            xN(zero<double,columnVector>(dim_val)),
            x0(zero<double,columnVector>(2*dim_val)),
            k_steps(k_steps_val)
    {
        create_dT(t_h);
//            std::cout << "dT created" << std::endl;
        m = dT.size();
        num_var = dim*m;
        construct_dynamic_related_matrices();
//            std::cout << "dynamic matrices created" << std::endl;
        u_max = zero <double, columnVector>(dim_val);
        u_min = zero <double, columnVector>(dim_val);
        set_u_constraints(u_max, u_min);
        is_u_lim_set = false;
        is_continuity_set = false;
        set_u_continuity_constraint();
//            std::cout << "continuity matrices created" << std::endl;
        Ac_und = zero<double, columnMajor> (m*MAX_NUM_CVX_CONS_BLAZE, 2*num_var);
        bc_und = zero <double, columnVector>(m*MAX_NUM_CVX_CONS_BLAZE);
        A3 = zero<double, columnMajor> (m*MAX_NUM_CVX_CONS_BLAZE, num_var);
        u3 = zero <double, columnVector>(m*MAX_NUM_CVX_CONS_BLAZE);

        Q1 = zero<double, columnMajor> (num_var, num_var);
        c1 = zero <double, columnVector>(num_var);

        V_h = zero <double, columnVector>(m*dim);

        Q2 = zero<double, columnMajor> (num_var, num_var);
        Q2 = 2*declid(Q2);

        Q3 = zero<double, columnMajor> (num_var, num_var);
        c3 = zero <double, columnVector>(num_var);
        create_utility_matrices(dim_val, m);
//            std::cout << "utility matrices created" << std::endl;
        Qo = zero<double, columnMajor> (num_var, num_var);
        c = zero <double, columnVector>(num_var);
        cf = 0.0;
        create_msk_const_indices(dim_val, m);
//            std::cout << "MSK const indices created" << std::endl;
        mpc_sol_u_k = zero<double, columnMajor> (dim_val, k_steps_val);
        mpc_sol_x_k = zero<double, columnMajor> (2*dim_val, k_steps_val);

    }
    void create_utility_matrices(const int &dim, const int &m);
    void create_dT(const double t_h);
    std::vector <double> get_dT() const;
    void construct_dynamic_related_matrices();
    void set_u_constraints(DynamicVector <double, columnVector> u_max1, DynamicVector <double, columnVector> u_min1); /* maximum-minimum values for virtual forces */
    void set_u_continuity_constraint();
    void set_x0(const DynamicVector <double, columnVector> &init_pos, const DynamicVector <double, columnVector> &init_vel);
    void reset_cvx_constraint();
    void update_cvx_constraints(const DynamicMatrix <double, columnMajor> &A_c, const DynamicVector <double, columnVector> &b_c);

    void update_obj_matrices(const DynamicVector <double, columnVector> &v_h, double w1 = 0.1, double w2 = 0.5, double w3 = 4);
//    void update_obj_matrices(const DynamicVector <double, columnVector> &v_h, double w1 = 0, double w2 = 0, double w3 = 1.0);
    void create_msk_const_indices(const int &dim, const int &m);
    void set_k_steps(int k_step);
    int get_k_steps() const;
    bool solve_mpc(MSKenv_t *existing_en=NULL);

    bool do_mpc_mosek_xN(const DynamicVector <double, columnVector> &initial_position,
                         const DynamicVector <double, columnVector> &initial_velocity,
                         const DynamicVector <double, columnVector> &next_way_pt,
                         const DynamicVector <double, columnVector> &v_h,
                         const DynamicMatrix <double, columnMajor> &A_c,
                         const DynamicVector <double, columnVector> &b_c);

    ~MPC_Mosek_xN_blaze()
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

#endif //MEBTEST_MPC_MOSEK_BLAZE_H
