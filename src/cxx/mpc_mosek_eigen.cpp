//
// Created by AliGriv on 2020-09-18.
//

#include "mpc_mosek_eigen.h"
#include <iostream>
#include <chrono>

#include "math.h"
#include "mpc_mosek_utility_eigen.h"

void check_res(MSKrescodee res) {
    if (res != MSK_RES_OK) {
        throw MPCMosekError(res);
    }
}
/* This function prints log output from MOSEK to the terminal. */
static void MSKAPI printstr(void *handle,
                            MSKCONST char str[])
{
    printf("%s",str);
} /* printstr */

void MPC_Mosek_xN::create_dT(const double t_h){
    double t {0.01};
    while (t <= t_h)
    {
        if (t <= 0.1)
        {
            this->dT.push_back(0.01);
            t += 0.01;
        }
        else if (t <= 0.2)
        {
            this->dT.push_back(0.02);
            t += 0.02;
        }
        else if (t <= 0.6)
        {
            this->dT.push_back(0.04);
            t += 0.04;
        }
        else if (t <= 1.0)
        {
            this->dT.push_back(0.08);
            t += 0.08;
        }
        else
        {
            this->dT.push_back(0.1);
            t += 0.1;
        }
    }
    if (t < t_h)
    {
        this->dT.push_back(t_h - t);
    }
}
std::vector <double> MPC_Mosek_xN::get_dT() const {
    return this->dT;
}
void MPC_Mosek_xN::construct_dynamic_related_matrices() {

    /* create A_d and B_d_temp */
    if (this->dim == 3)
    {
        this->A_d = std::vector <Eigen::MatrixXd> (this->m, Eigen::MatrixXd::Zero(2*(this->dim), 2*(this->dim)));
        this->B_d = std::vector <Eigen::MatrixXd> (this->m, Eigen::MatrixXd::Zero(2*(this->dim), (this->dim)));
        for (int i {0}; i < this->m ; ++i)
        {
            this->A_d.at(i) << 1, this->dT.at(i), 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0,
                    0, 0, 1, this->dT.at(i), 0, 0,
                    0, 0, 0, 1, 0, 0,
                    0, 0, 0, 0, 1, this->dT.at(i),
                    0, 0, 0, 0, 0, 1;
            this->B_d.at(i) << pow(this->dT.at(i),2)/2, 0, 0,
                    this->dT.at(i), 0, 0,
                    0, pow(this->dT.at(i),2)/2, 0,
                    0, this->dT.at(i), 0,
                    0, 0, pow(this->dT.at(i),2)/2,
                    0, 0, this->dT.at(i);
            this->B_d.at(i) = (1/this->mass)*(this->B_d.at(i));
        }
        this->A_und = Eigen::MatrixXd::Zero((this->m)*2*(this->dim), 2*(this->dim));
        Eigen::MatrixXd temp(2*(this->dim), 2*(this->dim));
        temp = A_d.at(0);
        for (int i {0}; i < this->m; ++i)
        {
            this->A_und.block(2*i*(this->dim), 0, 2*(this->dim), 2*(this->dim)) = temp;
            if (i < (this->m-1))
                temp = temp * A_d.at(i+1);
        }
        this->B_und = Eigen::MatrixXd::Zero((this->m)*2*(this->dim), (this->m)*(this->dim));
        for (int i {0}; i < this->m; ++i)
        {
            this->B_und.block(i*2*(this->dim), i*(this->dim), 2*(this->dim), (this->dim)) = B_d.at(i);
            if (i > 0)
            {
                Eigen::MatrixXd temp_B_und_block(2*(this->dim), (this->dim));
                for (int k {i-1}; k >= 0; --k)
                {
                    temp_B_und_block = B_d.at(k);
                    for (int t {k+1}; t <= i; ++t)
                    {
                        temp_B_und_block = A_d.at(t)*temp_B_und_block;
                    }
                    this->B_und.block(i*2*(this->dim), k*(this->dim), 2*(this->dim), (this->dim)) = temp_B_und_block;
                }
            }

        }
    }
    else if (this->dim == 2)
    {
        throw std::runtime_error("mpc not designed for this dim yet");
    }
    else
    {
        throw std::runtime_error("mpc not designed for this dim yet");
    }

//    std::cout << "A_und is\n" << this->A_und << std::endl;
//    std::cout << "B_und is\n" << this->B_und << std::endl;
}

void MPC_Mosek_xN::set_u_constraints(Eigen::VectorXd u_max1, Eigen::VectorXd u_min1) {
    /* maximum-minimum values for virtual forces */
    this->u_max = u_max1;
    this->u_min = u_min1;
    /* now we need to construct A1, l1, u1 */
    this->A1 = Eigen::MatrixXd::Identity(this->num_var, this->num_var);
    this->l1 = Eigen::VectorXd::Zero(this->num_var);
    this->u1 = Eigen::VectorXd::Zero(this->num_var);
    for (int k {0}; k < this->m; ++k)
    {
        for (int d {0}; d < this->dim; ++d)
        {
            this->l1((this->dim)*k + d) = this->u_min(d);
            this->u1((this->dim)*k + d) = this->u_max(d);
        }
    }
    this->create_msk_const_indices(this->dim, this->m);
    this->is_u_lim_set = true;
}

void MPC_Mosek_xN::set_u_continuity_constraint() {
    /* now we need to construct A2, l2, u2 */
    this->A2 = Eigen::MatrixXd::Identity((this->num_var - this->dim), this->num_var);
    this->l2 = Eigen::VectorXd::Zero(this->num_var - this->dim);
    this->u2 = Eigen::VectorXd::Zero(this->num_var - this->dim);
    for (int k {0}; k < (this->m - 1); ++k)
    {
        for (int d {0}; d < this->dim; ++d)
        {
            this->l2(k*(this->dim) + d) = -10*this->dT.at(k);
            this->u2(k*(this->dim) + d) = 10*this->dT.at(k);
            this->A2(k*(this->dim) + d, (k+1)*(this->dim) + d) = 1;
            this->A2(k*(this->dim) + d, (k)*(this->dim) + d) = -1;
        }
    }
    this->is_continuity_set = true;
}

void MPC_Mosek_xN::set_x0(const Eigen::VectorXd &init_pos, const Eigen::VectorXd &init_vel) {
    assert(this->dim == init_pos.size());
    assert(this->dim == init_vel.size());
    for (int d{0}; d < this->dim; ++d)
    {
        this->x0(2*d) = init_pos(d);
        this->x0(2*d+1) = init_vel(d);
    }
    this->is_cvx_const_update = false;
}

void MPC_Mosek_xN::reset_cvx_constraint() {
    this->Ac_und = Eigen::MatrixXd::Zero(m*MAX_NUM_CVX_CONS, 2*num_var);
    this->bc_und = Eigen::VectorXd::Zero(m*MAX_NUM_CVX_CONS);
    this->A3 = Eigen::MatrixXd::Zero(m*MAX_NUM_CVX_CONS, num_var);
//    this->l3 = Eigen::VectorXd::Zero(m*MAX_NUM_CVX_CONS);
    this->u3 = Eigen::VectorXd::Zero(m*MAX_NUM_CVX_CONS);
}

void MPC_Mosek_xN::update_cvx_constraints(const Eigen::MatrixXd &A_c, const Eigen::VectorXd &b_c) {
//    this->reset_cvx_constraint();
//    assert(m*MAX_NUM_CVX_CONS >= m*A_c.rows());
    this->num_cvx_const = this->m * A_c.rows();
//    assert(A_c.rows()==b_c.size());
    int s {A_c.rows()};
//    std::cout << "A_c.cols() is: " << A_c.cols() << std::endl;

    for (int k {0}; k < this->m ; ++k)
    {
        this->Ac_und.block(k*s, 2*(this->dim)*k, s, 2*(this->dim)) = A_c*I_p(this->dim);
//        std::cout << "Ac_und.block operation done!" << std::endl;
        this->bc_und.segment(k*s, s) = b_c;
    }

    this->A3.block(0,0,this->num_cvx_const,this->num_var) = this->Ac_und.block(0,0,this->num_cvx_const,2*(this->num_var)) * this->B_und;
    this->u3.head(this->num_cvx_const) = this->bc_und.head(this->num_cvx_const) - this->Ac_und.block(0,0,this->num_cvx_const,2*(this->num_var)) * this->A_und * this->x0;

    this->is_cvx_const_update = true;
}

void MPC_Mosek_xN::create_utility_matrices(const int &dim, const int &m) {
    this->I_v_mask = Eigen::MatrixXd::Zero(dim*m, 2*dim*m);
    for (int k {0}; k < m; ++k)
    {
        this->I_v_mask.block(k*dim, 2*dim*k, I_v(dim).rows(), I_v(dim).cols()) = I_v(dim);
    }

}

void MPC_Mosek_xN::update_obj_matrices(const Eigen::VectorXd &v_h, double w1, double w2, double w3) {

//    std::cout << "entering update obj matrices" << std::endl;
//    this->Q1 = this->I_v_mask;
//    std::cout << "Q1 is:\n" << this->Q1 << std::endl;

    for (int k {0}; k < this->m; ++k)
    {
        this->V_h.segment(k*(this->dim), (this->dim)) = v_h;
    }

//    std::cout << "update_obj_matrices check point 1" << std::endl;
//    std::cout << "x0 is:\n" << this->x0 << std::endl;

    this->V_h.noalias() = this->V_h - (this->I_v_mask)*(this->A_und * this->x0);
    this->Q1.noalias() = this->I_v_mask * this->B_und;
    this->Q1.noalias() = this->Q1.transpose() * this->Q1;
//    this->Q1 = this->B_und.transpose()*this->I_v_mask.transpose()*this->I_v_mask * this->B_und;
    this->c1.noalias() = -this->Q1.transpose()*(this->V_h);
//    std::cout << "update_obj_matrices check point 2" << std::endl;

    this->Q3.block((this->m - 1)*(this->dim), 0, this->dim , (this->m)*(this->dim)) =
            I_p(this->dim)*((this->B_und).block(2*(this->dim)*(this->m - 1), 0, 2*(this->dim), (this->m)*(this->dim)));
    this->c3 = (this->B_und).block(2*(this->dim)*(this->m - 1), 0, 2*(this->dim), (this->m)*(this->dim)).transpose()*I_p(this->dim).transpose()* (I_p(this->dim)*(this->A_d.at(m-1))*(this->x0) - this->xN);
    this->Q3 = 2 * (this->Q3).transpose() * (this->Q3);
    this->c3 = 2 * this->c3;

    this->Qo = w1*(this->Q1) + w2*(this->Q2) + w3*(this->Q3);

    /* need to make sure it is symmetric */
    Eigen::MatrixXd Qo_transpose(this->num_var, this->num_var);
    Qo_transpose = (this->Qo).transpose();
    this->Qo = Qo_transpose + (this->Qo);
    this->c = w1*(this->c1) + w3*(this->c3);

}

void MPC_Mosek_xN::create_msk_const_indices(const int &dim, const int &m) {
    /* objective paramters */
    /* quadratic part */
//    std::cout <<"entering create_msk_const_indices" << std::endl;
    this->numqonz = (pow(dim*m,2) + dim*m)/2;
    this->qosubi = new MSKint32t [this->numqonz];
    this->qosubj = new MSKint32t [this->numqonz];
    this->qovals = new MSKrealt [this->numqonz];
//    std::cout << "create_msk_const_indices check point 1" << std::endl;
    int counter {0};
    for (int i {0}; i < dim*m; ++i)
    {
        for (int j {0}; j <= i; ++j)
        {
            this->qosubi[counter] = i;
            this->qosubj[counter] = j;
//            this->qovals[counter] = this->Qo(i,j);
            ++counter;
        }
    }
//    std::cout << "create_msk_const_indices check point 2" << std::endl;
    counter = 0;
    /* c */
    this->num_clist = dim*m;
    this->subj_clist = new MSKint32t [dim*m];
    for (int i {0}; i < dim*m; ++i)
    {
        this->subj_clist[counter] = i;
        ++counter;
    }
//    std::cout << "create_msk_const_indices check point 3" << std::endl;
    counter = 0;
    /* boundary on variables */
    this->num_var_bounds = dim*m;
    this->var_bounds_subi = new MSKint32t [dim*m];
    this->var_bound_keys = new MSKboundkeye [dim*m];
    this->blx = new MSKrealt [dim*m];
    this->bux = new MSKrealt [dim*m];
//    std::cout << "create_msk_const_indices check point 4" << std::endl;
//    std::cout << "size of l1 and u1 are: " << this->l1.size() << " and " << this->u1.size() << std::endl;
//    std::cout << "l1 now isL\n" << this->l1 << std::endl;
//    std::cout << "u1 now isL\n" << this->u1 << std::endl;
    for (int k {0}; k < m; ++k)
    {
        for (int d {0}; d < dim; ++d)
        {
            this->var_bound_keys[k*dim + d] = MSK_BK_RA;
            this->var_bounds_subi[k*dim + d] = k*dim + d;
//            std::cout << "k*dim + d is now: " << k*dim+d << std::endl;
            this->blx[k*dim + d] = this->l1(k*dim+d);
            this->bux[k*dim + d] = this->u1(k*dim+d);
        }
    }

//    std::cout << "create_msk_const_indices check point 5" << std::endl;
    /* l<=A<=b */
    counter = 0;
    this->num_alist_1 = (m-1)*dim*2;
    this->subi_alist_1 = new MSKint32t [(m-1)*dim*2];
    this->subj_alist_1 = new MSKint32t [(m-1)*dim*2];
    this->a_vals_1 = new MSKrealt [(m-1)*dim*2];

    this->num_bound_list_1 = (m-1)*dim;
    this->bound_sub_list_1 = new MSKint32t [(m-1)*dim];
    this->bound_key_list_1 = new MSKboundkeye [(m-1)*dim];
    this->blc_list_1 = new MSKrealt [(m-1)*dim];
    this->buc_list_1 = new MSKrealt [(m-1)*dim];
    int counter2 {0};
//    std::cout << "create_msk_const_indices check point 6" << std::endl;
    for (int k {0}; k < (this->m - 1); ++k)
    {
        for (int d {0}; d < this->dim; ++d)
        {
//            this->l2(k*(this->dim) + d) = -10*this->dT.at(k);
//            this->u2(k*(this->dim) + d) = 10*this->dT.at(k);
//            this->A2(k*(this->dim) + d, (k+1)*(this->dim) + d) = 1;
//            this->A2(k*(this->dim) + d, (k)*(this->dim) + d) = -1;
            this->subi_alist_1[counter] = k*dim+d;
            this->subj_alist_1[counter] = (k+1)*(dim) + d;
            this->a_vals_1[counter] = 1;
            ++counter;
            this->subi_alist_1[counter] = k*dim+d;
            this->subj_alist_1[counter] = (k+1)*(dim) + d;
            this->a_vals_1[counter] = -1;
            ++counter;
            this->bound_sub_list_1[counter2] = counter2;
            this->bound_key_list_1[counter2] = MSK_BK_RA;
            this->blc_list_1[counter2] = -10*this->dT.at(k);
            this->buc_list_1[counter2] = 10*this->dT.at(k);
            ++counter2;
        }
    }
//    std::cout << "create_msk_const_indices last check point" << std::endl;
//    std::cout << "Exiting create_msk_const_indices" << std::endl;
}

void MPC_Mosek_xN::set_k_steps(int k) {
    this->k_steps = k;
}

int MPC_Mosek_xN::get_k_steps() const {
    return this->k_steps;
}
bool MPC_Mosek_xN::solve_mpc(MSKenv_t *existing_env) {
    /* returns false if it was not feasible */
    bool return_value = false;
    /* This function should only be called when everything is updated */
    MSKenv_t *env;
    if (existing_env) {
//        printf("env existed!\n");
        env = existing_env;
    } else {
//        printf("env not existed!\n");
        env = (MSKenv_t*) malloc(sizeof(MSKenv_t));
        check_res(MSK_makeenv(env, NULL));
    }
    MSKtask_t task = NULL;
    /* this->num_var */
    /* we need to calculate the number of constraints */
    int num_con = this->A2.rows() + this->num_cvx_const;
    check_res(MSK_maketask(*env, num_con, this->num_var, &task));
    check_res(MSK_appendcons(task, num_con));
    check_res(MSK_appendvars(task, this->num_var));

    /* Limits on variables */
    check_res(MSK_putvarboundlist(task, this->num_var_bounds, this->var_bounds_subi,
                                  this->var_bound_keys, this->blx, this->bux));
    /* Creating Objective */

    int counter {0};
    for (int i {0}; i < dim*m; ++i)
    {
        for (int j {0}; j <= i; ++j)
        {
            this->qovals[counter] = this->Qo(i,j);
            ++counter;
        }
    }
    check_res(MSK_putqobj(task, this->numqonz, this->qosubi, this->qosubj, this->qovals));
    check_res(MSK_putclist(task, this->num_clist, this->subj_clist, this->c.data()));
    check_res(MSK_putcfix(task, (this->V_h).dot(this->V_h) + (this->xN).dot(this->xN)));

    /* constraints */
    check_res(MSK_putaijlist(task, this->num_alist_1, this->subi_alist_1,
                             this->subj_alist_1, this->a_vals_1));
    check_res(MSK_putconboundlist(task, this->num_bound_list_1,
                                  this->bound_sub_list_1, this->bound_key_list_1,
                                  this->blc_list_1, this->buc_list_1));

//    MSKint32t first = 0;
//    MSKint32t last = this->A3.cols();
//    MSKint32t *ptrb[A3.cols()];
//    MSKint32t *ptre[A3.cols()];
    MSKint32t asub[this->num_cvx_const];
//    MSKrealt *aval = A3.data();
    MSKboundkeye bkc[this->num_cvx_const];
    MSKrealt blc[this->num_cvx_const];
    for (int i {0}; i < this->num_cvx_const; ++i)
    {
        asub[i] = (m-1)*dim + i;
        bkc[i] = MSK_BK_UP;
        blc[i] = -MSK_INFINITY;
    }
//    for (int j {0}; j < A3.cols(); ++j)
//    {
//        ptrb[j] = A3.col(j).data();
//        ptrb[j] = A3.col(j).data() + this->num_cvx_const + 1;
//    }
//    check_res(MSK_putacolslice(task, first, last, ptrb, ptre, asub, aval));

    MSKint32t num_aij_2 = (this->num_cvx_const)*(this->num_var);
    MSKint32t subi[num_aij_2];
    MSKint32t subj[num_aij_2];
    counter = 0;
    for (int i {0}; i < this->num_cvx_const; ++i)
    {
        for (int j {0}; j < A3.cols(); ++j)
        {
            subi[counter] = i + this->A2.rows();
            subj[counter] = j;
            ++counter;
        }
    }
    check_res(MSK_putaijlist(task, num_aij_2, subi, subj, A3.data()));
    MSKint32t num_bound = this->num_cvx_const;
    check_res(MSK_putconboundlist(task, num_bound, asub, bkc, blc, this->u3.data()));

    check_res(MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE));
    MSKrescodee trmcode;
    check_res(MSK_optimizetrm(task, &trmcode));
    MSKsolstae solsta;
    MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
//    std::cout << "solution status: " << solsta << std::endl;
    double sol_xx[this->num_var];

    switch(solsta) {
        case MSK_SOL_STA_OPTIMAL:
        case MSK_SOL_STA_NEAR_OPTIMAL:
            check_res(MSK_getxx(task, MSK_SOL_ITR, sol_xx));
            for (int i {0}; i < this->k_steps; ++i)
            {
                for (int d {0}; d < this->dim; ++d)
                {
                    this->mpc_sol_u_k(d, i) = sol_xx[i*(this->dim) + d];
                }
                if (i==0)
                {
                    this->mpc_sol_x_k.col(0) = this->A_d.at(i)*this->x0 + this->B_d.at(i)*this->mpc_sol_u_k.col(i);
                }
                else
                {
                    this->mpc_sol_x_k.col(i) = this->A_d.at(i)*this->mpc_sol_x_k.col(i-1) + this->B_d.at(i)*this->mpc_sol_u_k.col(i);
                }
            }

            return_value = true;
            break;
        case MSK_SOL_STA_DUAL_INFEAS_CER:
        case MSK_SOL_STA_PRIM_INFEAS_CER:
        case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
        case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
            printf("Primal or dual infeasibility certificate found in mpc_mosek_xN.\n");
            break;
        default:
            std::cout << "Other solution status!" << solsta << std::endl;
    }
    if (!existing_env) {
//        printf("It's here, before return\n");
        MSK_deleteenv(env);
        free(env);
    }
    else if (!env) {
        MSK_deleteenv(env);
    }
    return return_value;
}

bool MPC_Mosek_xN::do_mpc_mosek_xN(const Eigen::VectorXd &initial_position,
                                   const Eigen::VectorXd &initial_velocity,
                                   const Eigen::VectorXd &next_way_pt,
                                   const Eigen::VectorXd &v_h,
                                   const Eigen::MatrixXd &A_c,
                                   const Eigen::VectorXd &b_c) {

//    std::cout << "entering do_mpc_mosek_xN" << std::endl;
    this->set_x0(initial_position, initial_velocity);
//    std::cout << "set_x0 done" << std::endl;
    this->xN = next_way_pt;
//    std::cout << "setting xN done" << std::endl;
    if (this->is_continuity_set)
    {
//        auto begin = std::chrono::high_resolution_clock::now();
        this->set_u_constraints(Eigen::VectorXd::Constant(this->dim, 10),
                                Eigen::VectorXd::Constant(this->dim, -10));
//        auto end = std::chrono::high_resolution_clock::now();
//        auto elapsed = std::chrono::duration_cast<std::chrono::duration<float>>(end - begin);
//        std::cout << "It took " << elapsed.count() << " seconds!" << std::endl;
//        std::cout << "set_u_constraints done" << std::endl;
    }

    this->reset_cvx_constraint();
//    std::cout << "reset_cvx done" << std::endl;
    this->update_cvx_constraints(A_c, b_c);
//    std::cout << "update_cvx done" << std::endl;
//    auto begin = std::chrono::high_resolution_clock::now();
    this->update_obj_matrices(v_h);
//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::duration<float>>(end - begin);
//    std::cout << "It took " << elapsed.count() << " seconds!" << std::endl;
//    std::cout << "number of threads " << Eigen::nbThreads() << std::endl;
//    std::cout << "update_obj_matrices done" << std::endl;


//    auto begin = std::chrono::high_resolution_clock::now();
    bool solve_mpc_result {solve_mpc()};
//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::duration<float>>(end - begin);
//    std::cout << "It took " << elapsed.count() << " seconds!" << std::endl;
    return solve_mpc_result;
}