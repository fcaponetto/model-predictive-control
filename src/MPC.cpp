#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;
using Eigen::VectorXd;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
//   simulator around in a circle with a constant steering angle and velocity on
//   a flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
//   presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

const int N = 10; // how many states we "lookahead" in the future
const double dt = 0.1; // how much time we expect environment changes

const double VELOCITY_MAX = 100.0; // this is what we ideally want our speed to always be

const int NUMBER_OF_STATES = 6; // px, py, psi, v, cte, epsi
const int NUMBER_OF_ACTUATIONS = 2; // steering angle, acceleration


// Mark when each variable starts for convenience
// since state and actuator variables are stored in one vector
// in the format 'x...(N)x y...(N)y...'
const int ID_FIRST_px = 0;
const int ID_FIRST_py = ID_FIRST_px + N;
const int ID_FIRST_psi = ID_FIRST_py + N;
const int ID_FIRST_v = ID_FIRST_psi + N;
const int ID_FIRST_cte = ID_FIRST_v + N;
const int ID_FIRST_epsi = ID_FIRST_cte + N;
const int ID_FIRST_delta = ID_FIRST_epsi + N;
const int ID_FIRST_a = ID_FIRST_delta + N - 1;

// weights for cost computations
const double W_cte = 1500.0;
const double W_epsi = 1500.0;
const double W_v = 1.0;
const double W_delta = 10.0;
const double W_a = 10.0;
const double W_ddelta = 150.0; // weight cost for high difference between consecutive steering actuations
const double W_da = 15.0; // weight cost for high difference between consecutive acceleration actuations

class FG_eval
{
public:
    // Fitted polynomial coefficients
    VectorXd coeffs;

    FG_eval(VectorXd coeffs) { this->coeffs = coeffs; }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void operator()(ADvector &fg, const ADvector &vars)
    {
        /// fg a vector containing the cost and all constraints
        /// var is a vector containing all states and actuations for N "lookahead" states and actuations.

        ///*********************************************************
        ///* COST DEFINED HERE
        ///*********************************************************

        fg[0] = 0.0;

        for (int i = 0; i < N; ++i) {

            const auto cte = vars[ID_FIRST_cte + i];
            const auto epsi = vars[ID_FIRST_epsi + i];
            const auto v = vars[ID_FIRST_v + i] - VELOCITY_MAX;

            fg[0] += (W_cte * cte * cte + W_epsi * epsi * epsi + W_v * v * v);
        }

        for (int i = 0; i < N - 1; ++i) {

            const auto delta = vars[ID_FIRST_delta + i];
            const auto a = vars[ID_FIRST_a + i];

            fg[0] += (W_delta * delta * delta + W_a * a * a);
        }

        for (int i = 0; i < N - 2; ++i) {

            const auto ddelta = vars[ID_FIRST_delta + i + 1] - vars[ID_FIRST_delta + i];
            const auto da = vars[ID_FIRST_a + i + 1] - vars[ID_FIRST_a + i];

            fg[0] += (W_ddelta * ddelta * ddelta + W_da * da * da);
        }

        ///*********************************************************
        ///* CONSTRAINTS DEFINED HERE
        ///*********************************************************

        // given state does not vary
        fg[ID_FIRST_px + 1] = vars[ID_FIRST_px];
        fg[ID_FIRST_py + 1] = vars[ID_FIRST_py];
        fg[ID_FIRST_psi + 1] = vars[ID_FIRST_psi];
        fg[ID_FIRST_v + 1] = vars[ID_FIRST_v];
        fg[ID_FIRST_cte + 1] = vars[ID_FIRST_cte];
        fg[ID_FIRST_epsi + 1] = vars[ID_FIRST_epsi];

        // constraints based on our kinematic model
        for (int i = 0; i < N - 1; ++i) {

            // where the current state variables of interest are stored
            // stored for readability
            const int ID_CURRENT_px = ID_FIRST_px + i;
            const int ID_CURRENT_py = ID_FIRST_py + i;
            const int ID_CURRENT_psi = ID_FIRST_psi + i;
            const int ID_CURRENT_v = ID_FIRST_v + i;
            const int ID_CURRENT_cte = ID_FIRST_cte + i;
            const int ID_CURRENT_epsi = ID_FIRST_epsi + i;
            const int ID_CURRENT_delta = ID_FIRST_delta + i;
            const int ID_CURRENT_a = ID_FIRST_a + i;

            //current state and actuations
            const auto px0 = vars[ID_CURRENT_px];
            const auto py0 = vars[ID_CURRENT_py];
            const auto psi0 = vars[ID_CURRENT_psi];
            const auto v0 = vars[ID_CURRENT_v];
            const auto cte0 = vars[ID_CURRENT_cte];
            const auto epsi0 = vars[ID_CURRENT_epsi];
            const auto delta0 = vars[ID_CURRENT_delta];
            const auto a0 = vars[ID_CURRENT_a];

            // next state
            const auto px1 = vars[ID_CURRENT_px + 1];
            const auto py1 = vars[ID_CURRENT_py + 1];
            const auto psi1 = vars[ID_CURRENT_psi + 1];
            const auto v1 = vars[ID_CURRENT_v + 1];
            const auto cte1 = vars[ID_CURRENT_cte + 1];
            const auto epsi1 = vars[ID_CURRENT_epsi + 1];

            // desired py and psi
            const auto py_desired = coeffs[3] * px0 * px0 * px0 + coeffs[2] * px0 * px0 + coeffs[1] * px0 + coeffs[0];
            const auto psi_desired = CppAD::atan(3.0 * coeffs[3] * px0 * px0 + 2.0 * coeffs[2] * px0 + coeffs[1]);

            // relationship of current state + actuations and next state
            // based on our kinematic model
            const auto px1_f = px0 + v0 * CppAD::cos(psi0) * dt;
            const auto py1_f = py0 + v0 * CppAD::sin(psi0) * dt;
            const auto psi1_f = psi0 + v0 * (-delta0) / Lf * dt;
            const auto v1_f = v0 + a0 * dt;
            const auto cte1_f = py_desired - py0 + v0 * CppAD::sin(epsi0) * dt;
            const auto epsi1_f = psi0 - psi_desired + v0 * (-delta0) / Lf * dt;

            // store the constraint expression of two consecutive states
            fg[ID_CURRENT_px + 2] = px1 - px1_f;
            fg[ID_CURRENT_py + 2] = py1 - py1_f;
            fg[ID_CURRENT_psi + 2] = psi1 - psi1_f;
            fg[ID_CURRENT_v + 2] = v1 - v1_f;
            fg[ID_CURRENT_cte + 2] = cte1 - cte1_f;
            fg[ID_CURRENT_epsi + 2] = epsi1 - epsi1_f;
        }
    }
};

//
// MPC class definition implementation.
//
MPC::MPC()
{
    /**
    * TODO: Set the number of model variables (includes both states and inputs).
    * For example: If the state is a 4 element vector, the actuators is a 2
    *   element vector and there are 10 timesteps. The number of variables is:
    *   4 * 10 + 2 * 9
    */
    // State: [x,y,ψ,v,cte,eψ]
    // Actuators: [δ,a]

    size_t n_vars = NUMBER_OF_STATES * N + NUMBER_OF_ACTUATIONS * (N-1);
    size_t n_constraints = NUMBER_OF_STATES * N;

    //**************************************************************
    //* SET INITIAL VALUES OF VARIABLES
    //**************************************************************
    vars.resize(n_vars);

    // all states except the ID_FIRST are set to zero
    // the aforementioned states will be initialized when solve() is called

    for (int i = 0; i < n_vars; ++i) {
        vars[i] = 0.0;
    }

    //**************************************************************
    //* SET UPPER AND LOWER LIMITS OF VARIABLES
    //**************************************************************

    vars_lowerbound.resize(n_vars);
    vars_upperbound.resize(n_vars);

    // all other values large values the computer can handle
    for (int i = 0; i < ID_FIRST_delta; ++i) {
        vars_lowerbound[i] = -1.0e10;
        vars_upperbound[i] = 1.0e10;
    }

    // all actuation inputs (steering, acceleration) should have values between [-1, 1]
    for (int i = ID_FIRST_delta; i < ID_FIRST_a; ++i) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }

    for (int i = ID_FIRST_a; i < n_vars; ++i) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }

    //**************************************************************
    //* SET UPPER AND LOWER LIMITS OF CONSTRAINTS
    //**************************************************************
    constraints_lowerbound.resize(n_constraints);
    constraints_upperbound.resize(n_constraints);

    // the first constraint for each state veriable
    // refer to the initial state conditions
    // this will be initialized when solve() is called
    // the succeeding constraints refer to the relationship
    // between succeeding states based on our kinematic model of the system

    for (int i = 0; i < n_constraints; ++i) {
        constraints_lowerbound[i] = 0.0;
        constraints_upperbound[i] = 0.0;
    }
}

MPC::~MPC()
{

}

void MPC::Solve(const VectorXd &state, const VectorXd &coeffs) {
    typedef CPPAD_TESTVECTOR(double) Dvector;

    const double px = state[0];
    const double py = state[1];
    const double psi = state[2];
    const double v = state[3];
    const double cte = state[4];
    const double epsi = state[5];

    vars[ID_FIRST_px] = px;
    vars[ID_FIRST_py] = py;
    vars[ID_FIRST_psi] = psi;
    vars[ID_FIRST_v] = v;
    vars[ID_FIRST_cte] = cte;
    vars[ID_FIRST_epsi] = epsi;

    constraints_lowerbound[ID_FIRST_px] = px;
    constraints_lowerbound[ID_FIRST_py] = py;
    constraints_lowerbound[ID_FIRST_psi] = psi;
    constraints_lowerbound[ID_FIRST_v] = v;
    constraints_lowerbound[ID_FIRST_cte] = cte;
    constraints_lowerbound[ID_FIRST_epsi] = epsi;

    constraints_upperbound[ID_FIRST_px] = px;
    constraints_upperbound[ID_FIRST_py] = py;
    constraints_upperbound[ID_FIRST_psi] = psi;
    constraints_upperbound[ID_FIRST_v] = v;
    constraints_upperbound[ID_FIRST_cte] = cte;
    constraints_upperbound[ID_FIRST_epsi] = epsi;

    //**************************************************************
    //* SOLVE
    //**************************************************************

    // object that computes objective and constraints
    FG_eval fg_eval(coeffs);

    // NOTE: You don't have to worry about these options
    // options for IPOPT solver
    std::string options;
    // Uncomment this if you'd like more print information
    options += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    //   of sparse routines, this makes the computation MUCH FASTER. If you can
    //   uncomment 1 of these and see if it makes a difference or not but if you
    //   uncomment both the computation time should go up in orders of magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    // Change this as you see fit.
    options += "Numeric max_cpu_time          0.5\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
            options,
            vars,
            vars_lowerbound,
            vars_upperbound,
            constraints_lowerbound,
            constraints_upperbound,
            fg_eval,
            solution);

    // comment out the lines below to debug!
    bool ok = true;
    auto cost = solution.obj_value;

    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
    if (ok) {
      std::cout << "OK! Cost:" << cost << std::endl;
    } else {
      std::cout << "SOMETHING IS WRONG!" << cost << std::endl;
    }


    /**
     * TODO: Return the first actuator values. The variables can be accessed with
     *   `solution.x[i]`.
     *
     * {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
     *   creates a 2 element double vector.
     */

    steer = solution.x[ID_FIRST_delta];
    throttle = solution.x[ID_FIRST_a];

    future_xs = {};
    future_ys = {};

    for (int i = 0; i < N; ++i) {

        const double px = solution.x[ID_FIRST_px + i];
        const double py = solution.x[ID_FIRST_py + i];

        future_xs.emplace_back(px);
        future_ys.emplace_back(py);
    }
}