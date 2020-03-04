#ifndef MPC_H
#define MPC_H

#include <vector>
#include <cppad/cppad.hpp>
#include "Eigen-3.3/Eigen/Core"

typedef CPPAD_TESTVECTOR(double) Dvector;

class MPC {
public:

    double steer;
    double throttle;

    Dvector vars; // where all the state and actuation variables will be stored
    Dvector vars_lowerbound; //lower limit for each corresponding variable in vars
    Dvector vars_upperbound; //upper limit for each corresponding variable in vars
    Dvector constraints_lowerbound; // value constraint for each corresponding constraint expression
    Dvector constraints_upperbound; // value constraint for each corresponding constraint expression

    std::vector<double> future_xs;
    std::vector<double> future_ys;

    MPC();

    virtual ~MPC();

    // Solve the model given an initial state and polynomial coefficients.
    // Return the first actuations.
    void Solve(const Eigen::VectorXd &state,
                              const Eigen::VectorXd &coeffs);
};

#endif  // MPC_H
