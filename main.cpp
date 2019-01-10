#include <vector>
#include <iostream>
#include <cmath>
#include "writer.hpp"

double f(double t, double u) {
    return std::exp(-2 * t) - 2 * u;
}


/// Uses the SSP RK3 method to compute u from time 0 to time T
/// for the ODE $u'=e^{-2t}-2u$
///
/// @param[out] u at the end of the call (solution at all time steps)
/// @param[out] time contains the time levels
/// @param[in] u0 the initial data
/// @param[in] dt the step size
/// @param[in] T the final time up to which to compute the solution.
///

//----------------SSPRK3Begin----------------
void SSPRK3(std::vector<double> &u, std::vector<double> &time,
            const double &u0, double dt, double T) {
    const unsigned int nsteps = std::round(T / dt);
    u.resize(nsteps + 1);
    time.resize(nsteps + 1);
    u[0] = u0;
    time[0] = 0;

    for (int i = 1; i < nsteps + 1; ++i) {
        time[i] = i * dt;
    }

    for (int j = 0; j < nsteps; ++j) {
        double k1 = f(time[j] + dt, u[j]);
        double k2 = f(time[j] + dt, u[j] + dt * k1);
        double k3 = f(time[j] + dt * 0.5, u[j] + 0.25 * dt * (k1 + k2));
        double a = 1.0 / 6;
        double b = 2.0 / 3.0;
        u[j + 1] = u[j] + dt * (a * (k1 + k2) + b * k3);
    }
    // Write your SSPRK3 code here


}
//----------------SSPRK3End----------------

int main(int argc, char **argv) {

    double T = 10.0;
    std::vector<double> dt(8);
    std::vector<double> error(8);
    const double u0 = 0.;
    for (int i = 0; i < 8; i++) {
        dt[i] = std::pow(0.5, i);
        std::vector<double> time;
        std::vector<double> u;
        SSPRK3(u, time, u0, dt[i], T);
        double uex = T * std::exp(-2. * T);
        error[i] = std::abs(u.back() - uex);
    }

    writeToFile("dt.txt", dt);
    writeToFile("error.txt", error);
    return 0;

}
