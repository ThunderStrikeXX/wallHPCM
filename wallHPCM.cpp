#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <array>

#include "tdma.h"
#include "steel.h"

int main() {

    // ===================================================================
    //                    NUMERICAL / PHYSICAL SETUP
    // ===================================================================

    constexpr int N = 100;                      // Number of wall cells
    constexpr double L = 1.0;                   // Wall length [m]
    constexpr double dz = L / N;                // Cell size [m]
    constexpr double dt = 1e-3;                 // Time step [s]
    constexpr int time_iter = 1000;             // Number of time iterations

    std::vector<double> T_w(N, 300.0);          // Initial wall temperature [K]
    std::vector<double> T_w_old;                // Old temperature [K]
    std::vector<double> Q(N);                   // Heat pipe power [W]

    // Tridiagonal matrix coefficients
    std::vector<double> aTW(N, 0.0);
    std::vector<double> bTW(N, 0.0);
    std::vector<double> cTW(N, 0.0);
    std::vector<double> dTW(N, 0.0);

    constexpr double Q_tot = 1e3;                           // Total heat input [W]
    constexpr double z_evap_start = 0.0;                    // Evaporator start [m]
    constexpr double z_evap_end = 0.3;                      // Evaporator end [m]
    constexpr double z_cond_start = 0.7;                    // Condenser start [m]
    constexpr double z_cond_end = 1.0;                      // Condenser end [m]
    constexpr double A_wall = 1.0e-4;                       // Wall cross-section area [m2]
    constexpr double T_amb = 300.0;                         // Ambient temperature [K]
    constexpr double L_evap = z_evap_end - z_evap_start;    // Evaporator length [m]
    constexpr double q_vol = Q_tot / (A_wall * L_evap);     // Heat volumetric source [W/m3]

    for (int i = 0; i < N; ++i) {

        const double z = i * dz;
        if (z >= z_evap_start && z <= z_evap_end) Q[i] = q_vol;
        if (z >= z_cond_start && z <= z_cond_end) Q[i] = -q_vol;
    }

    // Output file
    std::ofstream file("wall_solution.dat");

    // Time loop
    for (int n = 0; n < time_iter; ++n) {

        const double time = (n + 1) * dt;

        // store previous time level
        T_w_old = T_w;

        // Assembly loop
        for (int i = 1; i < N - 1; ++i) {

            const double cp = steel::cp(T_w[i]);
            const double rho = steel::rho(T_w[i]);

            const double k_l = 0.5 * (steel::k(T_w[i - 1]) + steel::k(T_w[i]));
            const double k_r = 0.5 * (steel::k(T_w[i]) + steel::k(T_w[i + 1]));

            aTW[i] = -k_l / (rho * cp * dz * dz);
            bTW[i] = 1.0 / dt + (k_l + k_r) / (rho * cp * dz * dz);
            cTW[i] = -k_r / (rho * cp * dz * dz);

            dTW[i] =
                T_w_old[i] / dt
                + Q[i] / (rho * cp);
        }

        // BC: zero gradient (adiabatic)
        aTW[0] = 0.0;
        bTW[0] = 1.0;
        cTW[0] = -1.0;
        dTW[0] = 0.0;

        aTW[N - 1] = -1.0;
        bTW[N - 1] = 1.0;
        cTW[N - 1] = 0.0;
        dTW[N - 1] = 0.0;

        // Solve
        T_w = tdma::solve(aTW, bTW, cTW, dTW);

        // Output
        for (int i = 0; i < N; ++i)
            file << T_w[i] << " ";

        file << "\n";
    }

	file.flush();
    file.close();

    return 0;
}