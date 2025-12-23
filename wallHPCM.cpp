#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <array>

#include "steel.h"

// =======================================================================
//                        [VARIOUS ALGORITHMS]
// =======================================================================

/**
 * @brief Solves a tridiagonal system of linear equations A*x = d (TDMA).
 */
std::vector<double> solveTridiagonal(const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d) {

    const int n = static_cast<int>(b.size());

    std::vector<double> c_star(n, 0.0);
    std::vector<double> d_star(n, 0.0);
    std::vector<double> x(n, 0.0);

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        const double m = b[i] - a[i] * c_star[i - 1];
        c_star[i] = c[i] / m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m;
    }

    x[n - 1] = d_star[n - 1];

    for (int i = n - 2; i >= 0; --i)
        x[i] = d_star[i] - c_star[i] * x[i + 1];

    return x;
}

int main() {

    // ===================================================================
    //                    NUMERICAL / PHYSICAL SETUP
    // ===================================================================

    constexpr int N = 20;
    constexpr double L = 1.0;
    constexpr double dz = L / N;     // [m]
    constexpr double dt = 1.0e-2;     // [s]
    constexpr int time_iter = 10000;

    std::vector<double> T_w_bulk(N,800.0);
    std::vector<double> T_w_bulk_old(N, 800.0);

    constexpr double Q_tot = 1000.0;      // Total heat input [W]
    constexpr double z_evap_start = 0.0; // [m]
    constexpr double z_evap_end = 0.3; // [m]
    constexpr double z_cond_start = 0.7; // [m]
    constexpr double z_cond_end = 1.0; // [m]
    constexpr double A_wall = 1.0e-4;     // Wall cross-section area [m2]

    // Cell volume (1D slab)
    const double Vcell = A_wall * dz;

    // Evaporator length
    const double L_evap = z_evap_end - z_evap_start;

    // Uniform volumetric heat source in evaporator
    const double q_vol = Q_tot / (A_wall * L_evap); // [W/m3]

    // Initialize heat source
    std::vector<double> Q_ow(N, 0.0);

    for (int i = 0; i < N; ++i) {
        const double z = i * dz;
        if (z >= z_evap_start && z <= z_evap_end) {
            Q_ow[i] = q_vol;
        }

        if (z >= z_cond_start && z <= z_cond_end) {
            Q_ow[i] = -q_vol;
        }
    }

    std::ofstream file("T_wall.dat");

    for (int n = 0; n < time_iter; ++n) {

        const double time = (n + 1) * dt;

        // store previous time level
        T_w_bulk_old = T_w_bulk;

        // ===================================================================
        //                  TRIDIAGONAL COEFFICIENTS
        // ===================================================================

        std::vector<double> aTW(N, 0.0);
        std::vector<double> bTW(N, 0.0);
        std::vector<double> cTW(N, 0.0);
        std::vector<double> dTW(N, 0.0);

        // BC: zero gradient (adiabatic)
        aTW[0] = 0.0;
        bTW[0] = 1.0;
        cTW[0] = -1.0;
        dTW[0] = 0.0;

        aTW[N - 1] = -1.0;
        bTW[N - 1] = 1.0;
        cTW[N - 1] = 0.0;
        dTW[N - 1] = 0.0;

        // ===================================================================
        //                      ASSEMBLY LOOP
        // ===================================================================

        for (int i = 1; i < N - 1; ++i) {

            const double cp = steel::cp(T_w_bulk[i]);
            const double rho = steel::rho(T_w_bulk[i]);

            const double kL = steel::k(T_w_bulk[i - 1]);
            const double kR = steel::k(T_w_bulk[i + 1]);

            aTW[i] = -kL / (rho * cp * dz * dz);
            bTW[i] = 1.0 / dt + (kL + kR) / (rho * cp * dz * dz);
            cTW[i] = -kR / (rho * cp * dz * dz);

            dTW[i] =
                T_w_bulk_old[i] / dt
                + Q_ow[i] / (rho * cp);
        }

        // ===================================================================
        //                          SOLVE
        // ===================================================================

        T_w_bulk = solveTridiagonal(aTW, bTW, cTW, dTW);

        // ===================================================================
        //                          OUTPUT
        // ===================================================================


        for (int i = 0; i < N; ++i)
            file << T_w_bulk[i] << " ";

        file << "\n";
    }

    file.close();

    return 0;
}
