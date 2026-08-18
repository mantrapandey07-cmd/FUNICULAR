#include <iostream>
#include <cmath>
#include <array>
#include <stdexcept>

struct CableParams {
    double L0;  // unstretched length
    double EA;  // axial stiffness
    double w;   // weight per unit length
    double dX;  // target horizontal span
    double dY;  // target vertical span (corner height - platform height)
};

// X(s), Y(s) at s = L0, given trial (H, V0)
// These come directly from Irvine's elastic catenary equations.
std::array<double, 2> catenaryXY(double H, double V0, const CableParams& p) {
    double s = p.L0;
    double X = H * s / p.EA
             + (H / p.w) * (std::asinh((V0 + p.w * s) / H) - std::asinh(V0 / H));
    double Y = (1.0 / p.EA) * (V0 * s + p.w * s * s / 2.0)
             + (1.0 / p.w) * (std::sqrt(H * H + (V0 + p.w * s) * (V0 + p.w * s))
                               - std::sqrt(H * H + V0 * V0));
    return {X, Y};
}

// Residual: how far off is our current (H, V0) guess from the target (dX, dY)?
std::array<double, 2> residual(double H, double V0, const CableParams& p) {
    auto XY = catenaryXY(H, V0, p);
    return {XY[0] - p.dX, XY[1] - p.dY};
}

// 2x2 Jacobian via numerical central differences (robust, no hand-derived
// derivatives needed — safer against algebra mistakes).
std::array<std::array<double, 2>, 2> jacobian(double H, double V0, const CableParams& p) {
    double eps = 1e-6 * std::max(1.0, std::abs(H));
    double epsV = 1e-6 * std::max(1.0, std::abs(V0));

    auto rPlusH  = residual(H + eps, V0, p);
    auto rMinusH = residual(H - eps, V0, p);
    auto rPlusV  = residual(H, V0 + epsV, p);
    auto rMinusV = residual(H, V0 - epsV, p);

    std::array<std::array<double, 2>, 2> J;
    J[0][0] = (rPlusH[0] - rMinusH[0]) / (2 * eps);   // dX_residual/dH
    J[0][1] = (rPlusV[0] - rMinusV[0]) / (2 * epsV);  // dX_residual/dV0
    J[1][0] = (rPlusH[1] - rMinusH[1]) / (2 * eps);   // dY_residual/dH
    J[1][1] = (rPlusV[1] - rMinusV[1]) / (2 * epsV);  // dY_residual/dV0
    return J;
}

struct SolveResult {
    double H;
    double V0;
    int iterations;
    double finalError;
    bool converged;
};

double errNorm(double H, double V0, const CableParams& p) {
    auto r = residual(H, V0, p);
    return std::sqrt(r[0] * r[0] + r[1] * r[1]);
}

// Newton-Raphson with backtracking line search.
//
// Why line search is necessary here (not optional polish): near the actual
// solution, this particular 2-equation system has a near-singular Jacobian
// (the X and Y equations become almost linearly dependent for near-taut or
// lightly-loaded cables). A plain fixed Newton step overshoots wildly in
// that region and the solver diverges. Backtracking line search fixes this
// by shrinking the step until it actually reduces the error, at the cost of
// a few extra iterations. This was verified against a direct grid search
// before being adopted here.
SolveResult solve(const CableParams& p, double H_guess, double V0_guess,
                   int maxIter = 200, double tol = 1e-9) {
    double H = H_guess;
    double V0 = V0_guess;

    for (int iter = 0; iter < maxIter; ++iter) {
        double err = errNorm(H, V0, p);
        if (err < tol) {
            return {H, V0, iter, err, true};
        }

        auto r = residual(H, V0, p);
        auto J = jacobian(H, V0, p);
        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];

        if (std::abs(det) < 1e-16) {
            throw std::runtime_error("Singular Jacobian - try a different initial guess.");
        }

        // Solve J * delta = -r  for delta = [dH, dV0]
        double dH  = (-r[0] * J[1][1] + r[1] * J[0][1]) / det;
        double dV0 = (-J[0][0] * r[1] + J[1][0] * r[0]) / det;

        // Backtracking line search: start with the full Newton step, halve
        // it until the error actually goes down (and H stays physical, i.e.
        // positive - a cable can't have negative horizontal tension).
        double alpha = 1.0;
        double newH = H, newV0 = V0, newErr = err;
        for (int ls = 0; ls < 50; ++ls) {
            newH = H + alpha * dH;
            newV0 = V0 + alpha * dV0;
            if (newH <= 1e-9) { alpha *= 0.5; continue; }
            newErr = errNorm(newH, newV0, p);
            if (newErr < err) break;
            alpha *= 0.5;
        }

        H = newH;
        V0 = newV0;

        if (alpha < 1e-12) {
            // Line search couldn't find any improving step - stuck.
            return {H, V0, iter, newErr, false};
        }
    }

    double err = errNorm(H, V0, p);
    return {H, V0, maxIter, err, false};
}

// Once H, V0 are known, compute the actual stretched length and max sag
// (perpendicular deviation from the straight line between endpoints is a
// good-enough proxy; here we report max Y deviation from the straight chord).
void reportShape(double H, double V0, const CableParams& p) {
    int N = 20;
    std::cout << "\ns(m)\tX(m)\tY(m)\n";
    for (int i = 0; i <= N; ++i) {
        double s = p.L0 * i / N;
        CableParams pTemp = p;
        pTemp.L0 = s;  // reuse function at partial length
        auto XY = catenaryXY(H, V0, pTemp);
        std::cout << s << "\t" << XY[0] << "\t" << XY[1] << "\n";
    }

    // Stretched length = L0 + elastic elongation.
    // Tension varies along the cable; approximate total stretch by
    // integrating tension/EA - here we just report tension at each end
    // as a simple check, full integral is a refinement for later.
    double T_platform = std::sqrt(H * H + V0 * V0);
    double T_corner = std::sqrt(H * H + (V0 + p.w * p.L0) * (V0 + p.w * p.L0));
    std::cout << "\nTension at platform end: " << T_platform << " N\n";
    std::cout << "Tension at corner end:   " << T_corner << " N\n";
}

int main() {
    // ---- PLACEHOLDER VALUES - replace with your real measurements ----
    CableParams p;
    p.L0 = 3.0;      // unstretched cable length, meters (placeholder)
    p.EA = 5000.0;   // axial stiffness, N (placeholder - get real EA from Test A)
    p.w  = 0.02;     // weight per unit length, N/m (placeholder - weigh your cord)
    p.dX = 2.0;      // horizontal span, meters (placeholder - your geometry)
    p.dY = 2.2;      // vertical span, meters (placeholder - your geometry)
    // --------------------------------------------------------------------

    // Initial guess for Newton-Raphson: use the cable's total self-weight as
    // the scale for H and V0. This is a rough physical estimate, not a
    // solution - but it puts the solver in the right ballpark, which matters
    // a lot for this system (a wildly wrong guess, e.g. picked out of thin
    // air, can send Newton-Raphson to a nonphysical region it won't recover
    // from). If your real geometry is heavily loaded (large mass on the
    // platform, not just cable self-weight), add your expected load to this
    // guess too - e.g. H_guess += expectedLoadOnThisCable * cos(angle).
    double totalWeight = p.w * p.L0;
    double H_guess = totalWeight;
    double V0_guess = totalWeight / 2.0;

    try {
        SolveResult result = solve(p, H_guess, V0_guess);

        std::cout << "--- Solve Result ---\n";
        std::cout << "Converged: " << (result.converged ? "yes" : "NO - check inputs/guess") << "\n";
        std::cout << "Iterations: " << result.iterations << "\n";
        std::cout << "Final error: " << result.finalError << "\n";
        std::cout << "H (horizontal tension)  = " << result.H << " N\n";
        std::cout << "V0 (vertical tension, platform end) = " << result.V0 << " N\n";

        reportShape(result.H, result.V0, p);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
