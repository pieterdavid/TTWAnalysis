#pragma once

#include <vector>
#include <utility>

#include <Math/Vector4D.h>

#define SQ(x) (x*x)
#define CB(x) (x*x*x)
#define QU(x) (x*x*x*x)

inline double cosXpm2PI3(const double x, const double pm){
    return -0.5 * (std::cos(x) + pm * std::sin(x) * std::sqrt(3.));
}

bool solveQuadratic(const double a, const double b, const double c, std::vector<double>& roots);
bool solveCubic(const double a, const double b, const double c, const double d, std::vector<double>& roots);
bool solveQuartic(const double a, const double b, const double c, const double d, const double e, std::vector<double>& roots);
bool solve2Quads(const double a20, const double a02, const double a11, const double a10, const double a01, const double a00, const double b20, const double b02, const double b11, const double b10, const double b01, const double b00, std::vector<double>& E1, std::vector<double>& E2);
bool solve2QuadsDeg(const double a11, const double a10, const double a01, const double a00, const double b11, const double b10, const double b01, const double b00, std::vector<double>& E1, std::vector<double>& E2);
bool solve2Linear(const double a10, const double a01, const double a00, const double b10, const double b01, const double b00, std::vector<double>& E1, std::vector<double>& E2);

class NeutrinosSolver {
    public:
        using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;

        NeutrinosSolver(float top_mass, float w_mass):
            t_mass(top_mass), w_mass(w_mass) {
            // Empty
        }

        std::vector<std::pair<LorentzVector, LorentzVector>> getNeutrinos(const LorentzVector& lepton1_p4, 
                const LorentzVector& lepton2_p4, 
                const LorentzVector& bjet1_p4, 
                const LorentzVector& bjet2_p4,
                const LorentzVector& met);

    private:
        float t_mass = 172.5;
        float w_mass = 80.4;
};
