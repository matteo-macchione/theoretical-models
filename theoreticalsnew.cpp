#include <iostream>
#include <cmath>
#include <fstream>
#include "Func.h"
#include "GSLwrapper.h"
using namespace cbl;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Costanti
const double G = 6.674e-11;               // [m^3 kg^-1 s^-2]
const double M_sun = 1.989e30;            // [kg]
const double Mpc = 3.086e22;              // [m]
const double c = 3*pow(10, 8);            // [m/s]

// Densità critica
double rho_crit(double H0) {
    double H0_SI = H0 * 1e3 / Mpc;   // km/s/Mpc -> s^-1
    return 3 * H0_SI * H0_SI / (8 * M_PI * G);  // [kg/m^3]
}

// Relazione massa–concentrazione (Duffy+08, z=costante)
double c500_duffy(double M500, double z) {
    double A = 0.63, B = -0.08, C = -0.47;
    double Mpivot = 2e12; // pivot mass [M_sun/h], assumo h=1
    return A * pow(M500 / Mpivot, B) * pow(1 + z, C);
}

// Profilo densità NFW
double rho_NFW(double x, double R500_SI, double M500_SI, double c500, double g500) {
    //return rho_s / ((r / r_s) * pow(1 + r / r_s, 2));
    return (M500_SI*pow(c500, 2)*g500)/(4*M_PI*pow(R500_SI, 3)*x*pow(1+c500*x, 2));
}

// Potenziale
double potential(double x, double M500, double R500, double g500, double c500) {
    return -g500*(G*M500/R500)*(log(1+c500*x)/x);
}

// Sigma tilde (adimensionale) usando integrazione lungo x
double Sigma_tilde_x(double x_perp, double c500, double R500_SI, double M500_SI, int Ntheta = 1000)
{
    double g500 = 1.0 / (std::log(1 + c500) - c500 / (1 + c500));

    // massimo raggio di integrazione (adimensionale)
    double x_max = 100.0; 

    auto integrand = [&](double x) -> double {
        if (x <= x_perp) return 0.0; // evita radice negativa
        double rho_r = (M500_SI * std::pow(c500, 2) * g500) /
                       (4.0 * M_PI * std::pow(R500_SI, 3) * x * std::pow(1 + c500 * x, 2));
        double jac = x / std::sqrt(x*x - x_perp*x_perp); // dx/dtheta equivalente
        return rho_r * jac;
    };

    double epsabs = 1e-6;
    double epsrel = 1e-6;
    int limit_size = 5000;
    int rule = 6; // GSL_INTEG_GAUSS61

    double integral = wrapper::gsl::GSL_integrate_qag(
        integrand, x_perp, x_max, epsabs, epsrel, limit_size, rule
    );

    return 2.0 * integral; // fattore 2 per simmetria sopra e sotto il piano del cielo
}

// Sigma fisico
double sigma(double x_perp, double rho_s, double c500, double R500_SI, double M500_SI)
{
    return Sigma_tilde_x(x_perp, c500, R500_SI, M500_SI) * 1.0; // kg/m^2 già incluso in rho_r
}


//calcolo del velocity offset per un singolo cluster
double delta(double x_perp, double R500_SI, double phi0,
             double rho_s, double M500_SI, double c500)
{
    // g500 per NFW
    double g500 = 1.0 / (std::log(1 + c500) - c500 / (1 + c500));

    // integratore adattivo lungo x: da x_perp a un massimo (es. 100 R500)
    double x_max = 100.0; // adimensionale, puoi aumentare se serve
    auto integrand = [&](double x) -> double {
        if (x <= x_perp) return 0.0; // evita radice negativa
        double rho_r = (M500_SI * std::pow(c500, 2) * g500) /
                       (4.0 * M_PI * std::pow(R500_SI, 3) * x * std::pow(1 + c500 * x, 2));
        double phi_r = -g500 * (G * M500_SI / R500_SI) * (std::log(1 + c500 * x) / x);

        // jacobiano della trasformazione dx/dtheta = x / sqrt(x^2 - x_perp^2)
        double jac = x / std::sqrt(x*x - x_perp*x_perp);

        return rho_r * (phi0 - phi_r) * jac;
    };

    // Parametri GSL
    double epsabs = 1e-6;
    double epsrel = 1e-6;
    int limit_size = 5000;
    int rule = 6; // GSL_INTEG_GAUSS61

    // integrazione
    double integral = wrapper::gsl::GSL_integrate_qag(
        integrand, x_perp, x_max, epsabs, epsrel, limit_size, rule
    );

    // Surface density Σ(x⊥)
    double Sigma = sigma(x_perp, rho_s, c500, R500_SI, M500_SI);

    return (2.0 / (c * Sigma)) * integral; // [m/s]
}


int main() {
    int N = 1000;
    // Input
    double H0 = 70.0;       // km/s/Mpc
    double z = 0.5;         // redshift
    double M500 = 1e14;     // [M_sun]

    // Conversione massa in SI
    double M500_SI = M500 * M_sun;

    // Calcolo densità critica
    double rhoc = rho_crit(H0);  // [kg/m^3]

    // Calcolo R500 in metri
    double R500_SI = pow( (3 * M500_SI) / (4 * M_PI * 500 * rhoc), 1.0/3.0 );

    // Conversione in Mpc
    double R500 = R500_SI / Mpc;

    // Calcolo concentrazione con Duffy
    double c500 = c500_duffy(M500, z);

    // Parametri NFW
    double r_s = R500 / c500;   // [Mpc]
    double r_s_SI = r_s * Mpc;  // [m]
    double factor = std::log(1 + c500) - c500 / (1 + c500);
    double rho_s = M500_SI / (4 * M_PI * pow(r_s_SI, 3) * factor);  // [kg/m^3]

    // Output
    std::cout << "=== Parametri NFW da M500 ===\n";
    std::cout << "M500   = " << M500 << " M_sun\n";
    std::cout << "R500   = " << R500 << " Mpc\n";
    std::cout << "c500   = " << c500 << "\n";
    std::cout << "r_s    = " << r_s << " Mpc\n";
    std::cout << "rho_s  = " << rho_s << " kg/m^3\n\n";

    //creo una griglia
    double r[N];
    double delta_r = 3.0/N;
    for (int i=0; i<N-1; i++) {
        r[i+1]=r[i]+delta_r;
    };

    //calcolo NFW
    double NFW[N];
    double r_tilde[N]; 
    double g500 = pow(log(1+c500)-(c500/(1+c500)), -1);
    for (int i=1; i<N-1; i++) {
        r_tilde[i] = r[i]/R500;
        NFW[i] = rho_NFW(r_tilde[i], R500_SI, M500_SI, c500, g500);
    }

   
    double phi[N];
    phi[0] = -g500*c500*G*M500_SI/R500_SI;

    // calcolo il potenziale
    for (int i=1; i<N-1; i++) {
        phi[i] = potential (r_tilde[i], M500_SI, R500_SI, g500, c500);
        //std::cout << phi[0]-phi[i] << "\n";
    }

    //calcolo la surface density
    std::ofstream fout("sigma_profile_gsl.dat");
    fout << "x_perp,Sigma_kgm2\n";

    double Sigma[N];
    for (int i = 0; i < N-1; i++) {
        double x_perp = r_tilde[i];  
        double sigma_val = sigma(x_perp, rho_s, c500, R500_SI, M500_SI);
        fout << x_perp << "   " << sigma_val << "\n";
        Sigma[i] = sigma_val;
    }

    fout.close();

    //calcolo il velocity offset di un singolo cluster
    std::ofstream fout1("delta_single_cluster_gsl.dat");
    fout1 << "x_perp,Delta\n";

    double Delta[N];
    
    for (int i=1; i<N-1; i++) {
        double x_perp = r_tilde[i];
        double delta_val = delta(x_perp, R500_SI, phi[0], rho_s, M500_SI, c500);
        fout1 << x_perp << "   " << delta_val/1000 << "\n";   //il /1000 converte in Km/s
        Delta[i] = delta_val;
        Delta[0] = 0;
    }

    fout1.close();
    
    //calcolo il velocity offset per una popolazione di cluster
   
    return 0;
};