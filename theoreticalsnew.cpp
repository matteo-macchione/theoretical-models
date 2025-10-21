#include <iostream>
#include <cmath>
#include <fstream>
#include <memory>
#include <omp.h>
#include "Func.h"
#include "GSLwrapper.h"
#include "LCDM.h"
#include "MassFunction.h"
using namespace cbl;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Costanti
const double G = 6.674e-11;               // [m^3 kg^-1 s^-2]
const double M_sun = 1.989e30;            // [kg]
const double Mpc = 3.086e22;              // [m]
const double c = 3*pow(10, 8);            // [m/s]
const double H0 = 70.0;       // km/s/Mpc

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
double potential(double x, double M500, double R500, double g500, double c500, double k) {
    return -g500*(k*G*M500/R500)*(log(1+c500*x)/x);
}

// Sigma tilde (adimensionale) usando integrazione lungo x
double Sigma_tilde_x(double x_perp, double c500, double R500_SI, double M500_SI)
{
    double g500 = 1.0 / (std::log(1 + c500) - c500 / (1 + c500));

    // massimo raggio di integrazione (adimensionale)
    //double x_max = 100.0; 

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
    //int rule = 6; // GSL_INTEG_GAUSS61

    double integral = wrapper::gsl::GSL_integrate_qagiu(
        integrand, x_perp, epsabs, epsrel, limit_size
    );

    return 2.0 * integral * R500_SI; // fattore 2 per simmetria sopra e sotto il piano del cielo
}

// Sigma fisico
double sigma(double x_perp, double c500, double R500_SI, double M500_SI)
{
    return Sigma_tilde_x(x_perp, c500, R500_SI, M500_SI) * 1.0; // kg/m^2 già incluso in rho_r
}

//calcolo del velocity offset per un singolo cluster
double delta(double x_perp, double R500_SI, double phi0,
             double M500_SI, double c500, int model)
{
    // g500 per NFW
    double g500 = 1.0 / (std::log(1 + c500) - c500 / (1 + c500));

    
    auto integrand = [=](double x) -> double {
        if (x <= x_perp) return 0.0; // evita radice negativa
        double rho_r = rho_NFW(x, R500_SI, M500_SI, c500, g500);
        double phi_r;
        if (model == 0) {
           phi_r = potential(x, M500_SI, R500_SI, g500, c500, 1.0); // GR
        }
        if (model == 1) { // f(R)
           phi_r = potential(x, M500_SI, R500_SI, g500, c500, 4.0/3.0); // f(R)
        }
        if (model == 2) { // DGP
            // densità NFW in unità di rho_s, quindi adimensionale
            auto integrand_mass = [&](double xp) -> double {
            return rho_NFW(xp, R500_SI, M500_SI, c500, g500) / rho_NFW(1.0, R500_SI, M500_SI, c500, g500) * xp * xp;
            };

            double epsabs = 1e-6, epsrel = 1e-6;
            int limit_size = 5000, rule = 6;

            double integral_mass = wrapper::gsl::GSL_integrate_qag(integrand_mass, 0.0, x, epsabs, epsrel, limit_size, rule);

            // massa relativa (adimensionale)
            double M_in_x_over_M500 = integral_mass; 
            if (x > 10.0) x = 10.0; // evita problemi numerici per x molto grandi
            wrapper::gsl::GSL_integrate_qag(integrand_mass, 0.0, 1.0, epsabs, epsrel, limit_size, rule);

            // DGP factor
            double y = pow((0.32 * 500), -1.0/3.0) * x * pow(M_in_x_over_M500, -1.0/3.0);
            double g_y = pow(y, 3) * (sqrt(1.0 + pow(y, -3)) - 1.0);
            double k_model = 1.0 + (2.0 / (3.0 * (1.15))) * g_y;
            phi_r = potential(x, M500_SI, R500_SI, g500, c500, k_model);
        }

        // jacobiano della trasformazione dx/dtheta = x / sqrt(x^2 - x_perp^2)
        double jac = x / std::sqrt(x*x - x_perp*x_perp);

        return rho_r*(phi0 - phi_r)*jac;
        
    };

    // Parametri GSL
    double epsabs = 1e-6;
    double epsrel = 1e-6;
    int limit_size = 5000;
    //int rule = 6; // GSL_INTEG_GAUSS61

    // integrazione
    double integral = wrapper::gsl::GSL_integrate_qagiu(integrand, x_perp, epsabs, epsrel, limit_size);

    // Surface density Σ(x⊥)
    double Sigma = sigma(x_perp, c500, R500_SI, M500_SI);

    return (2.0 * R500_SI/ (c * Sigma)) * integral; // [m/s]
    }

// ---------- helper: dati cluster da M500 (M_sun) ----------
struct ClusterParams {
    double M500;        // M_sun
    double M500_SI;     // kg
    double R500_SI;     // m
    double R500;        // Mpc
    double c500;
    double r_s_SI;      // m
    double rho_s;       // kg/m^3
    double phi0;        // m^2/s^2
};

ClusterParams params_from_M(double M500, double H0, double z) {
    ClusterParams p;
    p.M500 = M500;
    p.M500_SI = M500 * M_sun;
    double rhoc = rho_crit(H0); // [kg/m^3]
    p.R500_SI = pow( (3.0 * p.M500_SI) / (4.0 * M_PI * 500.0 * rhoc), 1.0/3.0 );
    p.R500 = p.R500_SI / Mpc;
    p.c500 = c500_duffy(M500, z); // input M500 in M_sun
    p.r_s_SI = (p.R500 / p.c500) * Mpc;
    double factor = std::log(1.0 + p.c500) - p.c500 / (1.0 + p.c500);
    p.rho_s = p.M500_SI / (4.0 * M_PI * pow(p.r_s_SI, 3) * factor);
    double g500 = 1.0 / factor;
    p.phi0 = - g500 * (G * p.M500_SI / p.R500_SI) * (std::log(1.0 + p.c500 * 1e-8) / 1e-8);
    // use asymptotic limit phi(0) ~ -g500 * (G M500 / R500) * c500 / f_c,
    // but the line above avoided division by 0 by using small x; better explicit:
    p.phi0 = - g500 * (G * p.M500_SI / p.R500_SI) * ( (std::log(1.0 + p.c500) / (1e-8)) ); // temporary
    // more robust explicit phi0 (limit x->0): log(1+c*x)/x -> c for x->0
    p.phi0 = - g500 * (G * p.M500_SI / p.R500_SI) * p.c500; // correct limit
    return p;
}

// ---------- wrapper che calcola Sigma fisico e Delta fisico per un M500 e x_perp ----------
double Sigma_phys_for_M(double x_perp, const ClusterParams &p) {
    // usa la funzione sigma_x che avevamo definito (o sigma) adattata:
    // se hai sigma_x(x_perp, rho_s, c500, R500_SI, M500_SI)
    return sigma(x_perp, p.c500, p.R500_SI, p.M500_SI); // [kg/m^2]
}

double Delta_phys_for_M(double x_perp, const ClusterParams &p, int model) {
    // usa la funzione delta(x_perp, R500_SI, phi0, rho_s, M500_SI, c500)
    double phi0;
    if (model == 1) {
        phi0 = (4.0/3.0)*p.phi0; // f(R)
    }
    else {
        phi0 = p.phi0; // DGP
    }
    
    return delta(x_perp, p.R500_SI, phi0, p.M500_SI, p.c500, model); // [m/s]
}
auto cosmo = std::make_shared<cbl::cosmology::LCDM>("Planck18");
        cbl::cosmology::MassFunction mf(cosmo);
        const double h = H0/100.0; 
// ---------- integrale su massa in variabile y = log(M) ----------
double delta_population_avg(double x_phys, double Mmin, double Mmax,
                            double H0, double z, int model) {
    // integrale su y = ln M: M = exp(y), dM = exp(y) dy
    //double y0 = std::log(Mmin);
    //double y1 = std::log(Mmax);
    
    auto integrand_num = [&](double M) -> double {
        //double M = std::exp(y);
        ClusterParams p = params_from_M(M, H0, z);
        double x_perp = x_phys/p.R500; // adimensionale
        double Sigma = Sigma_phys_for_M(x_perp, p); // kg/m^2
        double Delta = Delta_phys_for_M(x_perp, p, model); // m/s
        
        double weight = mf.mass_function(M/h, z, "Tinker", "EisensteinHu", false, "none", 500.0); 
        //double dM_dy = M; // dM = M dy
        return Sigma * Delta * weight; // units: (kg/m^2)*(m/s)*M_sun^{-alpha} * M_sun
    };

    auto integrand_den = [&](double M) -> double {
        //double M = std::exp(y);
        ClusterParams p = params_from_M(M, H0, z);
        double x_perp = x_phys/p.R500; // adimensionale
        double Sigma = Sigma_phys_for_M(x_perp, p); // kg/m^2
        double weight = mf.mass_function(M/h, z, "Tinker", "EisensteinHu", false, "none", 500.0);
        //double dM_dy = M;
        return Sigma * weight; // units: kg/m^2 * M_sun^{1-alpha}
    };

    // wrapper GSL integration parameters (tolleranze moderate)
    double epsabs = 1e-4;
    double epsrel = 1e-4;
    int limit_size = 5000;
    int rule = 6;

    double num = wrapper::gsl::GSL_integrate_qag(integrand_num, Mmin, Mmax, epsabs, epsrel, limit_size, rule);
    double den = wrapper::gsl::GSL_integrate_qag(integrand_den, Mmin, Mmax, epsabs, epsrel, limit_size, rule);

    if (den == 0.0) return 0.0;
    return num / den; // result in [m/s]
}

// ==========================================================
//  Transverse Doppler effect Δ_TD(x_perp)
// ==========================================================
// x_perp = r_perp / R500 (adimensionale)
// restituisce Δ_TD in [m/s]
double delta_TD(double x_perp, double R500_SI,
                double M500_SI, double c500, int model)
{
    // g(c500) = [ln(1+c500) - c500/(1+c500)]^{-1}
    double g500 = 1.0 / (std::log(1 + c500) - c500 / (1 + c500));

    // Funzione dphi/dr del potenziale NFW calcolata con GSL
    auto dphi_dr = [&](double x) -> double {
        auto phi_func = [&](double x_val) -> double {
            if (x_val <= x_perp) return 0.0; // evita radice negativa
            double phi_r;
            if (model == 0) {
                phi_r = potential(x_val, M500_SI, R500_SI, g500, c500, 1.0); // GR
            }
            if (model == 1) { // f(R)
                phi_r = potential(x_val, M500_SI, R500_SI, g500, c500, 4.0/3.0); // f(R)
            }
            if (model == 2) { // DGP
                // densità NFW in unità di rho_s, quindi adimensionale
                auto integrand_mass = [&](double xp) -> double {
                return rho_NFW(xp, R500_SI, M500_SI, c500, g500) / rho_NFW(1.0, R500_SI, M500_SI, c500, g500) * xp * xp;
                };

                double epsabs = 1e-6, epsrel = 1e-6;
                int limit_size = 5000, rule = 6;

                double integral_mass = wrapper::gsl::GSL_integrate_qag(integrand_mass, 0.0, x_val, epsabs, epsrel, limit_size, rule);

                // massa relativa (adimensionale)
                double M_in_x_over_M500 = integral_mass; 
                
                wrapper::gsl::GSL_integrate_qag(integrand_mass, 0.0, 1.0, epsabs, epsrel, limit_size, rule);

                // DGP factor
                double y = pow((0.32 * 500), -1.0/3.0) * x_val* pow(M_in_x_over_M500, -1.0/3.0);
                double g_y = pow(y, 3) * (sqrt(1.0 + pow(y, -3)) - 1.0);
                double k_model = 1.0 + (2.0 / (3.0 * (1.15))) * g_y;
                phi_r = potential(x_val, M500_SI, R500_SI, g500, c500, k_model);
            }
            return phi_r;
        };
    

    // Tolleranze numeriche
    double epsrel = 1e-1;
    double step = 0.0001;  

    // Calcola la derivata numerica con GSL
    double deriv = cbl::wrapper::gsl::GSL_derivative(phi_func, x, step, epsrel);

    return deriv; // [m/s²]
    };


    // Densità NFW
    auto rho_r = [&](double x) -> double {
        return (M500_SI * std::pow(c500, 2) * g500) /
               (4.0 * M_PI * std::pow(R500_SI, 3) * x * std::pow(1 + c500 * x, 2));
    };

    // Integrando in unità adimensionali (x = r/R500)
    auto integrand = [&](double x) -> double {
        if (x <= x_perp) return 0.0;
        double num = (x*x - x_perp*x_perp);
        double denom = std::sqrt(x*x - x_perp*x_perp);
        return num/denom * (+dphi_dr(x)) * rho_r(x);
    };

    // Parametri di integrazione GSL
    double epsabs = 1e-6;
    double epsrel = 1e-6;
    int limit_size = 5000;

    double integral = wrapper::gsl::GSL_integrate_qagiu(
        integrand, x_perp, epsabs, epsrel, limit_size
    );

    // Calcolo della superficie Σ(x_perp)
    double Sigma = sigma(x_perp, c500, R500_SI, M500_SI);

    // Formula finale
    double delta_TD = (3.0 * R500_SI/ (c * Sigma)) * integral; // [m/s]
    return delta_TD;
}

// ==========================================================
//  Media su un intervallo di massa [Mmin, Mmax]
// ==========================================================
double deltaTD_avg(double x_phys, double Mmin, double Mmax, double H0, double z, int model)
{
    // Limiti in y = ln M
    //double y0 = std::log(Mmin);
    //double y1 = std::log(Mmax);
    
    // INTEGRANDO NUMERATORE: Sigma(M) * Delta_TD(M) * M^{-alpha} dM
    auto integrand_num = [&](double M) -> double {
        //double M = std::exp(y);              // M in M_sun
        ClusterParams p = params_from_M(M, H0, z);
        double x_perp = x_phys/p.R500; // adimensionale
        double Sigma = Sigma_phys_for_M(x_perp, p); // kg/m^2
        double DeltaTD = delta_TD(x_perp, p.R500_SI, p.M500_SI, p.c500, model); // m/s
        // dM = M dy -> integrand in y: Sigma * DeltaTD * M^{-alpha} * M
        double weight = mf.mass_function(M/h, z, "Tinker", "EisensteinHu", false, "none", 500.0);
        return Sigma * DeltaTD * weight;
    };

    // INTEGRANDO DENOMINATORE: Sigma(M) * M^{-alpha} dM
    auto integrand_den = [&](double M) -> double {
        //double M = std::exp(y);              // M in M_sun
        ClusterParams p = params_from_M(M, H0, z);
        double x_perp = x_phys/p.R500; // adimensionale
        double Sigma = Sigma_phys_for_M(x_perp, p); // kg/m^2
        double weight = mf.mass_function(M/h, z, "Tinker", "EisensteinHu", false, "none", 500.0);
        return Sigma * weight;
    };

    // Parametri GSL (tolleranze)
    double epsabs = 1e-4;   // se vuoi più preciso riduci (es. 1e-6)
    double epsrel = 1e-4;
    int limit_size = 5000;
    int rule = 6;

    // Calcolo integrali (nota: wrapper prende funzioni (double)->double)
    double num = wrapper::gsl::GSL_integrate_qag(integrand_num, Mmin, Mmax, epsabs, epsrel, limit_size, rule);
    double den = wrapper::gsl::GSL_integrate_qag(integrand_den, Mmin, Mmax, epsabs, epsrel, limit_size, rule);

    if (den == 0.0) {
        // evita divisione per zero
        return 0.0;
    }
    return num / den; // [m/s]
}
// ============================================
// Surface Brightness Effect (Δ_SB)
// ============================================
//
// Formula: Δ_SB = -(10/3) * <δ(z)> * Δ_TD
// con <δ(z)> = [3 + α(z)] * β
//
// Parametri tipici:
// α(z) ≈ 2.0   → slope luminosity function
// β ≈ 1.5      → slope della funzione n(>L)
// ============================================

double delta_SB(double Delta_TD) //double alpha_z = 2.0, double beta = 1.5)
{
    double delta_z_mean = 0.516;              //(3.0 + alpha_z) * beta;   calcolo diretto del <δ(z)>
    return - (10.0 / 3.0) * delta_z_mean * Delta_TD;
}




int main() {
    int N = 1000;
    // Input
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

    //creo una griglia
    double r[N];
    double delta_r = 3.0/N;
    for (int i=0; i<N-1; i++) {
        r[i+1]=r[i]+delta_r;
    };

    double r_tilde[N]; 
    double g500 = pow(log(1+c500)-(c500/(1+c500)), -1);
    for (int i=1; i<N-1; i++) {
        r_tilde[i] = r[i]/R500;
    }

    //Calcolo il velocity offset di un singolo cluster (M=10^14 M_sun) con diversi modelli di gravità
    std::ofstream fout1("delta_single_cluster_gsl.dat");

    // Potenziali centrali per i diversi modelli
    double phi0_GR = -g500*c500*G*M500_SI/R500_SI;
    double phi0_fR = (4.0/3.0)*phi0_GR;
    double phi0_DGP = phi0_GR;

    // Ciclo sui raggi
    #pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < N - 1; i++) {
        double x_perp = r_tilde[i];
        
        // Calcolo per ciascun modello di gravità
        double delta_GR  = delta(x_perp, R500_SI, phi0_GR, M500_SI, c500, 0);
        double delta_fR  = delta(x_perp, R500_SI, phi0_fR, M500_SI, c500, 1);
        double delta_DGP = delta(x_perp, R500_SI, phi0_DGP, M500_SI, c500, 2);

        // Conversione in km/s
        #pragma omp critical
        {
        fout1 << x_perp << "   "
            << delta_GR/1000.0 << "   "
            << delta_fR/1000.0 << "   "
            << delta_DGP/1000.0 << "\n";
        }

    };
    // chiudo il file
    fout1.close();
    std::cout << "File 'delta_single_cluster_gsl.dat' scritto con successo." << std::endl;
    
    //calcolo il velocity offset per una popolazione di cluster
    std::ofstream fout2("delta_population_avg_gsl.dat");
    //fout2 << "x_perp  Delta_population_GR  Delta_population_f(R)  Delta_population_DGP\n";
    double Mmin = 1e13; // M_sun
    double Mmax = 1e16; // M_sun
    //double alpha = 1.0; // slope della funzione di massa
    double Delta_pop_GR[N];
    double Delta_pop_fR[N];
    double Delta_pop_DGP[N];

    #pragma omp parallel for schedule(dynamic)
    for (int i=1; i<N-1; i++) {
        double x_perp = r[i];
        Delta_pop_GR[i] = delta_population_avg(x_perp, Mmin, Mmax, H0, z, 0);
        Delta_pop_fR[i] = delta_population_avg(x_perp, Mmin, Mmax, H0, z, 1);
        Delta_pop_DGP[i] = delta_population_avg(x_perp, Mmin, Mmax, H0, z, 2);
    }

    for (int i = 1; i < N - 1; i++) {
        double x_perp = r[i];
        fout2 << x_perp << "   "
          << Delta_pop_GR[i]/1000.0 << "   "
          << Delta_pop_fR[i]/1000.0 << "   "
          << Delta_pop_DGP[i]/1000.0 << "\n"; //il /1000 converte in Km/s
    
    }


    fout2.close();
    std::cout << "File 'delta_population_avg_gsl.dat' scritto con successo." << std::endl;
    // ==========================
    // TD EFFECT (parallelizzato)
    // ==========================
    std::vector<double> DeltaTD_GR(N), DeltaTD_fR(N), DeltaTD_DGP(N);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < N - 1; i++) {
        double x_perp = r[i];
        DeltaTD_GR[i]  = deltaTD_avg(x_perp, Mmin, Mmax, H0, z, 0);
        DeltaTD_fR[i]  = deltaTD_avg(x_perp, Mmin, Mmax, H0, z, 1);
        DeltaTD_DGP[i] = deltaTD_avg(x_perp, Mmin, Mmax, H0, z, 2);
    }   

    // scrittura su file (sequenziale, ordinata)
    std::ofstream fout3("delta_TD.dat");
    for (int i = 1; i < N - 1; i++) {
        fout3 << r_tilde[i] << "   "
          << DeltaTD_GR[i]/1000.0 << "   "
          << DeltaTD_fR[i]/1000.0 << "   "
          << DeltaTD_DGP[i]/1000.0 << "\n";
    } 
    fout3.close();
    std::cout << "File 'delta_TD.dat' scritto con successo." << std::endl;

    // ==========================
    // LC + SB EFFECT (parallelizzato)
    // ==========================
    std::vector<double> DeltaLC_GR(N), DeltaLC_fR(N), DeltaLC_DGP(N);
    std::vector<double> DeltaSB_GR(N), DeltaSB_fR(N), DeltaSB_DGP(N);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < N - 1; i++) {
        DeltaLC_GR[i]  = 2.0/3.0 * DeltaTD_GR[i];
        DeltaLC_fR[i]  = 2.0/3.0 * DeltaTD_fR[i];
        DeltaLC_DGP[i] = 2.0/3.0 * DeltaTD_DGP[i];

        DeltaSB_GR[i]  = delta_SB(DeltaTD_GR[i]);
        DeltaSB_fR[i]  = delta_SB(DeltaTD_fR[i]);
            DeltaSB_DGP[i] = delta_SB(DeltaTD_DGP[i]);
    }

    // scrittura su file
    std::ofstream fout4("delta_SB.dat");
    for (int i = 1; i < N - 1; i++) {
        fout4 << r_tilde[i] << "   "
          << (DeltaSB_GR[i]+DeltaLC_GR[i])/1000.0 << "   "
          << (DeltaSB_fR[i]+DeltaLC_fR[i])/1000.0 << "   "
          << (DeltaSB_DGP[i]+DeltaLC_DGP[i])/1000.0 << "\n";
    }
    fout4.close();
    std::cout << "File 'delta_SB.dat' scritto con successo." << std::endl;

    // ==========================
        // COMPLETE VELOCITY OFFSET
    // ==========================
    std::ofstream fout5("delta_complete_population_gsl.dat");
    for (int i = 1; i < N - 1; i++) {
        fout5 << r_tilde[i] << "   "
          << (Delta_pop_GR[i] + DeltaTD_GR[i] + DeltaLC_GR[i] + DeltaSB_GR[i]) / 1000.0 << "   "
          << (Delta_pop_fR[i] + DeltaTD_fR[i] + DeltaLC_fR[i] + DeltaSB_fR[i]) / 1000.0 << "   "
          << (Delta_pop_DGP[i] + DeltaTD_DGP[i] + DeltaLC_DGP[i] + DeltaSB_DGP[i]) / 1000.0 << "\n";
    }
    fout5.close();
    std::cout << "File 'delta_complete_population_gsl.dat' scritto con successo." << std::endl;
    return 0;
};