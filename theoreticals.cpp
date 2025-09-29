#include <iostream>
#include <cmath>

int main () {
    double Msol = 1.989*pow(10, 30);
    double M = pow(10., 14);
    double pi = M_PI;
    double Mpc_m = 3.0856*pow(10, 22);
    double H_0 = (70*Mpc_m)/1000;  // così l'unità di misura è 1/s
    double G = 6.674*pow(10,-11);
    double delta_r = (3.0/1000.0)*Mpc_m; //cambia le unità di misura (ora è in m)
    double r[1000];
    double dV[1000];
    double rho[1000];
    double rho_c = (3*H_0)/(8*pi*G);
    double r_500;
    double r_200;
    r[0] = 0;
    dV[0] = 0;
    rho[0] = M/(4/3*pi*pow(delta_r, 3));

    for (int i=0; i<1000; i++) {       //this for cycle is just to calculate a density to obtain a value of r200 and r500
        r[i+1] = r[i]+delta_r;
        dV[i+1] = dV[i]+(4/3*pi*pow(r[i], 3));
        rho[i+1] = M/dV[i+1];
        if (rho[i]>450*rho_c && rho[i]<550*rho_c) {
                r_500 = r[i];
        };
        if (rho[i]>150*rho_c && rho[i]<250*rho_c) {
            r_200 = r[i];
        };
    };
    
    double M_500 = 4/3*pi*500*rho_c*pow(r_500, 3);
    double M_200 = 4/3*pi*200*rho_c*pow(r_200, 3);
    double c_200 = 0.63*pow(M_200/M, -0.08);    //Duffy et al using constant redshift z=1
    double c_500 = c_200/pow(5/2*(M_200/M_500), 1/3);
    double nfw[1000];
    double r_tilde[1000];
    double g_500 = pow(log(1+c_500)-(c_500/(1+c_500)), -1);
    double r200[10];
    double r500[10];
    double M200[10];
    double M500[10];
    double c200[10];    
    double c500[10];
    double NFW[1000]; //redefine the nfw profile using the recalculated values

    for (int j=0; j<10; j++) {
        for (int i=0; i<1000; i++) {
            r_tilde[i] = r[i]/r_500;
            nfw[i] = (M_500*pow(c_500, 2)*g_500)/(4*pi*pow(r_500, 3)*r_tilde[i]*(1+r_tilde[i]));
        //recalculating r200 and r500 using the new density profile
            if (nfw[i]>170*rho_c && nfw[i]<230*rho_c) {
            r200[j] = r_tilde[i]*r_500;
            };
            if (nfw[i]>470*rho_c && nfw[i]<530*rho_c) {
            r500[j] = r_tilde[i]*r_500;
            };
        };
    //redefining M200, M500, c200, c500 using the new values of r200 and r500 
        M200[j] = 4.0/3.0*pi*200*rho_c*pow(r200[j], 3);
        M500[j] = 4.0/3.0*pi*500*rho_c*pow(r500[j], 3);
        c200[j] = 0.63*pow(M200[j]/M, -0.08);    //Duffy et al using again constant redshift z=1
        c500[j] = c200[j]/pow(5/2*(M200[j]/M500[j]), 1/3);
        
        for (int i=0; i<1000; i++) {
            NFW[i] = (M500[j]*pow(c500[j], 2)*g_500)/(4*pi*pow(r500[j], 3)*r_tilde[i]*(1+r_tilde[i]));
        };
    };

    std::cout << "rho_c = " << rho_c << " Kg/m^3" << "\n";
    std::cout << "rho1 = " << rho[1] << " Kg/m^3" << "\n";
    std::cout << "r_200 = " << r_200 << " Mpc" << "\n";
    std::cout << "r_500 = " << r_500 << " Mpc" << "\n";

}