#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

double assocLegendre(int l, int m, double x) {
    double pmm = 1.0;
    if (m > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m) return pmm;
    double pmmp1 = x * (2 * m + 1) * pmm;
    if (l == m + 1) return pmmp1;
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ll++) {
        pll = ((2 * ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

complex<double> sphericalHarmonic(int l, int m, double theta, double phi) {
    int am = abs(m);
    double x = cos(theta);
    double Plm = assocLegendre(l, am, x);
    double norm = sqrt((2.0 * l + 1.0) / (4.0 * M_PI) *
    tgamma(l - am + 1.0) / tgamma(l + am + 1.0));
    complex<double> Y = norm * Plm * polar(1.0, am * phi);
    if (m < 0) Y = pow(-1.0, am) * conj(Y);
    return Y;
}

double radialF(double r, double E, int l, double Z) {
    double V = -Z / r;
    return 2.0 * (V - E) + l * (l + 1) / (r * r);
}

vector<double> numerovRadial(int Npoints, double h, double rmin, int l, double Z, double E) {
    vector<double> u(Npoints, 0.0);
    double kappa = sqrt(-2.0*E);
    double rN1 = rmin + (Npoints-1)*h;
    double rN2 = rmin + (Npoints-2)*h;
    u[Npoints-1] = exp(-kappa*rN1);
    u[Npoints-2] = exp(-kappa*rN2);

    for (int i = Npoints-2; i >= 1; i--) {
        double r_im1 = rmin + (i-1)*h;
        double r_i   = rmin + i*h;
        double r_ip1 = rmin + (i+1)*h;
        double f_im1 = radialF(r_im1, E, l, Z);
        double f_i   = radialF(r_i, E, l, Z);
        double f_ip1 = radialF(r_ip1, E, l, Z);
        double num = 2.0*(1.0 - 5.0/12.0*h*h*f_i)*u[i] - (1.0 + h*h/12.0*f_ip1)*u[i+1];
        double den = 1.0 + h*h/12.0*f_im1;
        u[i-1] = num/den;
    }
    return u;
}


int main() {
    int n = 3;        
    int l = 2;        
    int m = 1;          
    double Z = 2;  

    cout << "Specify quantum numbers" << endl;
    cout << "Principal quantum number (n):" << endl;
    cin >> n;
    cout << "Angular momentum quantum number (l):" << endl;
    cin >> l;
    cout << "Magnetic quantum number (m):" << endl;
    cin >> m;

    if (l >= n || abs(m) > l) {
        cerr << "Incorrect quantum numbers (required: l<n, |m|<=l)\n";
        return 1;
    }

    double E = -(Z * Z) / (2.0 * n * n);

    double rmin = 1e-4;
    double rmax = 20.0 / Z * n * n;
    int Nr = 4000;
    double h = (rmax - rmin) / (Nr - 1);

    vector<double> u = numerovRadial(Nr, h, rmin, l, Z, E);

    double normConst = 0.0;
    for (int i = 0; i < Nr; i++) normConst += u[i] * u[i] * h;
    normConst = sqrt(normConst);
    for (int i = 0; i < Nr; i++) u[i] /= normConst;

    vector<double> R(Nr);
    for (int i = 0; i < Nr; i++) {
        double r = rmin + i * h;
        R[i] = u[i] / r;
    }

    int Ntheta = 60, Nphi = 60;
    double dtheta = M_PI / (Ntheta - 1);
    double dphi = 2.0 * M_PI / (Nphi - 1);

    double totalProb = 0.0;
    for (int i = 0; i < Nr; i++) {
        double r = rmin + i * h;
        double Rr2 = R[i] * R[i];
        double shellFactor = r * r * h;
        double angularSum = 0.0;
        for (int it = 0; it < Ntheta; it++) {
            double theta = it * dtheta;
            double s = sin(theta);
            for (int ip = 0; ip < Nphi; ip++) {
                double phi = ip * dphi;
                complex<double> Ylm = sphericalHarmonic(l, m, theta, phi);
                angularSum += norm(Ylm) * s * dtheta * dphi;
            }
        }
        totalProb += Rr2 * shellFactor * angularSum;
    }

    cout.setf(ios::fixed);
    cout.precision(6);
    cout << "=== Hydrogenlike model (He, one electron approximation) ===\n";
    cout << "n=" << n << " l=" << l << " m=" << m << "  Z_eff=" << Z << "\n";
    cout << "Energia E = " << E << " Hartree\n";
    cout << "Entire probability sum (should be around ~1.0): "
    << totalProb << "\n\n";

    double r0 = 1.0;
    double Pcum = 0.0;
    for (int i = 0; i < Nr; i++) {
        double r = rmin + i * h;
        if (r > r0) break;
        Pcum += R[i] * R[i] * r * r * h;
    }
    cout << "Probability of finding an electron in r < a0: " << Pcum << "\n\n";

    ofstream fout("HeProbability.csv");
    fout << "r,theta,phi,prob,energy" << "\n";
    for (int i = 0; i < Nr; i += 4) {
        double r = rmin + i * h;
        for (int it = 0; it<Ntheta; it++){
            double theta = it * dtheta;
            for (int ip=0; ip<Nphi; ip++){
                double phi = ip *dphi;

                complex<double> Ylm = sphericalHarmonic(l, m, theta, phi);
                double Y_Sq = norm(Ylm);
                double Prob = (R[i] * R[i]) * Y_Sq;

                fout << r << "," << theta << "," << phi << "," << Prob << "," << E << "\n";
            }
        }
    }

    fout.close();
    cout << "Saved to HeProbability.csv\n";

    return 0;
}
