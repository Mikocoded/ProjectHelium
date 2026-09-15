#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

// =====================================================================
// Jednostki atomowe: hbar = m_e = e = 4*pi*eps0 = 1
// Model: przyblizenie jednoelektronowe (wodoropodobne) z efektywnym
// ladunkiem jadra Z_eff. Dla helu:
//   Z = 1.6875   -> standardowy wariacyjny ladunek efektywny dla 1s He
// Prawdziwy atom helu (2 elektrony, oddzialywanie e-e) wymaga metod
// wieloelektronowych (Hartree-Fock, CI, QMC) - patrz komentarz na koncu.
// =====================================================================

// ---------- Stowarzyszone wielomiany Legendre'a P_l^m(x), m>=0 ----------
// (konwencja z Condon-Shortley phase, zgodna z definicja harmonik sferycznych)
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

// ---------- Znormalizowane harmoniki sferyczne Y_l^m(theta,phi) ----------
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

// ---------- Czesc radialna: metoda Numerova dla u(r) = r*R(r) ----------
// Rownanie:  u''(r) = f(r) * u(r)
// f(r) = 2*(V(r) - E) + l(l+1)/r^2,   V(r) = -Z/r
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
    // ---------------- Parametry orbitalu ----------------
    int n = 3;          // glowna liczba kwantowa
    int l = 2;          // orbitalna liczba kwantowa (0=s,1=p,2=d,...)
    int m = 1;          // magnetyczna liczba kwantowa
    double Z = 1.6875;  // efektywny ladunek jadra dla 1s helu (wariacyjny)

    if (l >= n || abs(m) > l) {
        cerr << "Niepoprawne liczby kwantowe (wymagane: l<n, |m|<=l)\n";
        return 1;
    }

    // Energia stanu wodoropodobnego (dokladna dla czystego potencjalu Coulomba)
    double E = -(Z * Z) / (2.0 * n * n);

    // ---------------- Siatka radialna ----------------
    double rmin = 1e-4;
    double rmax = 20.0 / Z * n * n; 
    int Nr = 4000;
    double h = (rmax - rmin) / (Nr - 1);

    vector<double> u = numerovRadial(Nr, h, rmin, l, Z, E);

    // Normalizacja: całka |u(r)|^2 dr = 1  <=>  całka |R(r)|^2 r^2 dr = 1
    double normConst = 0.0;
    for (int i = 0; i < Nr; i++) normConst += u[i] * u[i] * h;
    normConst = sqrt(normConst);
    for (int i = 0; i < Nr; i++) u[i] /= normConst;

    vector<double> R(Nr);
    for (int i = 0; i < Nr; i++) {
        double r = rmin + i * h;
        R[i] = u[i] / r;
    }

    // ---------------- Sprawdzenie normalizacji calej funkcji falowej ----------------
    // Psi(r,theta,phi) = R(r) * Y_l^m(theta,phi)
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
    cout << "=== Model wodoropodobny (He, przyblizenie 1-elektronowe) ===\n";
    cout << "n=" << n << " l=" << l << " m=" << m << "  Z_eff=" << Z << "\n";
    cout << "Energia E = " << E << " Hartree\n";
    cout << "Calkowite prawdopodobienstwo (powinno wynosic ~1.0): "
         << totalProb << "\n\n";

    // ---------------- Przykladowe obliczenie prawdopodobienstwa ----------------
    // P(r < r0) = całka_0^r0 |R(r)|^2 r^2 dr  (gestosc radialna)
    double r0 = 1.0; // promien w jednostkach a0
    double Pcum = 0.0;
    for (int i = 0; i < Nr; i++) {
        double r = rmin + i * h;
        if (r > r0) break;
        Pcum += R[i] * R[i] * r * r * h;
    }
    cout << "Prawdopodobienstwo znalezienia elektronu w r < " << r0
         << " a0: " << Pcum << "\n\n";

    // ---------------- Zapis gestosci radialnej r^2|R(r)|^2 do pliku ----------------
    ofstream fout("Heprobability.csv");
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
    cout << "Zapisano prawdopodobienstwo do pliku Heprobability.csv\n";

    return 0;
}