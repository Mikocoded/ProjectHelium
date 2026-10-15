#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <QApplication>
#include <QHBoxLayout>
#include <QWidget>
#include <QtDataVisualization/Q3DScatter>
#include <QtDataVisualization/QScatter3DSeries>
#include <QtDataVisualization/QScatterDataProxy>
#include <QtDataVisualization/QValue3DAxis>
#include <QColor>
#include <QVector3D>
#include <algorithm>
#include <QtDataVisualization/Q3DTheme> 
#include <random>

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

QColor colormap(double t) { 
    static const double stops[5][3] = {
        {0, 0, 4}, {87, 16, 110}, {188, 55, 84}, {249, 142, 9}, {252, 255, 164}
    };
    t = clamp(t, 0.0, 1.0) * 4.0;
    int i = min(int(t), 3);
    double f = t - i;
    return QColor(int(stops[i][0] + f * (stops[i+1][0] - stops[i][0])),
                  int(stops[i][1] + f * (stops[i+1][1] - stops[i][1])),
                  int(stops[i][2] + f * (stops[i+1][2] - stops[i][2])));
}


int main(int argc, char *argv[]) {
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
    int Nr = 1e5;
    double h = (rmax - rmin) / (20*Nr - 1);

    vector<double> u = numerovRadial(Nr, h, rmin, l, Z, E);

    double normConst = 0.0;
    for (int i = 0; i < Nr; i++) normConst += u[i] * u[i] * h;
    normConst = sqrt(normConst);
    for (int i = 0; i < Nr; i++) u[i] /= normConst;

    vector<double> R(Nr);
    for (int i = 0; i < Nr; i++) {
        double r = rmin + i * h;
        R[i] = u[i] / r;
    };

    struct Pt { float x, y, z, prob; };
    vector<Pt> pts;

    int N = 1e5;

    vector<double> cdf(Nr);
    double acc = 0.0;
    for (int i = 0; i < Nr; i++) {
        double r = rmin + i * h;
        acc += R[i] * R[i] * r * r;
        cdf[i] = acc;
    }
    for (double &c : cdf) c /= acc;


    double fmax = 0.0;
    for (int i=0; i < 2000; i++){
        double th = (i + 0.5) * M_PI/2000;
        fmax = max(fmax, norm(sphericalHarmonic(l, m, th, 0.0))*sin(th));
    }

    double radInt = 0.0;
    for (int i=0; i< Nr; i++){
        double r = rmin + i * h;
        radInt += R[i]*R[i]*r*r*h;
    }
    
    double angInt = 0.0;
    double Nth=2000;
    double dth = M_PI/Nth;
    for (int k=0; k<Nth; k++){
        double th = (k+0.5)*dth;
        angInt += norm(sphericalHarmonic(l,m,th,0.0)) * sin(th) * dth * 2*M_PI;    
    }
    double totalProb = radInt * angInt;

    double Pcum = 0.0;
    double r0 = 1.0;
    for (int i=0; i < Nr; i++){
        double r = rmin + i * h;
        if (r>r0) break;
        Pcum += R[i]*R[i]*r*r*h;
    }


    cout.setf(ios::fixed);
    cout.precision(6);
    cout << "=== Hydrogenlike model (He, one electron approximation) ===\n";
    cout << "n=" << n << " l=" << l << " m=" << m << "  Z_eff=" << Z << "\n";
    cout << "Energia E = " << E << " Hartree\n";
    cout << "Entire probability sum (should be around 1.0): "
    << totalProb << "\n\n";

    cout << "Probability of finding an electron in r < a0: " << Pcum << "\n\n";

    mt19937 rng(42);
    uniform_real_distribution<double> U(0.0, 1.0);

    double maxProb = 0.0;
    ofstream fout("HeProbability.csv");
    fout << "x,y,z,prob,energy" << "\n";
    pts.reserve(N);
    while((int)pts.size()< N){
        double q = U(rng);
        int i = lower_bound(cdf.begin(), cdf.end(), q) - cdf.begin();
        if (i==0) i=1;
        double f = (q- cdf[i-1])/(cdf[i] - cdf[i-1]);
        double r = rmin + (i - 1 + f)*h;

        double theta;
        do { theta = U(rng) * M_PI; }
        while (U(rng) * fmax > norm(sphericalHarmonic(l, m, theta, 0.0)) * sin(theta));
        double phi = U(rng) * 2.0 * M_PI;

        double Prob = R[i] * R[i] * norm(sphericalHarmonic(l, m, theta, phi));
        double x = r * sin(theta) * cos(phi);
        double y = r * sin(theta) * sin(phi);
        double z = r * cos(theta);

        pts.push_back({(float)x, (float)y, (float)z, (float)Prob});
        maxProb = max(maxProb, Prob);

        fout << x << "," << y << "," << z << "," << Prob << "," << E << "\n";
    }

    fout.close();
    cout << "Saved to HeProbability.csv\n";

    int vecsiz = pts.size();
    cout << vecsiz << endl;

    int plotmaybe;

    cout << "Would you like to open Qt plot? (press 0)" << endl;
    cin >> plotmaybe;

    if(plotmaybe == 0){
        QApplication app(argc, argv);

        const int nBins = 24;     

        vector<QScatterDataArray*> bins(nBins);
        for (auto &b : bins) { b = new QScatterDataArray; }

        vector<float> pv;
        pv.reserve(pts.size());
        for (const Pt &p : pts) pv.push_back(p.prob);
        size_t k = size_t(pv.size() * 0.99);
        nth_element(pv.begin(), pv.begin() + k, pv.end());
        double pRef = pv[k];
        pv.clear();
        pv.shrink_to_fit();

        const double gamma = 0.4;

        float extent = 0.0f;
        for (const Pt &p : pts) {
            double t = pow(min(p.prob / pRef, 1.0), gamma);
            int b = min(int(t * nBins), nBins - 1);
            bins[b]->append(QScatterDataItem(QVector3D(p.x, p.z, p.y)));
            extent = max({extent, fabs(p.x), fabs(p.y), fabs(p.z)});
        }

        auto *graph = new Q3DScatter();
        graph->setShadowQuality(QAbstract3DGraph::ShadowQualityNone);
        graph->activeTheme()->setBackgroundEnabled(true);

        graph->activeTheme()->setType(Q3DTheme::ThemeQt);
        graph->activeTheme()->setWindowColor(Qt::black);
        graph->activeTheme()->setBackgroundColor(Qt::black);

        graph->setAspectRatio(1.0);
        graph->setHorizontalAspectRatio(1.0);

        for (int b = 0; b < nBins; b++) {
            if (bins[b]->isEmpty()) { delete bins[b]; continue; }
            auto *series = new QScatter3DSeries;
            series->dataProxy()->resetArray(bins[b]); 
            series->setMesh(QAbstract3DSeries::MeshSphere);
            series->setItemSize(0.02f);
            series->setColorStyle(Q3DTheme::ColorStyleUniform);
            series->setBaseColor(colormap((b + 0.5) / nBins));
            graph->addSeries(series);
        }

        QWidget *container = QWidget::createWindowContainer(graph);
        QWidget window;
        window.resize(1600, 1200);
        auto *layout = new QHBoxLayout(&window);
        layout->addWidget(container);
        window.show();

        return app.exec();
    }
}
