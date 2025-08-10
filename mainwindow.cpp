#include "mainwindow.h"
#include <QVBoxLayout>
#include <cmath>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), chart(new QChart()), chartView(new QChartView(chart)) {
    setupUI();
}

MainWindow::~MainWindow() {}

void MainWindow::setupUI() {
    auto *centralWidget = new QWidget(this);
    auto *mainLayout = new QVBoxLayout(centralWidget);


    auto *formLayout = new QFormLayout();

    //fasta startvärden
    amplitudeInput = new QLineEdit("10");
    frequencyInput = new QLineEdit("50");
    capacitance1Input = new QLineEdit("0.001");
    capacitance2Input = new QLineEdit("0.001");
    resistanceInput = new QLineEdit("1000");
    diodeTypeInput = new QComboBox();
    diodeTypeInput->addItems({"Schottky", "Silicon"});

    formLayout->addRow("Amplitude (V):", amplitudeInput);
    formLayout->addRow("Frequency (Hz):", frequencyInput);
    formLayout->addRow("Capacitor C1 (F):", capacitance1Input);
    formLayout->addRow("Capacitor C2 (F):", capacitance2Input);
    formLayout->addRow("Resistor R (Ohm):", resistanceInput);
    formLayout->addRow("Diode Type:", diodeTypeInput);

    timeInput = new QLineEdit("0.1");
    formLayout->addRow("Analysis Time (s):", timeInput);


    mainLayout->addLayout(formLayout);

    // Button to start simulation
    auto *simulateButton = new QPushButton("Simulate");
    mainLayout->addWidget(simulateButton);
    connect(simulateButton, &QPushButton::clicked, this, &MainWindow::simulate);

    // Chart setup
    chartView->setRenderHint(QPainter::Antialiasing);
    mainLayout->addWidget(chartView);

    setCentralWidget(centralWidget);
}

void MainWindow::simulate() {
    const auto A = amplitudeInput->text().toDouble();
    const auto f = frequencyInput->text().toDouble();
    const auto C1 = capacitance1Input->text().toDouble();
    const auto C2 = capacitance2Input->text().toDouble();
    const auto R = resistanceInput->text().toDouble();

    const auto dt = 1e-5;// testa sänk höj tidssteget
    const auto T = timeInput->text().toDouble();

    auto [v1, v2] = solveODE(A, f, C1, C2, R, dt, T);

    //std::vector<double> v2 = solveODE(A, f, C1, C2, R, dt, T);
    //std::vector<double> v1 = solveODE(A, f, C1, C2, R, dt, T);

    // Plotta resulat
    auto *seriesv1 = new QLineSeries();
    auto *seriesv2 = new QLineSeries();
    auto *seriesvs = new QLineSeries();

    seriesv1->setName("v1 (V)");
    seriesv2->setName("v2 (V)");
    seriesvs->setName("vs (V)");

    for (size_t i = 0; i < v1.size(); ++i) {
        const auto t = i * dt;
        const auto Vs = A * std::sin(2 * M_PI * f * t);

        seriesv1->append(t, v1[i]);
        seriesv2->append(t, v2[i]);
        seriesvs->append(t, Vs);
    }

    chart->removeAllSeries();
    chart->addSeries(seriesvs);
    chart->addSeries(seriesv1);
    chart->addSeries(seriesv2);

    chart->createDefaultAxes();


    chart->axes(Qt::Horizontal).first()->setTitleText("Tid (s)");
    chart->axes(Qt::Vertical).first()->setTitleText("Spänning (V)");


    chart->setTitle("Resultat från simuleringen");


    chart->legend()->setVisible(true);
    chart->legend()->setAlignment(Qt::AlignBottom);
}


std::pair<std::vector<double>, std::vector<double>> MainWindow::solveODE(
    double A, double f, double C1, double C2, double R, double dt, double T)
{
    const auto N = static_cast<int>(T / dt);
    std::vector<double> v1(N, 0.0), v2(N, 0.0);

    // Diodparametrar
    struct D { double Is, n, Vt, Gmin; };
        const auto d1 = (diodeTypeInput->currentText()=="Schottky") ? D{1e-7, 1.05, 0.02585, 1e-8} : D{1e-12, 2.0,  0.02585, 1e-8};
        const auto d2 = d1;

    auto I = [](double vd, const D& d){
        return d.Is*(std::exp(vd/(d.n*d.Vt)) - 1.0) + d.Gmin*vd;
    };
    auto G = [](double vd, const D& d){
        return (d.Is/(d.n*d.Vt))*std::exp(vd/(d.n*d.Vt)) + d.Gmin;
    };

    for (int k=0; k<N-1; ++k) {
        // Backward Euler
        const double t_k  = k*dt;
        const double t_k1 = (k+1)*dt;
        const double vs_k  = A*std::sin(2*M_PI*f*t_k);
        const double vs_k1 = A*std::sin(2*M_PI*f*t_k1);

        auto v1n = v1[k], v2n = v2[k]; // startgissning

        for (int it=0; it<40; ++it) {

            const auto vd1 = -v1n;
            const auto Id1 = I(vd1, d1);
            const auto Gd1 = G(vd1, d1);


            const auto vd2 = v1n - v2n;
            const auto Id2 = I(vd2, d2);
            const auto Gd2 = G(vd2, d2);

            // KCL/Backward Euler

            const auto F1 = C1 * ( (v1n - vs_k1) - (v1[k] - vs_k) )/dt - Id1 + Id2;


            const auto F2 = (v2n/R) + C2 * ( (v2n - v2[k]) / dt ) - Id2;

            // Newton/Jacobian
            const auto dF1dv1 = C1/dt + Gd1 + Gd2;
            const auto dF1dv2 = -Gd2;
            const auto dF2dv1 = -Gd2;
            const auto dF2dv2 = (1.0/R) + C2/dt + Gd2;

            const auto det = dF1dv1*dF2dv2 - dF1dv2*dF2dv1;
            if (std::fabs(det) < 1e-20) break;

            const auto rhs1 = -F1, rhs2 = -F2;
            const auto dv1  = ( rhs1*dF2dv2 - dF1dv2*rhs2) / det;
            const auto dv2  = (-rhs1*dF2dv1 + dF1dv1*rhs2) / det;

            v1n += dv1; v2n += dv2;

        }

        v1[k+1] = v1n;
        v2[k+1] = v2n;
    }

    return {v1, v2};
}























