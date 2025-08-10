#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QLineEdit>
#include <QComboBox>
#include <QPushButton>
#include <QChartView>
#include <QLineSeries>
#include <QFormLayout>
#include <vector>


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE



class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void simulate();

private:
    QLineEdit *amplitudeInput;
    QLineEdit *frequencyInput;
    QLineEdit *capacitance1Input;
    QLineEdit *capacitance2Input;
    QLineEdit *resistanceInput;
    QComboBox *diodeTypeInput;
    QChart *chart;
    QChartView *chartView;

    QLineEdit *timeInput;


    void setupUI();
    std::pair<std::vector<double>, std::vector<double>> solveODE(double A, double f, double C1, double C2, double R, double dt, double T);

};

#endif // MAINWINDOW_H
