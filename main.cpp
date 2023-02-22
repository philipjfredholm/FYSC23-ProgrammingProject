#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include<unistd.h> //Replace with #include <windows.h> if not running Linux.

using namespace Eigen;

//ROOT headers
#include <TRint.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>

MatrixXd makeReal(MatrixXcd& complexMatrix) {
    int rows = complexMatrix.rows();
    int columns = complexMatrix.cols();
    MatrixXd realMatrix(rows, columns);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (complexMatrix(i,j).imag() != 0) {
                std::cout << "Warning: non-real values detected!" << std::endl;
            }

            realMatrix(i, j) = complexMatrix(i, j).real();

        }
    }

    return realMatrix;

}

MatrixXd hamiltonianConstructor(int size, double potential) {
    //Constructs the matrix from equation (28) in the manual
    MatrixXd matrix(size, size); 

    for (int n = 0; n < size; n++) {
        for (int m = 0; m < n; m++) {
            if (n == m+1) {
                matrix(n,m) = potential;
                matrix(m, n) = potential;
            }
            else {
                matrix(n, m) = 0;
                matrix(m, n) = 0;
            }
        }
    }
    
    return matrix;
}

std::complex<double> coefficients(int number, double time, 
                                EigenSolver<MatrixXd>& energyBasis, VectorXd& initialValues) {
    const double hbar = 1;
    const std::complex<double> i {0, 1}; //Imaginary unit
    std::complex<double> coefficient = {0, 0}; //Value to be returned
    const int problemSize = initialValues.size();

    //Eigenenergies
    const VectorXd eigenEnergies = makeReal(energyBasis.eigenvalues());

    //Gets the basis vectors and normalises them
    const MatrixXd basis = makeReal(energyBasis.eigenvectors());
    for (int n = 0; n < problemSize; n++) {
        basis.col(n).normalize();
    }

    //Calculates the double sum
    double internalSum = 0;
    double projection = 0;
    double secondProjection = 0;
    for (int lambdaNum = 0; lambdaNum < problemSize; lambdaNum++) {
        internalSum = 0;
        for (int m = 0; m < problemSize; m++) {
            projection = basis.col(lambdaNum)(m);   //c_m^lambdaNum  (m:th value in the lambda:th eigenvector)
            internalSum += projection * initialValues(m) ;      //c_m^lambdaNum * c_m(0)

        }

        std::complex<double> exponential = std::exp(-i* eigenEnergies(lambdaNum)*time/hbar);
        secondProjection = basis.col(lambdaNum)(number);

        coefficient += exponential * secondProjection * internalSum;

    }

    return coefficient;

}

void makePlotInteractive(EigenSolver<MatrixXd>& hamiltonian, VectorXd& initialValues,
            double timeInterval, double timeStepLength, int argc, char** argv) {

    int length = initialValues.size();
    const char* title = "Values of the coefficients at t = ";
    std::string stringTitle = "Values of the coefficients at t = ";

    //ROOT handles removing things from the heap for us
    TRint* application = new TRint("application", &argc, argv);
    TCanvas* canvas = new TCanvas("canvas", "Simulation Result", 0, 0, 800, 600);
    TH1D* hist = new TH1D("hist", title, length, 1, length);
    hist->GetXaxis()->SetTitle("Well Number");
    hist->GetYaxis()->SetTitle("#rho = |c_{n}|^{2}");
    hist->GetXaxis()->CenterTitle(true);
    hist->GetYaxis()->CenterTitle(true);
    hist->GetXaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->SetMaximum(0.5);
    hist->SetStats(0);
    hist->SetFillColor(kBlue-5);

    canvas->Show();

    double time = 0;
    const int million = 1000000;
    std::complex<double> coefficient;
    while (time <= timeInterval) {
        hist->Reset();
        for (int k = 1; k <= length; k++) { //Start at 1 to not go into underflow bin
            coefficient = coefficients(k-1, time, hamiltonian, initialValues);
            double magnitude = pow(std::abs(coefficient),2);
            hist->SetBinContent(k, magnitude);

        }


        std::string timeString = stringTitle + std::to_string(time);
        const char* newTitle = timeString.c_str();
        hist->SetTitle(newTitle);

        hist->Draw();
        canvas->Update();

        usleep(timeStepLength*million); //usleep() works in microseconds
        time += timeStepLength;
    }


    
    
    application->Run();


}


void makePlot(EigenSolver<MatrixXd>& hamiltonian, VectorXd& initialValues,
            double timeInterval, double timeStepLength) {

    const int length = initialValues.size();
    const char* title = "Values of the coefficients at t = ";
    const std::string stringTitle = "Values of the coefficients at t = ";

    //ROOT handles removing things from the heap for us
    TCanvas* canvas = new TCanvas("canvas", "Simulation Result", 0, 0, 800, 600);
    TH1D* hist = new TH1D("hist", title, length, 1, length);
    hist->GetXaxis()->SetTitle("Well Number");
    hist->GetYaxis()->SetTitle("#rho = |c_{n}|^{2}");
    hist->GetXaxis()->CenterTitle(true);
    hist->GetYaxis()->CenterTitle(true);
    hist->GetXaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->SetMaximum(0.5);
    hist->SetStats(0);
    hist->SetFillColor(kBlue-5);

    canvas->Show();

    double time = 0;
    const int million = 1000000;
    std::complex<double> coefficient;
    while (time <= timeInterval) {
        hist->Reset();
        for (int k = 1; k <= length; k++) { //Start at 1 to not go into underflow bin
            coefficient = coefficients(k-1, time, hamiltonian, initialValues);
            double magnitude = pow(std::abs(coefficient),2);
            hist->SetBinContent(k, magnitude);

        }


        std::string timeString = stringTitle + std::to_string(time);
        const char* newTitle = timeString.c_str();
        hist->SetTitle(newTitle);

        hist->Draw();
        canvas->Update();
        canvas->Print("TaskAPlus2.gif+");

        usleep(timeStepLength*million); //usleep() works in microseconds
        time += timeStepLength;
    }



}



int main(int argc, char** argv) {
    //Parameters for the simulation
    const int length = 6;
    const double potential = -1;
    const double timeInterval = 20;
    const double timeStepLength = 0.01;
    const double perturbation = 2;



    //Constructs the Hamiltonians.
    MatrixXd rawHamiltonian = hamiltonianConstructor(length, potential);     //Equation (28) in the manual
    MatrixXd plusHamiltonian = rawHamiltonian; 
    MatrixXd minusHamiltonian = rawHamiltonian;
    plusHamiltonian(0,0) += perturbation;                                    //Equation (29) in the manual
    minusHamiltonian(0,0) -= perturbation;

    //Introduces necessary things to calculate eigenvectors and eigenvalues.
    const EigenSolver<MatrixXd> raw(rawHamiltonian);
    const EigenSolver<MatrixXd> plus(plusHamiltonian);
    const EigenSolver<MatrixXd> minus(minusHamiltonian);

    //Finds the ground state values c_m(0) before the heavi-side function kicks in.
    const VectorXd rawEigenVals = makeReal(rawHamiltonian.eigenvalues());
    auto minValue = std::min_element(rawEigenVals.begin(), rawEigenVals.end());
    int minIndex = std::distance(rawEigenVals.begin(), minValue); 

    VectorXd initialValuesTemp = makeReal(raw.eigenvectors().col(minIndex)); //I am aware of /= notation, but the documentation says to do it this way.
    double norm = initialValuesTemp.squaredNorm();
    const VectorXd initialValues = initialValuesTemp/norm; //This gives the values c_m(0) in the \phi_n basis.

    makePlotInteractive(plus, initialValues, timeInterval, timeStepLength, argc, argv);

}

