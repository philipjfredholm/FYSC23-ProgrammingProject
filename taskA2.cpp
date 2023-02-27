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
#include <TSystem.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>


//Parameters for the simulation
const int length = 6;
const double potential = -1;
const double timeInterval = 20;
const double timeStepLength = 0.05;
const double perturbation = -2;


MatrixXd makeReal(const MatrixXcd& complexMatrix) {
    const int rows = complexMatrix.rows();
    const int columns = complexMatrix.cols();
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
                                const MatrixXd& basis,
                                const VectorXd& eigenEnergies,
                                const VectorXd& initialValues) {
    const double hbar = 1;
    const std::complex<double> i {0, 1}; //Imaginary unit
    std::complex<double> coefficient = {0, 0}; //Value to be returned

    //Calculates the double sum
    double internalSum = 0;
    double projection = 0;
    double secondProjection = 0;
    for (int lambdaNum = 0; lambdaNum < length; lambdaNum++) {
        internalSum = 0;
        for (int m = 0; m < length; m++) {
            projection = basis.col(lambdaNum)(m);   //c_m^lambdaNum  (m:th value in the lambda:th eigenvector)
            internalSum += projection * initialValues(m) ;      //c_m^lambdaNum * c_m(0)

        }

        std::complex<double> exponential = std::exp(-i* eigenEnergies(lambdaNum)*time/hbar);
        secondProjection = basis.col(lambdaNum)(number);

        coefficient += exponential * secondProjection * internalSum;

    }

    return coefficient;

}

void makePlotInteractive(const MatrixXd& hamiltonian, const VectorXd& energies, const VectorXd& initialValues,
            double timeInterval, double timeStepLength, int argc, char** argv) {

    const char* title = "Values of the coefficients at t = , #varepsilon_1 =";
    std::string stringTitle = "Values of the coefficients at t = ";

    TRint application ("application", &argc, argv);
    TCanvas canvas ("canvas", "Simulation Result", 0, 0, 800, 600);
    TH1D hist ("hist", title, length, 1, length);
    hist.GetXaxis()->SetTitle("Well Number");
    hist.GetYaxis()->SetTitle("#rho = |c_{n}|^{2}");
    hist.GetXaxis()->CenterTitle(true);
    hist.GetYaxis()->CenterTitle(true);
    hist.GetXaxis()->SetTitleSize(0.04);
    hist.GetYaxis()->SetTitleSize(0.045);
    hist.SetMaximum(0.5);
    hist.SetStats(0);
    hist.SetFillColor(kBlue-5);

    canvas.Show();

    double time = 0;
    int million = 1000000;
    std::complex<double> coefficient;
    while (time <= timeInterval) {
        hist.Reset();
        for (int k = 1; k <= length; k++) { //Start at 1 to not go into underflow bin
            coefficient = coefficients(k-1, time, hamiltonian, energies, initialValues);
            double magnitude = pow(std::abs(coefficient),2);
            hist.SetBinContent(k, magnitude);

        }


        std::string timeString = stringTitle + std::to_string(time)
                                +" , #varepsilon_1 = " + std::to_string(perturbation);
        const char* newTitle = timeString.c_str();
        hist.SetTitle(newTitle);

        hist.Draw();
        canvas.Update();

        usleep(timeStepLength*million); //usleep() works in microseconds
        time += timeStepLength;
    }


    
    
    application.Run();


}


void plotGraph(const MatrixXd& hamiltonian, const VectorXd& energies, VectorXd initialValues,
            double timeInterval, double timeStepLength) {


    std::string stringTitle = "Values of the coefficients #rho_{n} = |c_{n}|^{2} as a function of time for #varepsilon_{1} = " +
                             std::to_string(perturbation).substr(0,5); 
    const char* title = stringTitle.c_str();

    std::string ylabelString = "#rho = |c_{n}|^{2}";
    const char* ylabel = ylabelString.c_str();

    TCanvas* canvas = new TCanvas("canvas", "Simulation Result", 0, 0, 900, 600);
    TMultiGraph* graph = new TMultiGraph("test", "test");
    TLegend* myLegend = new TLegend(0.7, 0.7, .9, .9);
    
    graph->SetTitle(title);
    graph->GetXaxis()->SetTitle("Time (seconds #upoint #hbar)");
    graph->GetYaxis()->SetTitle(ylabel);
    graph->GetXaxis()->CenterTitle(true);
    graph->GetYaxis()->CenterTitle(true);
    graph->GetXaxis()->SetTitleSize(0.04);
    graph->GetYaxis()->SetTitleSize(0.04);
    gStyle->SetTitleSize(3);
    graph->SetMaximum(1);
    gPad->SetGrid();
    graph->GetXaxis()->SetNdivisions(20);
    graph->GetXaxis()->SetLabelOffset(0.01);
    graph->GetYaxis()->SetNdivisions(20);
    //graph->GetYaxis()->SetLabelOffset(0.01);
    graph->GetXaxis()->SetLimits(0,20);
    //graph->SetStats(0);
    gStyle->SetCanvasPreferGL(kTRUE);
    std::string legendNumber = "n = ";
    canvas->Show();

    double time = 0;

 

    for (int number = 0; number < length; number++) {
        

        std::string stringNumber = std::to_string(number+1);
        const char* legendName = (legendNumber+stringNumber).c_str();


        TGraph* newGraph = new TGraph();
        graph->SetTitle(legendName);
        graph->SetName(legendName);

        time = 0;
        std::complex<double> coefficient;
        while (time <= timeInterval) {
            coefficient = coefficients(number, time, hamiltonian, energies, initialValues);
            double magnitude = pow(std::abs(coefficient),2);
            newGraph->SetPoint(newGraph->GetN(), time, magnitude);

            

            

            time += timeStepLength;
        }
  
        myLegend->AddEntry(newGraph, legendName, "L");
        graph->Add(newGraph, "AL");

        
    }
     
    
    graph->Draw("AL PLC");
    myLegend->SetNColumns(2);
    myLegend->Draw();
    std::string filenameString ="TaskA2coefficients.pdf"; 
    const char* filename = filenameString.c_str();
    gPad->Print(filename);
        

    

 }





int main(int argc, char** argv) {
    //Constructs the Hamiltonians.
    MatrixXd rawHamiltonian = hamiltonianConstructor(length, potential);     //Equation (28) in the manual
    MatrixXd Hamiltonian = rawHamiltonian; 
    Hamiltonian(0,0) += perturbation;                                    //Equation (29) in the manual

    //Introduces necessary things to calculate eigenvectors and eigenvalues.
    EigenSolver<MatrixXd> raw(rawHamiltonian);
    EigenSolver<MatrixXd> final(Hamiltonian);

    //Finds the ground state values c_m(0) before the heavi-side function kicks in.
    const VectorXd rawEigenVals = makeReal(rawHamiltonian.eigenvalues());
    const auto minValue = std::min_element(rawEigenVals.begin(), rawEigenVals.end());
    const int minIndex = std::distance(rawEigenVals.begin(), minValue); 

    VectorXd initialValuesTemp = makeReal(raw.eigenvectors().col(minIndex)); //I am aware of /= notation, but the documentation says to do it this way.
    double norm = initialValuesTemp.squaredNorm();
    const VectorXd initialValues = initialValuesTemp/norm; //This gives the values c_m(0) in the \phi_n basis.


    //Finds the values for <\phi_n | \varphi_\lambda> and eigenvalues E_\lambda
    const VectorXd energies = makeReal(final.eigenvalues());

    MatrixXd basis= makeReal(final.eigenvectors());
    //Normalises the eigenvectors
    for (int n = 0; n < length; n++) {
        basis.col(n).normalize(); 
    }



    //Makes the plots

    //makePlotInteractive(basis, energies, initialValues, timeInterval, timeStepLength, argc, argv);
    //makePlot(basis, energies, initialValues, timeInterval, timeStepLength);
    
    

    plotGraph(basis, energies, initialValues, timeInterval, timeStepLength);
    
    (void)argc; //Just to turn of the -Werror flag in g++
    (void)argv;
}

