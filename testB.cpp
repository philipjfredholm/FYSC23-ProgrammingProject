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


MatrixXd makeReal(const MatrixXcd& complexMatrix) {
    int columns = complexMatrix.cols();
    int rows = complexMatrix.rows();
    MatrixXd returnMatrix(rows, columns);
    double threshold = std::pow(10,-8);

    for (int n = 0; n < columns; n++) {
        for (int k = 0; k < rows; k++) {
            if (std::abs(complexMatrix(k, n).imag()) > threshold) {
                std::cout << "Non-real values detected, i: " << complexMatrix(k,n).imag() << std::endl;
            }

            returnMatrix(k, n) = complexMatrix(k, n).real();

        }
    }

    return returnMatrix;

}


MatrixXd make1DHamiltonian(int size, double potential) {
    MatrixXd core(size, size);
    core.setZero();

    for (int n = 0; n < size; n++) { //Columns
        for (int k = 0; k < size; k++) { //Rows
            if (n == k+1) {
                core(k, n) = potential;

            }
            if (k == n+1) {
                core(k,n) = potential;
            }
        }
    }

    return core;
}



MatrixXd make2DHamiltonian(int size, double potential, const VectorXd& epsilonValues) {
    int squaredSize = std::pow(size, 2);
    MatrixXd core(squaredSize, squaredSize);
    core.setZero();

    MatrixXd oneDimensionEquivalent = make1DHamiltonian(size, potential);

    for (int i = 0; i < squaredSize; i++) { //Column
        int n = std::trunc(i/size);
        int nPrime = i%size;
        for (int j = 0; j < squaredSize; j++) { //Row
            int m = std::trunc(j/size);
            int mPrime = j%size;

            if (nPrime == mPrime) {
                core(j,i) += oneDimensionEquivalent(m, n); 
            }

            if (n == m) {
                core(j, i) += oneDimensionEquivalent(mPrime, nPrime);
            }


            if (n == m && nPrime == mPrime) {
                core(j, i) += epsilonValues(n) + epsilonValues(nPrime);
            }

        }
    }

    return core;

}


std::complex<double> coefficients(int number, double time, const MatrixXcd& basis, const VectorXcd& energies, const VectorXcd& initialValues) {
    int numberOfVectors = basis.cols();
    std::complex<double> i {0,1};
    double hbar = 1;

    std::complex<double> internalSum  {0, 0};
    std::complex<double> projectionFactor {0,0};
    std::complex<double> returnValue {0, 0};
    std::complex<double> externalProjectionFactor {0, 0};


    for (int lambdaNum = 0; lambdaNum < numberOfVectors; lambdaNum++) {
        internalSum = {0, 0};
        
        for (int m = 0; m < numberOfVectors; m++) {
            projectionFactor = {0, 0};
            projectionFactor = basis(m, lambdaNum);
            internalSum += projectionFactor*initialValues(m);
        }
        
        std::complex<double> exponentialFactor = exp(-i*energies(lambdaNum)*time/hbar);
        externalProjectionFactor = {0, 0};
        externalProjectionFactor = basis(number, lambdaNum);

        returnValue += exponentialFactor*externalProjectionFactor*internalSum;


    }

    return returnValue;


}

double coefficientsSingle(int number1, double time,
                                    const MatrixXcd& basis,
                                    const VectorXd& eigenEnergies, 
                                    const VectorXd& initialValues) {

                    
    double internalSum = 0;
    int length = std::sqrt(basis.cols());
    //int newLength = basis.cols();
    for (int n = 0; n < length; n++) {
        std::complex<double> value = coefficients(number1*length + n,
                         time, basis, eigenEnergies, initialValues);
        double magnitude = pow(std::abs(value),2);
        
        internalSum += magnitude;
    }

    return internalSum;
}




void plotDoubleGraph(const MatrixXcd& hamiltonian, const VectorXd& energies, VectorXd initialValues,
            double timeInterval, double timeStepLength) {



    int length = std::sqrt(hamiltonian.cols());
    std::string stringTitle = "Values of the coefficients #rho_{n} = |c_{n}|^{2} as a function of time for #varepsilon_{1} = " +
                             std::to_string(0).substr(0,5); 
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
        double coefficient;
        while (time <= timeInterval) {
            coefficient = coefficientsSingle(number, time, hamiltonian, energies, initialValues);
            //std::cout << coefficient << std::endl;
            newGraph->SetPoint(newGraph->GetN(), time, coefficient);

            

            

            time += timeStepLength;
        }
  
        myLegend->AddEntry(newGraph, legendName, "L");
        graph->Add(newGraph, "AL");

        
    }
     
    
    graph->Draw("AL PLC");
    myLegend->SetNColumns(2);
    myLegend->Draw();
    std::string filenameString ="TaskB1coefficients.pdf"; 
    const char* filename = filenameString.c_str();
    gPad->Print(filename);
        

    

 }




int main() {
    //Parameters
    const int length = 3;
    const double V = -1;
    const double initialEpsilon = -2;

    VectorXd epsilons(length);
    epsilons.setZero();
    VectorXd zeroInitialEnergy = epsilons;
    epsilons(0) = initialEpsilon;

    //Initial Values
    MatrixXd initialHamiltonian = make2DHamiltonian(length, V, zeroInitialEnergy);
    EigenSolver<MatrixXd> initialEigenSolver(initialHamiltonian);
    VectorXd initialEigenValues = makeReal(initialEigenSolver.eigenvalues());
    std::cout << initialEigenValues << std::endl;
    MatrixXcd initialEnergyBasis = initialEigenSolver.eigenvectors();

    const auto minValue = std::min_element(initialEigenValues.begin(), initialEigenValues.end());
    const int minIndex = std::distance(initialEigenValues.begin(), minValue);
    const VectorXd initialValues = makeReal(initialEnergyBasis.col(minIndex));



    //Simulation values
    MatrixXd hamiltonian = make2DHamiltonian(length, V, epsilons);
    EigenSolver<MatrixXd> eigenSolver(hamiltonian);
    VectorXd eigenValues = makeReal(eigenSolver.eigenvalues());
    MatrixXcd energyBasis = eigenSolver.eigenvectors();


    //Comparing if values match
    std::complex<double> coefficent = coefficients(1,0, energyBasis, eigenValues, initialValues);

    std::cout << "Calculated value " << coefficent << std::endl;
    std::cout << "Input value " << initialValues(1) << std::endl;

    plotDoubleGraph(energyBasis, eigenValues, initialValues, 20, 0.001);


    /*
    std::cout << "The ground state is " << std::endl;
    VectorXcd groundState = initialEnergyBasis.col(0);
    std::cout << groundState << std::endl;

    std::complex<double> singleOccupancy = 0;
    std::complex<double> temp = 0;
    for (int n = 0; n < 3; n++) {
        temp = std::pow( groundState(n)   , 2);
        singleOccupancy += temp;
    }


    std::cout << "Single occpancy probability" << std::endl;
    std::cout << singleOccupancy << std::endl;

    */




}