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



MatrixXd make2DHamiltonian(int size, double potential, const VectorXd& epsilonValues, const VectorXd& uValues) {
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


            if (n == nPrime && m == mPrime && n == m) {
                core(j, i) += uValues(n);
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

double coefficientsDouble(int number1, double time,
                                    const MatrixXcd& basis,
                                    const VectorXd& eigenEnergies, 
                                    const VectorXd& initialValues) {

                    
    int length = std::sqrt(basis.cols());
    std::complex<double> value = coefficients(number1*length + number1, time, basis, eigenEnergies, initialValues);
    double returnValue = pow(std::abs(value),2);
    


    return returnValue;
}


void plotDoubleGraph(const MatrixXcd& hamiltonian, const VectorXd& energies, VectorXd initialValues,
            double timeInterval, double timeStepLength, double delta) {



    int length = std::sqrt(hamiltonian.cols());
    std::string stringTitle = "Shows the probability of a spin up electron being at a specific site as a function of time for #Delta = " +
                             std::to_string(delta).substr(0,5); 
    const char* title = stringTitle.c_str();

    std::string ylabelString = "Probability";
    const char* ylabel = ylabelString.c_str();

    TCanvas* canvas = new TCanvas("canvas", "Simulation Result", 0, 0, 900, 600);
    TMultiGraph* graph = new TMultiGraph("test", "test");
    TLegend* myLegend = new TLegend(0.7, 0.5, .9, .7);
    
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
    std::string filenameString ="TaskB3coefficients.pdf"; 
    const char* filename = filenameString.c_str();
    gPad->Print(filename);
        

    

 }




void plotDoubleGraph2(const MatrixXcd& hamiltonian, const VectorXd& energies, VectorXd initialValues,
            double timeInterval, double timeStepLength, double delta) {



    int length = std::sqrt(hamiltonian.cols());
    std::string stringTitle = "Shows the double occupancy probability as a function of time for #Delta = " +
                             std::to_string(delta).substr(0,5); 
    const char* title = stringTitle.c_str();

    std::string ylabelString = "Probability";
    const char* ylabel = ylabelString.c_str();

    TCanvas* canvas = new TCanvas("canvas", "Simulation Result", 0, 0, 900, 600);
    TMultiGraph* graph = new TMultiGraph("test", "test");
    TLegend* myLegend = new TLegend(0.7, 0.5, .9, .7);
    
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
            coefficient = coefficientsDouble(number, time, hamiltonian, energies, initialValues);
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
    std::string filenameString ="TaskB4coefficients.pdf"; 
    const char* filename = filenameString.c_str();
    gPad->Print(filename);
        

    

 }



void plotDoubleGraph3(const MatrixXcd& hamiltonian, const VectorXd& energies, VectorXd initialValues,
            double timeInterval, double timeStepLength, double delta) {



    int length = std::sqrt(hamiltonian.cols());
    std::string stringTitle = "Shows |#rho_{n #uparrow } - #rho_{n}^{(2)}| for #Delta = " +
                             std::to_string(delta).substr(0,5); 
    const char* title = stringTitle.c_str();

    std::string ylabelString = "Probability";
    const char* ylabel = ylabelString.c_str();

    TCanvas* canvas = new TCanvas("canvas", "Simulation Result", 0, 0, 900, 600);
    TMultiGraph* graph = new TMultiGraph("test", "test");
    TLegend* myLegend = new TLegend(0.7, 0.5, .9, .7);
    
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
            double coefficient1 = coefficientsDouble(number, time, hamiltonian, energies, initialValues);
            double coefficient2 = coefficientsSingle(number, time, hamiltonian, energies, initialValues);
            coefficient = std::abs(coefficient1 - coefficient2);
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
    std::string filenameString ="TaskB5diffDelta" + std::to_string(delta).substr(0,4) + ".pdf"; 
    const char* filename = filenameString.c_str();
    gPad->Print(filename);
        

    

 }




void singleOccupancy(int length, double V, double initialEpsilon, double delta, double uValue) {


  //Sets up the initial values
    VectorXd initialEpsilons(length);
    initialEpsilons.setZero();
    //VectorXd zeroInitialEnergy = epsilons; //For debugging purposes
    VectorXd uValues = initialEpsilons;

    uValues(0) = 0;
    for (int n = 1; n < length; n++) {
        uValues(n) = uValue;
    }

    initialEpsilons(0) = initialEpsilon;
    VectorXd epsilons = initialEpsilons;
    epsilons(0) += delta;

    /*
    std::cout << "initialEpsilons:" << initialEpsilons << std::endl;
    std::cout << "uValues: " << uValues << std::endl;
    std::cout << "epsilons: " << epsilons << std::endl;
    */
    

    //Computes initial values for the simulation
    MatrixXd initialHamiltonian = make2DHamiltonian(length, V, initialEpsilons, uValues);
    SelfAdjointEigenSolver<MatrixXd> initialEigenSolver(initialHamiltonian);
    VectorXd initialEigenValues = makeReal(initialEigenSolver.eigenvalues());
    MatrixXcd initialEnergyBasis = initialEigenSolver.eigenvectors();
    

    const auto minValue = std::min_element(initialEigenValues.begin(), initialEigenValues.end());
    const int minIndex = std::distance(initialEigenValues.begin(), minValue);
    const VectorXd initialValues = makeReal(initialEnergyBasis.col(minIndex));



    //Computes the eigenenergies and eigenvectors to be used in the simulation
    MatrixXd hamiltonian = make2DHamiltonian(length, V, epsilons, uValues);
    SelfAdjointEigenSolver<MatrixXd> eigenSolver(hamiltonian);
    VectorXd eigenValues = makeReal(eigenSolver.eigenvalues());
    MatrixXcd energyBasis = eigenSolver.eigenvectors();


    plotDoubleGraph(energyBasis, eigenValues, initialValues, 20, 0.001, delta);


}




void doubleOccupancy(int length, double V, double initialEpsilon, double delta, double uValue) {


  //Sets up the initial values
    VectorXd initialEpsilons(length);
    initialEpsilons.setZero();
    //VectorXd zeroInitialEnergy = epsilons; //For debugging purposes
    VectorXd uValues = initialEpsilons;

    uValues(0) = 0;
    for (int n = 1; n < length; n++) {
        uValues(n) = uValue;
    }

    initialEpsilons(0) = initialEpsilon;
    VectorXd epsilons = initialEpsilons;
    epsilons(0) += delta;

    /*
    std::cout << "initialEpsilons:" << initialEpsilons << std::endl;
    std::cout << "uValues: " << uValues << std::endl;
    std::cout << "epsilons: " << epsilons << std::endl;
    */
    

    //Computes initial values for the simulation
    MatrixXd initialHamiltonian = make2DHamiltonian(length, V, initialEpsilons, uValues);
    SelfAdjointEigenSolver<MatrixXd> initialEigenSolver(initialHamiltonian);
    VectorXd initialEigenValues = makeReal(initialEigenSolver.eigenvalues());
    MatrixXcd initialEnergyBasis = initialEigenSolver.eigenvectors();
    

    const auto minValue = std::min_element(initialEigenValues.begin(), initialEigenValues.end());
    const int minIndex = std::distance(initialEigenValues.begin(), minValue);
    const VectorXd initialValues = makeReal(initialEnergyBasis.col(minIndex));



    //Computes the eigenenergies and eigenvectors to be used in the simulation
    MatrixXd hamiltonian = make2DHamiltonian(length, V, epsilons, uValues);
    SelfAdjointEigenSolver<MatrixXd> eigenSolver(hamiltonian);
    VectorXd eigenValues = makeReal(eigenSolver.eigenvalues());
    MatrixXcd energyBasis = eigenSolver.eigenvectors();


    plotDoubleGraph2(energyBasis, eigenValues, initialValues, 20, 0.001, delta);


}


void occupancyDiff(int length, double V, double initialEpsilon, double delta, double uValue) {


  //Sets up the initial values
    VectorXd initialEpsilons(length);
    initialEpsilons.setZero();
    //VectorXd zeroInitialEnergy = epsilons; //For debugging purposes
    VectorXd uValues = initialEpsilons;

    uValues(0) = 0;
    for (int n = 1; n < length; n++) {
        uValues(n) = uValue;
    }

    initialEpsilons(0) = initialEpsilon;
    VectorXd epsilons = initialEpsilons;
    epsilons(0) += delta;

    /*
    std::cout << "initialEpsilons:" << initialEpsilons << std::endl;
    std::cout << "uValues: " << uValues << std::endl;
    std::cout << "epsilons: " << epsilons << std::endl;
    */
    

    //Computes initial values for the simulation
    MatrixXd initialHamiltonian = make2DHamiltonian(length, V, initialEpsilons, uValues);
    SelfAdjointEigenSolver<MatrixXd> initialEigenSolver(initialHamiltonian);
    VectorXd initialEigenValues = makeReal(initialEigenSolver.eigenvalues());
    MatrixXcd initialEnergyBasis = initialEigenSolver.eigenvectors();
    

    const auto minValue = std::min_element(initialEigenValues.begin(), initialEigenValues.end());
    const int minIndex = std::distance(initialEigenValues.begin(), minValue);
    const VectorXd initialValues = makeReal(initialEnergyBasis.col(minIndex));



    //Computes the eigenenergies and eigenvectors to be used in the simulation
    MatrixXd hamiltonian = make2DHamiltonian(length, V, epsilons, uValues);
    SelfAdjointEigenSolver<MatrixXd> eigenSolver(hamiltonian);
    VectorXd eigenValues = makeReal(eigenSolver.eigenvalues());
    MatrixXcd energyBasis = eigenSolver.eigenvectors();


    plotDoubleGraph3(energyBasis, eigenValues, initialValues, 20, 0.001, delta);


}




int main() {
    //Parameters
    const int length = 6;
    const double V = -1;
    const double initialEpsilon = -15;
    //const double delta = 15;
    const double uValue = 15;


    //singleOccupancy(length, V, initialEpsilon, delta, uValue);
    //doubleOccupancy(length, V, initialEpsilon, delta, uValue);

    
    double deltaStart = 20;

    for (double n = deltaStart; n <= 25; n+=0.5) {
        std::cout << n << std::endl;
        occupancyDiff(length, V, initialEpsilon, n, uValue);
    }
    


    



}