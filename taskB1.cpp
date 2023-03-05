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
const double perturbation = 0;
const double delta = -2;
const double uValue = 0;



//For running the simulation
MatrixXd makeReal(const MatrixXcd& complexMatrix) {
    const double threshold = std::pow(10, -8);
    const int rows = complexMatrix.rows();
    const int columns = complexMatrix.cols();
    MatrixXd realMatrix(rows, columns);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (complexMatrix(i,j).imag() > threshold) {
                std::cout << "Warning: non-real values detected!" << std::endl;
                std::cout << "Imaginary Value: " <<  complexMatrix(i,j).imag() << std::endl;
            }

            realMatrix(i, j) = complexMatrix(i, j).real();

        }
    }

    return realMatrix;

}

MatrixXd HamiltonianConstructor(int size, double potential) {
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
        matrix(n,n) = 0;
    }

    
    return matrix;
}



MatrixXd hamiltonianConstructor(int size, double potential) {
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



MatrixXd doubleHamiltonianConstructor(int size, double potential, VectorXd uValues, VectorXd epsilonValues) {
    int newSize = std::pow(size,2);
    MatrixXd matrix(newSize, newSize);
    MatrixXd oneDimensionCase = hamiltonianConstructor(size, potential);
    for (int i = 0; i < newSize; i++) {
        //Each column of the Hamiltonian corresponds to a basis vector
        int basisNumberOne = std::trunc(i/size); //n
        int basisNumberTwo = i%size;                 //n'
        for (int j = 0; j < newSize; j++) {     
            int projectionBasisOne = std::trunc(j/size); // m
            int projectionBasisTwo = j%size ; // m'

            //Raw Structure
            if (basisNumberTwo == projectionBasisTwo) {
                matrix(i, j) = oneDimensionCase(basisNumberOne, projectionBasisOne);
            }
            else {
                matrix(i,j) = 0;
            }

            if (basisNumberOne == projectionBasisOne) {
                matrix(i,j) += oneDimensionCase(basisNumberTwo, projectionBasisTwo);
            }

            //U values
            if (basisNumberTwo == projectionBasisTwo && basisNumberOne == projectionBasisOne &&
                                                basisNumberTwo == projectionBasisOne) {
                matrix(i,j) += uValues(basisNumberOne);
            }

            //Epsilon values
            if (i == j) {
                matrix(i,j) += epsilonValues(basisNumberOne) + epsilonValues(basisNumberTwo);

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

    int newLength = basis.cols();

    //Calculates the double sum
    double internalSum = 0;
    double projection = 0;
    double secondProjection = 0;
    for (int lambdaNum = 0; lambdaNum < newLength; lambdaNum++) {
        internalSum = 0;
        for (int m = 0; m < newLength; m++) {
            projection = basis.col(lambdaNum)(m);   //c_m^lambdaNum  (m:th value in the lambda:th eigenvector)
            internalSum += projection * initialValues(m) ;      //c_m^lambdaNum * c_m(0)

        }

        std::complex<double> exponential = std::exp(-i* eigenEnergies(lambdaNum)*time/hbar);
        secondProjection = basis.col(lambdaNum)(number);

        coefficient += exponential * secondProjection * internalSum;

    }

    return coefficient;

}



std::complex<double> coefficients2D(int number1, int number2, double time,
                                    const MatrixXd& basis,
                                    const VectorXd& eigenEnergies, 
                                    const VectorXd& initialValues) {

                    
    //This is just a wrapper for the function coefficients() to
    //not have to calculate the correct indices every time.

    return coefficients(number1*length + number2, time, basis, eigenEnergies, initialValues);
}




std::complex<double> coefficientsDouble(int number1, double time,
                                    const MatrixXd& basis,
                                    const VectorXd& eigenEnergies, 
                                    const VectorXd& initialValues) {

                    


    return coefficients(number1*length + number1, time, basis, eigenEnergies, initialValues);
}



double coefficientsSingle(int number1, double time,
                                    const MatrixXd& basis,
                                    const VectorXd& eigenEnergies, 
                                    const VectorXd& initialValues) {

                    
    double internalSum = 0;
    //int newLength = basis.cols();
    for (int n = 0; n < length; n++) {
        std::complex<double> value = coefficients(number1*length + n,
                         time, basis, eigenEnergies, initialValues);
        double magnitude = pow(std::abs(value),2);
        
        internalSum += magnitude;
    }

    return internalSum;
}






//Plotting
void interactiveDoubleOccupancy(const MatrixXd& basis, const VectorXd& energies, const VectorXd& initialValues,
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
    hist.SetMaximum(0.005);
    hist.SetStats(0);
    hist.SetFillColor(kBlue-5);

    canvas.Show();

    double time = 0;
    int million = 1000000;
    std::complex<double> coefficient;
    while (time <= timeInterval) {
        hist.Reset();
        for (int k = 1; k <= length; k++) { //Start at 1 to not go into underflow bin
            coefficient = coefficients2D(k-1, k-1, time, basis, energies, initialValues);
            double magnitude = pow(std::abs(coefficient),2);
            if (hist.GetMaximum() < magnitude*1.25) {
                hist.SetMaximum(magnitude*1.25);
            }


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





void plotDoubleGraph(const MatrixXd& hamiltonian, const VectorXd& energies, VectorXd initialValues,
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



int main(int argc, char** argv) {
    std::cout.precision(16);
    //Sets up initial values for \varepsilon and the U-values
    VectorXd uValues(length);
    VectorXd timeZeroValues(length);
    uValues(0) = 0;
    timeZeroValues(0) = perturbation;

    for (int n = 1; n < length; n++) {
        uValues(n) = uValue;
        timeZeroValues(n) = 0;
    }

    VectorXd epsilonValues = timeZeroValues;
    epsilonValues(0) += delta;

    //Prepares matrices and eigenvalues
    const MatrixXd hamiltonian = doubleHamiltonianConstructor(length, potential, uValues, epsilonValues);
    const MatrixXd rawHamiltonian = doubleHamiltonianConstructor(length, potential, uValues, timeZeroValues);
    
    std::cout << "line" << std::endl;
    std::cout << hamiltonian << std::endl;

    EigenSolver<MatrixXd> eigenData(hamiltonian);
    EigenSolver<MatrixXd> rawEigenData(rawHamiltonian);

    //Finds the ground state
    std::cout << "Ground state calculation starts now" << std::endl;
    const VectorXd rawEigenvalues = makeReal(rawEigenData.eigenvalues());
    const auto minValue = std::min_element(rawEigenvalues.begin(), rawEigenvalues.end());
    const int minIndex = std::distance(rawEigenvalues.begin(), minValue); 

    
    std::cout << "mid" << std::endl;
    //std::cout << rawEigenData.eigenvectors().col(minIndex).real() << std::endl;
    VectorXd initialValuesTemp = makeReal(rawEigenData.eigenvectors().col(minIndex));
    const double norm = initialValuesTemp.norm();
    const VectorXd initialValues = initialValuesTemp/norm;


    std::cout << "Basis Calculation starts now" << std::endl;
    //Finds the basis and eigenenergies after the perturbation is introduced.
    VectorXd energies = makeReal(eigenData.eigenvalues());


    std::cout << energies << std::endl;

    MatrixXd basis = makeReal(eigenData.eigenvectors());
    std::cout << basis << std::endl;



    for (int n = 0; n < std::pow(length,2) ; n++) {
        basis.col(n).normalize();  //In place normalisation
    }

    //std::cout << basis << std::endl;
    VectorXd testvector = basis.col(0);
    double temp = testvector.norm();
    VectorXd testvector2 = testvector/temp;
    double myvalue = 0;
    for (int n = 0; n < 6; n++) {
        myvalue += pow(std::abs(testvector2[n]),2);
        
    }


    //std::cout << rawEigenvalues << std::endl;
    std::cout << "Plotting starts now" << std::endl;

    plotDoubleGraph(basis, energies, initialValues, timeInterval, timeStepLength);

    
    std::cout << coefficients(3, 0, basis, energies, initialValues) << "and " << initialValues[3] << std::endl;

    /*int n1 = 4;
    int n2 = 5;
    std::cout << coefficients(n1*length + n2, 5, basis, energies, initialValues) << std::endl;
    std::cout << coefficients(n2*length + n1, 5, basis, energies, initialValues) << std::endl;
    */

    (void)argc; //Just to turn of the -Werror flag in g++
    (void)argv;


   
}

