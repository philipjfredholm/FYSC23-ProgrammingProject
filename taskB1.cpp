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



//Parameters for the simulation
const int length = 6;
const double potential = -1;
const double timeInterval = 20;
const double timeStepLength = 0.05;
const double perturbation = -15;
const double delta = 10;
const double uValue = 15;



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
                std::cout << complexMatrix(i,j).imag() << std::endl;
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



std::complex<double> coefficients2D(int number1, int number2, double time,
                                    const MatrixXd& basis,
                                    const VectorXd& eigenEnergies, 
                                    const VectorXd& initialValues) {

                    
    //This is just a wrapper for the function coefficients() to
    //not have to calculate the correct indices every time.

    return coefficients(number1*length + number2, time, basis, eigenEnergies, initialValues);
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



/*

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
    hist.SetFillColor(kRed-5);

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



void makePlot(const MatrixXd& hamiltonian, const VectorXd& energies, VectorXd initialValues,
            double timeInterval, double timeStepLength) {

    gSystem->Unlink("TaskAPlus2.gif"); //Stops us from appending to an old gif. 

    const char* title = "Values of the coefficients at t = 0, #varepsilon_1 = ";
    std::string stringTitle = "Values of the coefficients at t = "; 

    TCanvas canvas("canvas", "Simulation Result", 0, 0, 800, 600);
    TH1D hist ("hist", title, length, 1, length);
    hist.GetXaxis()->SetTitle("Well Number");
    hist.GetYaxis()->SetTitle("#rho = |c_{n}|^{2}");
    hist.GetXaxis()->CenterTitle(true);
    hist.GetYaxis()->CenterTitle(true);
    hist.GetXaxis()->SetTitleSize(0.04);
    hist.GetYaxis()->SetTitleSize(0.04);
    hist.SetMaximum(0.5);
    hist.SetStats(0);
    hist.SetFillColor(kRed-5);
    gStyle->SetCanvasPreferGL(kTRUE);

    canvas.Show();

    double time = 0;
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
        canvas.Print("TaskA1.gif+5");

        time += timeStepLength;
    }



}


void plotGraph(const MatrixXd& hamiltonian, const VectorXd& energies, VectorXd initialValues,
            double timeInterval, double timeStepLength, int number) {


    std::string stringTitle = "Values of the coefficient c_" + std::to_string(number) +
                            " as a function of time for #varepsilon_1 = " +
                             std::to_string(perturbation); 
    const char* title = stringTitle.c_str();

    std::string ylabelString = "#rho = |c_{"+std::to_string(number) +"}|^{2}";
    const char* ylabel = ylabelString.c_str();

    TCanvas canvas("canvas", "Simulation Result", 0, 0, 900, 600);
    TGraph* graph = new TGraph();
    graph->SetTitle(title);
    graph->GetXaxis()->SetTitle("Time");
    graph->GetYaxis()->SetTitle(ylabel);
    graph->GetXaxis()->CenterTitle(true);
    graph->GetYaxis()->CenterTitle(true);
    graph->GetXaxis()->SetTitleSize(0.04);
    graph->GetYaxis()->SetTitleSize(0.04);
    gStyle->SetTitleSize(3);
    graph->SetMaximum(0.5);
    graph->SetStats(0);
    graph->SetFillColor(kRed-5);
    gStyle->SetCanvasPreferGL(kTRUE);
    graph->SetLineColor(kRed);
    graph->SetMarkerColor(kRed);

    canvas.Show();

    double time = 0;
    std::complex<double> coefficient;
    while (time <= timeInterval) {
        coefficient = coefficients(number-1, time, hamiltonian, energies, initialValues);
        double magnitude = pow(std::abs(coefficient),2);
        graph->SetPoint(graph->GetN(), time, magnitude);

        

        

        time += timeStepLength;
    }
    graph->Draw("AL");

    
    std::string filenameString ="TaskA1coefficient"+std::to_string(number)+".pdf"; 
    const char* filename = filenameString.c_str();
    gPad->Print(filename);
    
    delete graph;

 }



*/
int main(int argc, char** argv) {
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
    
    std::cout << hamiltonian << std::endl;

    EigenSolver<MatrixXd> eigenData(hamiltonian);
    EigenSolver<MatrixXd> rawEigenData(rawHamiltonian);

    //Finds the ground state
    const VectorXd rawEigenvalues = makeReal(rawEigenData.eigenvalues());
    const auto minValue = std::min_element(rawEigenvalues.begin(), rawEigenvalues.end());
    const int minIndex = std::distance(rawEigenvalues.begin(), minValue); 

    VectorXd initialValuesTemp = makeReal(rawEigenData.eigenvectors().col(minIndex));
    const double norm = initialValuesTemp.squaredNorm();
    const VectorXd initialValues = initialValuesTemp/norm;

    //Finds the basis and eigenenergies after the perturbation is introduced.
    VectorXd energies = makeReal(eigenData.eigenvalues());
    MatrixXd basis = makeReal(eigenData.eigenvectors());
    for (int n = 0; n < std::pow(length,2) ; n++) {
        basis.col(n).normalize();  //In place normalisation
    }

    std::cout << "Plotting starts now" << std::endl;
    interactiveDoubleOccupancy(basis, energies, initialValues, timeInterval, timeStepLength, argc, argv);



/*


    //Makes the plots

    //makePlotInteractive(basis, energies, initialValues, timeInterval, timeStepLength, argc, argv);
    //makePlot(basis, energies, initialValues, timeInterval, timeStepLength);
    for (int n = 1; n <= length; n++) {
        plotGraph(basis, energies, initialValues, timeInterval, timeStepLength, n);

    }

 */
    (void)argc; //Just to turn of the -Werror flag in g++
    (void)argv;


   
}

