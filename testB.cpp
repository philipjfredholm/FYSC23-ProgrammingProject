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


int main() {
    const int length = 3;
    const double V = -1;
    const double initialEpsilon = -2;

    VectorXd epsilons(length);
    epsilons.setZero();
    VectorXd zeroInitialEnergy = epsilons;
    epsilons(0) = initialEpsilon;

    MatrixXd initialHamiltonian = make2DHamiltonian(length, V, zeroInitialEnergy);
    EigenSolver<MatrixXd> initialEigenSolver(initialHamiltonian);
    VectorXcd initialEigenValues = initialEigenSolver.eigenvalues();
    std::cout << initialEigenValues << std::endl;
    MatrixXcd initialEnergyBasis = initialEigenSolver.eigenvectors();
    std::cout << initialEnergyBasis << std::endl;



    MatrixXd hamiltonian = make2DHamiltonian(length, V, epsilons);
    EigenSolver<MatrixXd> eigenSolver(hamiltonian);
    VectorXd eigenValues = makeReal(eigenSolver.eigenvalues());
    MatrixXd energyBasis = makeReal(eigenSolver.eigenvectors());


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