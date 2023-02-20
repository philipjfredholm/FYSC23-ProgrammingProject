#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <Eigen/Dense>


#include <TH1D.h>


/*
double coefficient(int N, double startValue, int nStates, double orbitalEnergy, double time) {
    std::complex<double> i {0, 1}; //Defines the imaginary unit
    double reducedPlanck = 1; //1.054571817*std::pow(10, -34);

    std::complex<double> value  {0, 0};         //Value to be returned
    std::complex<double> exponential {0, 0}; //Placeholder value
    for (int k = 0; int k > nStates; k++) {
        exponential = std::exp(-i*orbitalEnergy*time/reducedPlanck);
        


    }
}
*/

Eigen::MatrixXd hamiltonianConstructor(int size, double potential) {
    Eigen::MatrixXd matrix(size, size); //Defines the matrix
    
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



double energy(Eigen::VectorXcd energies, Eigen::VectorXd coefficients) {
    double totalEnergy = 0;
    for (int n = 0; n < energies.size(); n++) {
        totalEnergy += real(energies[n]) * abs(coefficients[n]);
        std::cout << real(energies[n]) << std::endl;

    }

    return totalEnergy;
}



int main() {

    //Parameters
    const int length = 6;
    const double potential = -1;
    //const double timeInterval = 20;
    const double perturbation = 2;

    //Constructs the Hamiltonians and finds the eigenvalues (tantamount to diagonalising them)
    Eigen::MatrixXd rawHamiltonian = hamiltonianConstructor(length, potential);
    Eigen::MatrixXd plusHamiltonian = rawHamiltonian;
    Eigen::MatrixXd minusHamiltonian = rawHamiltonian;
    plusHamiltonian(0,0) += perturbation;
    minusHamiltonian(0,0) -= perturbation;


    Eigen::VectorXcd rawEigenVals = rawHamiltonian.eigenvalues(); //Gives an error message if we do not accept complex values (Xcd instead of Xd)
    Eigen::VectorXcd plusEigenVals = plusHamiltonian.eigenvalues();
    Eigen::VectorXcd minusEigenVals = minusHamiltonian.eigenvalues();

    std::cout << rawEigenVals << std::endl;
    std::cout << plusEigenVals << std::endl;
    std::cout << minusEigenVals << std::endl;


    //Eigen::VectorXd initialCoefficients {{1, 2, 1, 9, 4, 5}};
    //std::cout << initialCoefficients << std::endl;
    //std::cout << energy(rawEigenVals, initialCoefficients) << std::endl;
    




}

