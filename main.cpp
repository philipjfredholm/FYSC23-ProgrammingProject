#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace Eigen;

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


std::complex<double> coefficient(int N, ) {
    
}


VectorXd projectCoefficients(VectorXd state);
    /* Projects the coefficients of the energy eigenstates to the
    desired basis. */

    //Normalisation
    VectorXd temp = state; //Necessary according to the documentation
    double length = state.squaredNorm();
    state = temp/length; 

    //No need to project this onto the basis vectors as they are
    //already just [1,0,...,0] , [0,1,0,...0] etc.

    return state;


MatrixXd hamiltonianConstructor(int size, double potential) {
    MatrixXd matrix(size, size); //Defines the matrix
    
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



double energy(VectorXcd energies, VectorXd coefficients) {
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
    MatrixXd rawHamiltonian = hamiltonianConstructor(length, potential);
    MatrixXd plusHamiltonian = rawHamiltonian;
    MatrixXd minusHamiltonian = rawHamiltonian;
    plusHamiltonian(0,0) += perturbation;
    minusHamiltonian(0,0) -= perturbation;

    EigenSolver< MatrixXd> raw(rawHamiltonian);
    EigenSolver<MatrixXd> plus(plusHamiltonian);
    EigenSolver<MatrixXd> minus(minusHamiltonian);

    VectorXcd rawEigenVals = rawHamiltonian.eigenvalues(); //Gives an error message if we do not accept complex values (Xcd instead of Xd)
    VectorXcd plusEigenVals = plusHamiltonian.eigenvalues();
    VectorXcd minusEigenVals = minusHamiltonian.eigenvalues();

    //MatrixXcd rawEigenVecs = raw.eigenvectors();
    //std::cout << rawEigenVecs.col(0) << std::endl;



    std::cout << raw.eigenvectors().col(0) << std::endl;
    //std::cout << plusEigenVals << std::endl;
    //std::cout << minusEigenVals << std::endl;


    //VectorXd initialCoefficients {{1, 2, 1, 9, 4, 5}};
    //std::cout << initialCoefficients << std::endl;
    //std::cout << energy(rawEigenVals, initialCoefficients) << std::endl;
    




}

