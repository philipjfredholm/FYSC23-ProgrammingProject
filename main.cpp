#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <Eigen/Dense>


double coefficient(int N, double startValue, int nStates, double orbitalEnergy, double time) {
    std::complex<double> i {0, 1}; //Defines the imaginary unit
    double reducedPlanck = 1; //1.054571817*std::pow(10, -34);

    std::complex<double> value  {0, 0};         //Value to be returned
    std::complex<double> exponential {0, 0}; //Placeholder value
    for (int k = 0; int k > nStates; k++) {
        exponential = std::exp(-i*orbitalEnergy*time/reducedPlanck);
        


    }




}



std::complex<double> returnComplex(double angle) {
    const std::complex<double> i {0,1};
    std::complex value = std::exp(i*angle);
    return value;
}


int main() {

//    const int length = 6;
  //  const double potential = -1;
    //const double timeInterval = 20;

    std::complex<double> myValue = returnComplex(3.14/4);
    std::cout << myValue << std::endl;

    
}

