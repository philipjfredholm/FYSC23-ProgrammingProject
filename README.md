# FYSC23-ProgrammingProject






    /*
    Now we want to find the initial values c_m(0) by minimising the energy.
    In general, this seems like a hard problem as we need to both minimise
    \sum E_n |c_n|^2 and keep the normalisation constraint on \sum |c_n|^2.
    However, that we have negative eigenvalues simplifies this greatly as the 
    minimum is then obviously when we set all other c_m to 0 and the one with the 
    most negative eigenvalue to magnitude 1 (or just 1 since we have global gauge invariance).
    */

    /*
    Since the energy-eigenstates basis \varphi is already given in the basis \phi,
    all we need to do is normalise the eigenvector correspond to the minimum (negative)
    eigenenergy to find c_m(0)
    */

    /*
    We now want c_m^\lambda. Since the energy-eigenbasis is already given in the desired basis,
    <\phi_m|\varphi_\lambda> is just the n:th row of the identity matrix mutliplied with 
    the \lambda:th eigen(column)vector. (In other words we just set all elements except m to 0).

    We need to keep in mind that since epsilon is not small (\pm 2) compared to V which is -1,
    this does not just peturb the energy-eigenstate. We get the initial c_m() from the
    matrix in equation (28) in the manual, but we need to project it onto the new energy-eigen basis once
    the heaviside theta-function kicks in. This is what equation (11) in the manual does. To do this, we
    use that both of these are given in the \phi basis, which is orthogonal. Hence,
    all we need to do is take the dot product.
    */
VectorXd makeReal(VectorXcd complexVector) {
    /* 
    The Eigen software gives back a vector
    of possibly complex eigenvalues. This just 
    makes that vector real and double checks so that
    they are not complex 
    */

    int length = complexVector.size();
    VectorXd realVector(length);

    for (int n = 0; n < length; n++) {
        if (complexVector(n).imag() != 0) {
            std::cout << "Warning: non-real eigenvalues detected!" << std::endl;
        }

        realVector(n) =  complexVector(n).real();
    }

    return realVector;

}



void makePlot(EigenSolver<MatrixXd> hamiltonian, VectorXd initialValues,
            double timeInterval, double timeStepLength, int argc, char** argv) {

    int length = initialValues.size();
    const char* title = "Values of the coefficients at #t = ";
    std::string stringTitle = "Values of the coefficients at #t = ";

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
    int million = 1000000;
    std::complex<double> coefficient;
    while (time <= timeInterval) {
        hist->Reset();
        for (int k = 1; k <= length; k++) { //Start at 1 to not go into underflow bin
            coefficient = coefficients(k-1, time, hamiltonian, initialValues);
            double magnitude = pow(std::abs(coefficient),2);
            //magnitude *= 1000; //Histogram only acccepts integer entries.
            //magnitude = std::ceil(magnitude);
            hist->SetBinContent(k, magnitude);

        }


        //std::string timeString = stringTitle + std::to_string(time);
        //const char* newTitle = timeString.c_str();
        //hist->SetTitle(newTitle);

        hist->Draw();
        canvas->Update();

        usleep(timeStepLength*million); //usleep() works in microseconds
        time += timeStepLength;
    }


    
    
    application->Run();


}