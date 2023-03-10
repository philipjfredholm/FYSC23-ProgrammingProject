This code is for the computer project "Time dependent quantum transport in 1D nanoarrays:
one and two electrons, and the role of electronic correlations" in the course FYSC23 Solid State Physics at Lund University.

To run the code, you must have Cern's ROOT installed. The code has been tested with ROOT 6.26. You must also have sourced 'thisroot.sh'.

To generate the plots for Task A, type "make taskA1 taskA2" while having set the pwd to the same as the Makefile. Then, type "./taskA1 && ./taskA2" which should generate the plots. 

To generate plots for Task B, the code might have to be adjusted somewhat. Most plots should however be generated by first typing "make taskB2" and then "./taskB2". 


Note also that taskB1.cpp is an old version of taskB2.cpp and will not produce the correct results.


Finally, it should also be remarked upon that it is important that SelfAdjointEigenSolver\<MatrixXd\> is used in taskB2.cpp and not EigenSolver\<MatrixXd\> or incorrect results will be produced.
