#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <Eigen/UmfPackSupport>

#include "lnmatrxe.h"
#include "lnmatrxf.h"
#include "arlnsmat.h"
#include "arlgnsym.h"
#include "lnsymsol.h"

void mysolve(Eigen::VectorXd &u, Eigen::SparseMatrix<double> &AA, Eigen::VectorXd &rhs)
{
    // Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solve(AA);
    // rhs=-1.0*rhs;
    // Eigen::VectorXd p = solve.solve(rhs);
    // solve.printUmfpackInfo();
    // u=p;

    //   umfpack_di_report_info(0, 0);
}
void findev(Eigen::SparseMatrix<double> &AA, Eigen::SparseMatrix<double> &BB)
{
}

int main()
{
    Eigen::VectorXd b(2);
    b << 1, 2;

    Eigen::SparseMatrix<double> AA(2, 2);

    AA.insert(0, 0) = 1;
    AA.insert(1, 1) = 3;

    Eigen::SparseMatrix<double> BB(2, 2);

    BB.insert(0, 0) = 1;
    BB.insert(1, 1) = 1;

    Eigen::VectorXd x(2);
    mysolve(x, AA, b);

    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solve(AA);

    std::cout << x << std::endl;

    //Defining variables;

    int n;               // Dimension of the problem.
    int nnza, nnzb;      // Number of nonzero elements in A and B.
    int *irowa, *irowb;  // pointers to arrays that store the row
                         // indices of the nonzeros in A and B.
    int *pcola, *pcolb;  // pointers to arrays of pointers to the
                         // beginning of each column of A and B in
                         // valA and valB.
    double *valA, *valB; // pointers to arrays that store the
                         // nonzero elements of A and B.

    // Creating matrices A and B.

    n = 100;
    //n = AA.rows();
    NonSymMatrixE(n, nnza, valA, irowa, pcola, AA);
    ARluNonSymMatrix<double, double> A(n, nnza, valA, irowa, pcola);

    NonSymMatrixF(n, nnzb, valB, irowb, pcolb, BB);
    ARluNonSymMatrix<double, double> B(n, nnzb, valB, irowb, pcolb);

    // Defining what we need: the four eigenvectors nearest to 0.4 + 0.6i.

    ARluNonSymGenEig<double> dprob(4L, A, B, 'R', 0.4, 0.0);

    // Finding eigenvalues and eigenvectors.

    dprob.FindEigenvectors();

    // Printing solution.

    Solution(A, B, dprob);

    return 0;
}
