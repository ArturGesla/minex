#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>

void mysolve(Eigen::VectorXd &u, Eigen::SparseMatrix<double> &AA, Eigen::VectorXd &rhs)
{
    Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solve(AA);
    rhs=-1.0*rhs;
    Eigen::VectorXd p = solve.solve(rhs);
	solve.printUmfpackInfo();
    u=p;

      umfpack_di_report_info(0, 0);

}

int main()
{
    Eigen::VectorXd b(2);
    b << 1, 2;

    Eigen::SparseMatrix<double> AA(2, 2);

    AA.insert(0, 0) = 1;
    AA.insert(1, 1) = 3;

    Eigen::VectorXd x(2);
    mysolve(x, AA, b);

    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solve(AA);

    std::cout << x << std::endl;
    return 0;
}
