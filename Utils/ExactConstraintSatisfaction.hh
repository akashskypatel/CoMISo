#ifndef EXACTCONSTRAINTSATISFACTION_HH
#define EXACTCONSTRAINTSATISFACTION_HH

#include <CoMISo/Config/config.hh>
#include <CoMISo/Config/CoMISoDefines.hh>

#include <CoMISo/NSolver/NProblemInterface.hh>
#include <vector>
#include <set>

class COMISODLLEXPORT ExactConstraintSatisfaction
{
public:
    ExactConstraintSatisfaction();

    //-----------------------helpfull variables-------------------------------

    int number_pivots = 0; //number of rows with a pivot;
    int largest_exponent = 0;
    double delta = 0;

    //-----------------------helpfull methods---------------------------------

    int gcd(const int a, const int b);

    int gcdRow(const Eigen::SparseMatrix<int>::RowXpr row, const int b);

    void swapRows(Eigen::SparseMatrix<int>* A, int row1, int row2);

    void printMatrix(Eigen::SparseMatrix<int> A);
    void printVector(Eigen::VectorXi b);

    int largestExponent(Eigen::VectorXd x);
    double F_delta(double x);
    int lcm(const int a, const int b);
    int lcm_Set(const std::set<int> D);
    int indexPivot(Eigen::SparseMatrix<int>::RowXpr row);

    //--------------------matrix transformation-------------------------------

    void IREF_Gaussian(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b);
    void IRREF_Jordan(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b);

    //-------------------Evaluation--------------------------------------------

    void evaluation(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b, Eigen::VectorXd x);
    double makeDiv(const std::set<int>& D, double x);
    double safeDot(const std::set<std::pair<double, double>>& S);
};

#endif // EXACTCONSTRAINTSATISFACTION_HH
