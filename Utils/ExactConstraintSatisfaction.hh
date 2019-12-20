#ifndef EXACTCONSTRAINTSATISFACTION_HH
#define EXACTCONSTRAINTSATISFACTION_HH

#include <CoMISo/Config/config.hh>
#include <CoMISo/Config/CoMISoDefines.hh>

#include <CoMISo/NSolver/NProblemInterface.hh>
#include <vector>
#include <list>

class COMISODLLEXPORT ExactConstraintSatisfaction
{
public:
    ExactConstraintSatisfaction();

    typedef Eigen::SparseVector<int>::InnerIterator iteratorV;
    typedef Eigen::SparseVector<int> sparsVec;

    //-----------------------helpfull methods---------------------------------

    int gcd(const int a, const int b);

    int gcdRow(const Eigen::SparseMatrix<int>::RowXpr row, const int b);

    void swapRows(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b,  int row1, int row2);

    void printMatrix(Eigen::SparseMatrix<int> A);
    void printVector(Eigen::VectorXi b);

    int largestExponent(const Eigen::SparseMatrix<int>* A, const Eigen::VectorXd* x);
    double F_delta(double x);
    int lcm(const int a, const int b);
    int lcm_list(const std::list<int> D);
    int indexPivot(const Eigen::SparseMatrix<int>* A, int row);
    double get_delta();

    //--------------------matrix transformation-------------------------------

    void IREF_Gaussian(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b);
    void IRREF_Jordan(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b);

    //-------------------Evaluation--------------------------------------------

    void evaluation(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b, Eigen::VectorXd* x);
    double makeDiv(const std::list<int>& D, double x);
    double safeDot(const std::list<std::pair<int, double>>& S);

private:

    //-----------------------helpfull variables-------------------------------

    int number_pivots = 0; //number of rows with a pivot;
    int largest_exponent = 0;
    double delta = 0;
};

#endif // EXACTCONSTRAINTSATISFACTION_HH
