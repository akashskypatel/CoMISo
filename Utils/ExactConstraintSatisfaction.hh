#ifndef EXACTCONSTRAINTSATISFACTION_HH
#define EXACTCONSTRAINTSATISFACTION_HH

#include <CoMISo/Config/config.hh>

#include <CoMISo/NSolver/NProblemInterface.hh>
#include <vector>
#include <set>

class ExactConstraintSatisfaction
{
public:
    ExactConstraintSatisfaction();

    //-----------------------helpfull variables-------------------------------

    int number_pivots = 0; //number of rows with a pivot;
    struct pair{
        int x;
        int y;
    };

    //-----------------------helpfull methods---------------------------------

    int gcd(const int a, const int b);

    int gcdRow(const Eigen::SparseMatrix<int>::RowXpr row, const int b);

    void swapRows(Eigen::SparseMatrix<int>* A, int row1, int row2);

    void printMatrix(Eigen::SparseMatrix<int> A);
    void printVector(Eigen::VectorXi b);

    //--------------------matrix transformation-------------------------------

    void IREF_Gaussian(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b);
    void IRREF_Jordan(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b);

    //-------------------Evaluation--------------------------------------------

    void evaluation(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b, Eigen::VectorXi x);
    void makeDiv(std::set<int> D, double x);
    void safeDot(std::set<pair> S);
};

#endif // EXACTCONSTRAINTSATISFACTION_HH
