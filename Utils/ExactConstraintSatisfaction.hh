#ifndef EXACTCONSTRAINTSATISFACTION_HH
#define EXACTCONSTRAINTSATISFACTION_HH

#include <CoMISo/Config/config.hh>
#include <CoMISo/Config/CoMISoDefines.hh>

#include <CoMISo/NSolver/NProblemInterface.hh>
#include <vector>
#include <list>
#include <time.h>

class COMISODLLEXPORT ExactConstraintSatisfaction
{
public:
    ExactConstraintSatisfaction();

    typedef Eigen::SparseVector<int>::InnerIterator iteratorV;
    typedef Eigen::SparseVector<int> sparseVec;
    typedef Eigen::SparseMatrix<int, Eigen::ColMajor> SP_Matrix_C;
    typedef Eigen::SparseMatrix<int, Eigen::RowMajor> SP_Matrix_R;

    //-----------------------helpfull methods---------------------------------
    void   print_matrix(const Eigen::SparseMatrix<int, Eigen::RowMajor> A);
    void   print_vector(Eigen::VectorXi b);


    int    gcd(const int a, const int b);
    int    gcd_row(const Eigen::SparseVector<int> row, const int b);

    int    lcm(const int a, const int b);
    int    lcm(const std::vector<int>& D);

    void   swap_rows(SP_Matrix_R& mat,  int row1, int row2);
    void   eliminate_row(SP_Matrix_R& mat, int row1, int row2, int pivot_column);
    int    largest_exponent(const Eigen::VectorXd& x);
    int    index_pivot(const sparseVec& row);
    double F_delta(double x);
    double get_delta();

    //--------------------matrix transformation-------------------------------

    void   IREF_Gaussian(SP_Matrix_R& A, Eigen::VectorXi& b, const Eigen::VectorXd x);
    void   IRREF_Jordan(SP_Matrix_R& A, Eigen::VectorXi& b);

    //-------------------Evaluation--------------------------------------------

    void   evaluation(SP_Matrix_R& _A, Eigen::VectorXi& b, Eigen::VectorXd& x, const Eigen::VectorXd values);
    double makeDiv(const std::vector<int>& D, double x);
    double safeDot(const std::vector<std::pair<int, double> >& S);

private:

    // for evaluation

    int get_pivot_row_student(const SP_Matrix_C& A, int col);
    int get_pivot_row_new(const SP_Matrix_C& A, const SP_Matrix_R& _A, int col);

    std::vector<int> get_divisors_student(const SP_Matrix_C& A, int col);
    std::vector<int> get_divisors_new(const SP_Matrix_C& A, const SP_Matrix_R& _A, int col);

    std::vector<std::pair<int, double>> get_dot_product_elements_student(const SP_Matrix_C& A, const Eigen::VectorXd& x, int k, int pivot_row);
    std::vector<std::pair<int, double>> get_dot_product_elements_new(const SP_Matrix_R& A, const Eigen::VectorXd& x, int pivot_row);

    //-----------------------helpfull variables-------------------------------

    int    number_pivots_ = 0; //number of rows with a pivot;
    int    largest_exponent_ = 0;
    double delta_ = 0;
};

#endif // EXACTCONSTRAINTSATISFACTION_HH
