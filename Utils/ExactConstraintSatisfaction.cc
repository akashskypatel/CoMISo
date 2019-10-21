#include "ExactConstraintSatisfaction.hh"

#include <CoMISo/Config/config.hh>

#include <CoMISo/NSolver/NProblemInterface.hh>
#include <vector>

ExactConstraintSatisfaction::ExactConstraintSatisfaction()
{

}


// ------------------------Helpfull Methods----------------------------------------
void ExactConstraintSatisfaction::swapRows(Eigen::SparseMatrix<int>* A, int row1, int row2){
    int temp = 0;
    for(int i = 0; i < A->cols(); i++){
        temp = A->coeffRef(row1,i);
        A->coeffRef(row1,i) = A->coeffRef(row2, i);
        A->coeffRef(row2, i) = temp;
    }
}

int ExactConstraintSatisfaction::gcd(const int a, const int b){
    if(b == 0 && a == 0) //to prevent dividing by 0 later
        return 1;
    if(b == 0){
        return a;
    }
    return gcd(b , a % b);
}

int ExactConstraintSatisfaction::gcdRow(const Eigen::SparseMatrix<int>::RowXpr row, const int b){
    int gcdValue = 0;
    for(int i = 0; i < row.cols(); i ++){
        gcdValue = gcd(gcdValue, row.coeff(i));
    }
    gcdValue = gcd(gcdValue, b);
    return gcdValue;
}

void ExactConstraintSatisfaction::printMatrix(Eigen::SparseMatrix<int> A){
    for(int i = 0; i < A.rows(); i++){
        for(int j = 0; j < A.cols(); j++){
            std::cout << A.coeffRef(i,j) << " ";
        }
        std::cout << std::endl;
    }
}

void ExactConstraintSatisfaction::printVector(Eigen::VectorXi b)
{
    std::cout << "the vector contains the elements: ";
    for(int i = 0; i < b.size(); i++){
        std::cout << b.coeffRef(i) << " ";
    }
    std::cout << std::endl;
}

int ExactConstraintSatisfaction::largestExponent(Eigen::VectorXd x)
{
    int expo = 0;
    for(int i = 0; i < x.size(); i++){
        expo = std::max(expo, (int)std::ceil(std::log2(std::abs(x.coeffRef(i)))) + 1);
    }
    largest_exponent = expo;
    delta = std::pow(2, expo);
    return expo;
}

double ExactConstraintSatisfaction::F_delta(double x)
{
    int sign = -1;
    if(x >= 0)
        sign = 1;
    double x_of_F = (x + sign * delta) - sign * delta;
    return x_of_F;
}

int ExactConstraintSatisfaction::lcm(const int a, const int b)
{
    return (a/gcd(a,b)) * b;
}

int ExactConstraintSatisfaction::lcm_Set(const std::set<int> D)
{
    int lcm_D = 1;
    for(int d : D){
        lcm_D = lcm(lcm_D, d);
    }
    return lcm_D;
}

int ExactConstraintSatisfaction::indexPivot(Eigen::SparseMatrix<int>::RowXpr row)
{
    for(int i = 0; i < row.size(); i++){
        if(row.coeffRef(i) != 0)
            return i;
    }
    return -1;
}

// ----------------------Matrix transformation-----------------------------------
void ExactConstraintSatisfaction::IREF_Gaussian(Eigen::SparseMatrix<int>* A, Eigen::VectorXi* b){
    number_pivots = 0;
    int rows = A->rows(); //number of rows
    int cols = A->cols(); //number of columns
    for (int k = 0; k < rows;k++) { //order Matrix after pivot

        //needed for the first line
        int gcdValue = gcdRow(A->row(k), b->coeffRef(k)); //compute the gcd to make the values as small as possible
        for(int x = 0; x < cols; x++){
            A->coeffRef(k, x) = A->coeffRef(k,x) / gcdValue;
        }
        b->coeffRef(k) = b->coeffRef(k) / gcdValue;


        if(k < cols){
            if(A->coeffRef(k,k) == 0){
                int l = -1;
                for (int i = k +1;i < rows;i++) { //find row with pivot in this column
                    if(A->coeffRef(i,k) != 0){
                        l = i;
                        number_pivots++;
                        break;
                    }
                }
                if(l == -1)
                    continue;
                swapRows(A,l,k); //swap rows so the pivot is in the right row
            }else{
                number_pivots++;
            }
        }
        for(int i = k+1; i < rows; i++){
            for(int j = k; j < cols; j ++){                                                                 //change k+1 to k to eliminate the elements below the pivot
                A->coeffRef(i,j) = A->coeffRef(k,k) * A->coeffRef(i,j) - A->coeffRef(i,k) * A->coeffRef(k,j); //eliminate the rows below the row with pivot, only one pivot per column
            }
            b->coeffRef(i) = A->coeffRef(k,k) * b->coeffRef(i) - A->coeffRef(i,k) * b->coeffRef(k);
            int gcdValue = gcdRow(A->row(i), b->coeffRef(i)); //compute the gcd to make the values as small as possible
            for(int x = 0; x < cols; x++){
                A->coeffRef(i, x) = A->coeffRef(i,x) / gcdValue;
            }
            b->coeffRef(i) = b->coeffRef(i) / gcdValue;
        }
    }
}

void ExactConstraintSatisfaction::IRREF_Jordan(Eigen::SparseMatrix<int> *A, Eigen::VectorXi *b)
{
    int cols = A->cols();
    for(int k = number_pivots - 1; k > 0; k--){
        int l = k;
        for(l = k; l < cols; l++){                //find pivot
            if(A->coeffRef(k,l) != 0)
                break;
        }
        for(int i = k-1; i >= 0; i--){      //eliminate row i with row k
            int c = A->coeffRef(i,l);
            for(int j = cols-1; j >= i; j--){
                if(j >= k){
                    A->coeffRef(i,j) = A->coeffRef(k,l) * A->coeffRef(i,j) - c * A->coeffRef(k,j);
                }else{
                    A->coeffRef(i,j) = A->coeffRef(k,l) * A->coeffRef(i,j);
                }
            }
            b->coeffRef(i) = A->coeffRef(k,l) * b->coeffRef(i) - c * b->coeffRef(k);
            int gcdValue = gcdRow(A->row(i), b->coeffRef(i)); //compute the gcd to make the values as small as possible
            for(int x = 0; x < cols; x++){
                A->coeffRef(i, x) = A->coeffRef(i,x) / gcdValue;
            }
            b->coeffRef(i) = b->coeffRef(i) / gcdValue;
        }
    }
}

//-------------------------------------Evaluation--------------------------------------------------------------------

void ExactConstraintSatisfaction::evaluation(Eigen::SparseMatrix<int> *A, Eigen::VectorXi *b, Eigen::VectorXd x)
{
    int cols = A->cols();
    largestExponent(x);
    for(int k = cols -1; k >= 0; k--){
        int pivot = -1;

        for(int i = 0; i <= k; i++){     //find row with pivot in this column
           int pivotIndex = indexPivot(A->row(i));
           if(pivotIndex == k){
               pivot = i;           //the row with the pivot
               break;
           }

        }

        if(pivot == -1){      //there is no pivot in this column
            std::set<int> D;
            D.clear();
            for(int i = 0; i <= std::min(k, number_pivots - 1); i ++){
                if(A->coeffRef(i,k) != 0){
                    D.insert(A->coeffRef(i, indexPivot(A->row(i))));
                }
            }
            x.coeffRef(k) = makeDiv(D, x.coeffRef(k)); //fix free variables so they are in F_delta
        }else{                                      //compute now the implied variables
            std::set<std::pair<double, double>> S;
            S.clear();
            for(int i = k+1; i < cols; i++){        //construct the Set S to do the dot Product
                std::pair<double, double> tuple;
                tuple.first = A->coeffRef(pivot,i);
                tuple.second = x.coeffRef(i) / A->coeffRef(pivot,k);
                S.insert(tuple);
            }
            x.coeffRef(k) = b->coeffRef(k) - safeDot(S);
        }
    }
}

double ExactConstraintSatisfaction::makeDiv(const std::set<int>& D, double x)
{
  if(D.empty()){
    return F_delta(x);
  }
  int d = lcm_Set(D);
  return F_delta(x/d) * d;

}

double ExactConstraintSatisfaction::safeDot(const std::set<std::pair<double, double>> &S)
{
  std::set<std::pair<double, double>> P, N;
  P.clear();
  N.clear();
  for(std::pair<double, double> element : S){
    if(element.first * element.second > 0)
    {
      P.insert(element);
    }
    else
    {
      N.insert(element);
    }
  }
  double r = 0;
  while(!P.empty() || !N.empty()){
    if(!P.empty() && (r < 0 || N.empty())){
      //std::pair<double, double> element = P.begin();

    }
  }
}

