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
    int rows = A->rows();
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

void ExactConstraintSatisfaction::evaluation(Eigen::SparseMatrix<int> *A, Eigen::VectorXi *b, Eigen::VectorXi x)
{

}

void ExactConstraintSatisfaction::makeDiv(std::set<int> D, double x)
{

}

void ExactConstraintSatisfaction::safeDot(std::set<ExactConstraintSatisfaction::pair> S)
{

}







/*
// ------------------------Helpfull Methods----------------------------------------
void swapRows(Eigen::SparseMatrix<int>& A, int row1, int row2){
    int temp = 0;
    for(int i = 0; i < A.cols(); i++){
        temp = A.coeffRef(row1,i);
        A.coeffRef(row1,i) = A.coeffRef(row2, i);
        A.coeffRef(row2, i) = temp;
    }
}

int gcd(const int a, const int b){
    if(b == 0){
        return a;
    }
    return gcd(b , a % b);
}

int gcdRow(const Eigen::SparseMatrix<int>::RowXpr row, const int b){
    int gcdValue = 0;
    for(int i = 0; i < row.cols(); i ++){
        gcdValue = gcd(gcdValue, row.coeff(i));
    }
    gcdValue = gcd(gcdValue, b);
    return gcdValue;
}

void printMatrix(Eigen::SparseMatrix<int> A){
    for(int i = 0; i < A.rows(); i++){
        for(int j = 0; j < A.cols(); j++){
            std::cout << A.coeffRef(i,j) << " ";
        }
        std::cout << std::endl;
    }
}

// ----------------------Matrix transformation-----------------------------------
void IREF_Gaussian(Eigen::SparseMatrix<int>& A, Eigen::VectorXi& b){
    int rows = A.rows(); //number of rows
    int cols = A.cols(); //number of columns
    for (int k = 0;k < rows;k++) { //order Matrix after pivot
        if(k < cols){
            if(A.coeffRef(k,k) == 0){
                int l = -1;
                for (int i = k;i < rows;i++) { //find row with pivot in this column
                    if(A.coeffRef(i,k) != 0){
                        l = i;
                        return;
                    }
                }
                if(l == -1)
                    continue;
                swapRows(A,l,k); //swap rows so the pivot is in the right row
            }
        }
        for(int i = k+1; i < rows; i++){
            for(int j = k+1; j < cols; j ++){
                A.coeffRef(i,j) = A.coeffRef(k,k) * A.coeffRef(i,j) - A.coeffRef(i,k) * A.coeffRef(k,j); //eliminate the rows below the row with pivot, only one pivot per column
            }
            b.coeffRef(i) = A.coeffRef(k,k) * b.coeffRef(i) - A.coeffRef(i,k) * b.coeffRef(k);
            int gcdValue = gcdRow(A.row(i), b.coeffRef(i)); //compute the gcd to make the values as small as possible
            for(int x = 0; x < cols; x++){
                A.coeffRef(i, x) = A.coeffRef(i,x) / gcdValue;
            }
            b.coeffRef(i) = b.coeffRef(i) / gcdValue;
        }
    }
}
*/
