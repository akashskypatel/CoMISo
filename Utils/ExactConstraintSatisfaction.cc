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
  //std::cout << "method largest expo, size of vector :" << x.size() << std::endl;
    int expo = 0;
    for(int i = 0; i < x.size(); i++){
        expo = std::max(expo, (int)std::ceil(std::abs(std::log2(std::abs(x.coeffRef(i))))) + 1); //added a norm around the log2 because it can get negative and ruin the max
        //std::cout << "vecor x at posion i: " << std::log2(x.coeffRef(i)) << std::endl;
        //std::cout << "expo = " << expo << std::endl;
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

int ExactConstraintSatisfaction::lcm_list(const std::list<int> D)
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
          int under_pivot = A->coeffRef(i,k);
            for(int j = k; j < cols; j ++){                                                                //change k+1 to k to eliminate the elements below the pivot
              std::cout << " outer for i = " << i << "inner for j = " << j << "set A(i, j) = " << A->coeffRef(k,k) << " * " << A->coeffRef(i,j) << " - " << A->coeffRef(i,k) << " * " << A->coeffRef(k,j) << std::endl;
                A->coeffRef(i,j) = A->coeffRef(k,k) * A->coeffRef(i,j) - under_pivot * A->coeffRef(k,j); //eliminate the rows below the row with pivot, only one pivot per column

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

void ExactConstraintSatisfaction::evaluation(Eigen::SparseMatrix<int> *A, Eigen::VectorXi *b, Eigen::VectorXd *x)
{
    int cols = A->cols();
    largestExponent(*x);
    for(int k = cols -1; k >= 0; k--){
      std::cout << "evaluation for column : " << k << std::endl;
        int pivot = -1;

        for(int i = 0; i <= k; i++){     //find row with pivot in this column
          if(i < A->rows()){
            std::cout << "find row in collumn, row : " << i << std::endl;
            int pivotIndex = indexPivot(A->row(i));
            if(pivotIndex == k){
            std::cout << "pivot found in : " << i << std::endl;
               pivot = i;           //the row with the pivot
               break;
           }
          }

        }

        if(pivot == -1){      //there is no pivot in this column
          std::cout << "no pivot" << std::endl;
            std::list<int> D;
            D.clear();
            for(int i = 0; i <= std::min(k, number_pivots - 1); i ++){
              std::cout << "collect divisors for : " << i << std::endl;
                if(A->coeffRef(i,k) != 0){
                    D.push_front(A->coeffRef(i, indexPivot(A->row(i))));
                }
            }
            std::cout << "fix free variables : "<< std::endl;
            x->coeffRef(k) = makeDiv(D, x->coeffRef(k)); //fix free variables so they are in F_delta
        }else{                                      //compute now the implied variables
          std::cout << "pivot found"<< std::endl;
            std::list<std::pair<double, double>> S;
            S.clear();
            for(int i = k+1; i < cols; i++){        //construct the list S to do the dot Product
              std::cout << "construct S"<< std::endl;
                std::pair<double, double> tuple;
                tuple.first = A->coeffRef(pivot,i);
                tuple.second = x->coeffRef(i) / A->coeffRef(pivot,k);
                std::cout << "add to S (c, x) : c = " << tuple.first << " x = " << tuple.second << std::endl;
                if(tuple.first != 0)
                  S.push_front(tuple);
            }
            std::cout << "compute implied vaiables"<< std::endl;
            double divided_B = b->coeffRef(k);
            divided_B = F_delta(divided_B/A->coeffRef(pivot, k));
            if( divided_B != ((double)b->coeffRef(k) / A->coeffRef(pivot, k)))
              std::cout << "WARNING: Can't handle the right hand side perfectly" << std::endl;
            std::cout << "divided value of b = " << divided_B << std::endl;
            x->coeffRef(k) = divided_B - safeDot(S);
        }
    }
}

double ExactConstraintSatisfaction::makeDiv(const std::list<int>& D, double x)
{
  std::cout << "make div "<< x << std::endl;
  if(D.empty()){
    return F_delta(x);
  }
  int d = lcm_list(D);
  double result = F_delta(x/d) * d;
  std::cout << "the result of make div is : " << result << std::endl;
  return result;

}

double ExactConstraintSatisfaction::safeDot(const std::list<std::pair<double, double>>& S)
{
  if(S.empty()) //case all Cij are zero after the pivot
    return 0;
  int safebreak = 5;
  std::cout << "safe dot"<< std::endl;
  std::cout << "delta is " << delta << " and max expo is " << largest_exponent << std::endl;
  std::list<std::pair<double, double>> P, N;
  P.clear();
  N.clear();
  int k = 0; //eventuelly a double not sure
  double r = 0;          //return value of the dot
  for(std::pair<double, double> element : S){
    if(element.first * element.second > 0)
    {
      element.first = std::abs(element.first);
      element.second = std::abs(element.second);
      P.push_front(element);
    }
    else
    {
      element.first = std::abs(element.first);
      element.second = - std::abs(element.second);
      N.push_front(element);
    }
  }
  while((!P.empty() || !N.empty()) && safebreak > 0){
    safebreak--; //break out of the while after an amount of time
    std::cout << "size of P, N : " << P.size() << N.size() << std::endl;
    if(!P.empty() && (r < 0 || N.empty())){
      const std::pair<double, double> element = P.front();
      P.pop_front();
      k = std::min(element.first, std::floor((delta - r)/element.second));
      r = r + k * element.second;
      if(k < element.first)
        P.push_front({element.first - k, element.second});
    }else{
      const std::pair<double, double> element = N.front();
      std::cout << "element x, y :" << element.first << ", " << element.second << std::endl;
      std::cout << "r = " << r << " k = " << k << std::endl;
      N.pop_front();
      std::cout << "size of P, N : " << P.size() << N.size() << std::endl;
      k = std::min(element.first, std::floor((-delta - r)/element.second));
      r = r + k * element.second;
      std::cout << "r = " << r << " k = "<< k << std::endl;
      if(k < element.first)
        N.push_front({element.first - k, element.second});
    }
  }
  std::cout << "value r = " << r << std::endl;
  return r;
}

