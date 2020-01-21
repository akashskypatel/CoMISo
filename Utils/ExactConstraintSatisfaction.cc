#include "ExactConstraintSatisfaction.hh"

#include <CoMISo/Config/config.hh>

#include <CoMISo/NSolver/NProblemInterface.hh>
#include <vector>

ExactConstraintSatisfaction::ExactConstraintSatisfaction()
{

}


// ------------------------Helpfull Methods----------------------------------------
void ExactConstraintSatisfaction::swapRows(Eigen::SparseMatrix<int>& A, Eigen::VectorXi& b, int row1, int row2){
  int temp = 0;
  for(int i = 0; i < A.cols(); i++){
    temp = A.coeffRef(row1,i);
    A.coeffRef(row1,i) = A.coeffRef(row2, i);
    A.coeffRef(row2, i) = temp;
  }

  auto b_Value = b.coeffRef(row1);
  b.coeffRef(row1) = b.coeffRef(row2);
  b.coeffRef(row2) = b_Value;
}

int ExactConstraintSatisfaction::gcd(const int a, const int b){
  if(b == 0){
    return std::abs(a);
  }
  return gcd(std::abs(b) , std::abs(a) % std::abs(b));
}

int ExactConstraintSatisfaction::gcdRow(const Eigen::SparseMatrix<int>::RowXpr row, const int b){
  int gcdValue = b;
  bool first = true;
  bool is_negativ = 0;
  for(int i = 0; i < row.cols(); i ++){
    if(row.coeff(i) != 0 && first){
      first = false;
      is_negativ = row.coeff(i) < 0;
    }
    gcdValue = gcd(gcdValue, row.coeff(i));
  }
  if(gcdValue == 0)
    return 1;
  if(is_negativ)
    gcdValue = std::abs(gcdValue) * -1;
  return gcdValue;
}

void ExactConstraintSatisfaction::printMatrix(Eigen::SparseMatrix<int> A){
  if(A.size() < 50000){
    for(int i = 0; i < A.rows(); i++){
      for(int j = 0; j < A.cols(); j++){
        std::cout << A.coeffRef(i,j) << " ";
      }
      std::cout << std::endl;
    }
  }else{
    for(int i = 0; i < A.rows(); i++){
      bool first = true;
      for(int j = 0; j < A.cols(); j++){
        auto val = A.coeffRef(i,j);
        if(val != 0 && first){
          std::cout << j << " x 0 ";
          first = false;
        }
        if(!first && val != 0)
          std::cout << val << " ";
      }
      if(!first){
        std::cout << std::endl;
        std::cout << std::endl;
      }

    }
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

int ExactConstraintSatisfaction::largestExponent(const Eigen::SparseMatrix<int>& A, const Eigen::VectorXd& x)
{
  int expo = -65;
  bool empty_col= true;
  for(int i = 0; i < x.size(); i++){
    for(int j = 0; j < A.rows(); j++){
      if(A.coeff(j,i) != 0){
        empty_col = false;
        break;
      }
    }
    if(!empty_col)
      expo = std::max(expo, static_cast<int>(std::ceil(std::log2(std::abs(x.coeffRef(i)))) + 2));
  }
  largest_exponent_ = expo;
  delta_ = std::pow(2, expo);
  return expo;
}

double ExactConstraintSatisfaction::F_delta(double x)
{
  int sign = -1;
  if(x >= 0)
    sign = 1;
  double x_of_F = (x + sign * delta_) - sign * delta_;
  return x_of_F;
}

int ExactConstraintSatisfaction::lcm(const int a, const int b)
{
  if(gcd(a,b) == 0)
    return 0;
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

int ExactConstraintSatisfaction::indexPivot(const Eigen::SparseMatrix<int>& A, int row)
{
  auto row1 = A.row(row);
  for(int i = 0; i < row1.size(); i++){
    if(row1.coeff(i) != 0)
      return i;
  }
  return -1;
}

double ExactConstraintSatisfaction::get_delta()
{
  return delta_;
}

// ----------------------Matrix transformation-----------------------------------
void ExactConstraintSatisfaction::IREF_Gaussian(Eigen::SparseMatrix<int>& A, Eigen::VectorXi& b){

  number_pivots_ = 0;
  int rows = A.rows();        //number of rows
  int cols = A.cols();        //number of columns
  int col_index = -1;         //save the last column where we found a pivot

  for (int k = 0; k < rows; k++) {                                        //order Matrix after pivot

    //needed for the first line
    if(k == 0){
      int gcdValue = gcdRow(A.row(k), b.coeffRef(k));                     //compute the gcd to make the values as small as possible
      for(int x = 0; x < cols; x++){
        A.coeffRef(k, x) = A.coeffRef(k,x) / gcdValue;
      }
      b.coeffRef(k) = b.coeffRef(k) / gcdValue;
    }

    if(k < cols){
      if(A.coeffRef(k,k) == 0){
        int pivot_row = -1;
        for(col_index += 1; col_index < cols; col_index++){              //find the smallest column with a pivot
          for (int i = k; i < rows; i++) {                               //find row with pivot in this column
            if(A.coeffRef(i,col_index) != 0){
              pivot_row = i;
              number_pivots_++;
              break;
            }
          }
          if(pivot_row != -1)
            break;
        }
        if(pivot_row == -1)
          continue;
        if(pivot_row != k)
          swapRows(A, b, pivot_row, k);                                  //swap rows so the pivot is in the right row
      }else{
        col_index ++;
        number_pivots_++;
      }
    }
    for(int i = k+1; i < rows; i++){
      int under_pivot = A.coeffRef(i,col_index);
      for(int j = col_index; j < cols; j++){                                                          //change k+1 to k to eliminate the elements below the pivot
        A.coeffRef(i,j) = A.coeffRef(k,col_index) * A.coeffRef(i,j) - under_pivot * A.coeffRef(k,j);  //eliminate the rows below the row with pivot, only one pivot per column
      }
      b.coeffRef(i) = A.coeffRef(k,col_index) * b.coeffRef(i) - under_pivot * b.coeffRef(k);
    }

    for(int i = k; i < rows; i++){
      int gcdValue = gcdRow(A.row(i), b.coeffRef(i));                                                 //compute the gcd to make the values as small as possible
      for(int j = 0; j < cols; j++){
        A.coeffRef(i, j) = A.coeffRef(i,j) / gcdValue;
      }
      b.coeffRef(i) = b.coeffRef(i) / gcdValue;
    }

  }
}

void ExactConstraintSatisfaction::IRREF_Jordan(Eigen::SparseMatrix<int> &A, Eigen::VectorXi &b)
{
  int cols = A.cols();
  for(int k = number_pivots_ - 1; k > 0; k--){
    int p_col = k;
    for(p_col = k; p_col < cols; p_col++){              //find pivot
      if(A.coeffRef(k,p_col) != 0)
        break;
    }
    for(int i = k-1; i >= 0; i--){                      //eliminate row i with row k
      int above_p = A.coeffRef(i,p_col);                //element in rows above the pivot
      for(int j = cols-1; j >= i; j--){
        if(j >= p_col){
          A.coeffRef(i,j) = A.coeffRef(k,p_col) * A.coeffRef(i,j) - above_p * A.coeffRef(k,j);
        }else if(j >= i){
          A.coeffRef(i,j) = A.coeffRef(k,p_col) * A.coeffRef(i,j);
        }
      }
      b.coeffRef(i) = A.coeffRef(k,p_col) * b.coeffRef(i) - above_p * b.coeffRef(k);

      int gcdValue = gcdRow(A.row(i), b.coeffRef(i));   //compute the gcd to make the values as small as possible
      for(int j = 0; j < cols; j++){
        A.coeffRef(i, j) = A.coeffRef(i,j) / gcdValue;
      }
      b.coeffRef(i) = b.coeffRef(i) / gcdValue;
    }

  }
}

//-------------------------------------Evaluation--------------------------------------------------------------------

void ExactConstraintSatisfaction::evaluation(Eigen::SparseMatrix<int>& A, Eigen::VectorXi& b, Eigen::VectorXd& x)
{
  IREF_Gaussian(A, b);
  IRREF_Jordan(A, b);

  int cols = A.cols();
  std::cout << "largest Expo" << std::endl;
  largestExponent(A, x);
  std::cout << "findPivo" << std::endl;
  for(int k = cols -1; k >= 0; k--){
    int pivot = -1;

//new faster Version test, iterates over the non zero (non empty) elements of one col
//    for(iteratorV it(static_cast<sparsVec>(A.col(k))); it; ++it){
//      int i = it.index();


      for(int i = 0; i <= k; i++){                      //find row with pivot in this column
      if(i < A.rows()){
        int pivotIndex = indexPivot(A,i);
        if(pivotIndex == k){
          pivot = i;                                    //the row with the pivot
          break;
        }
      }

    }

    if(pivot == -1){                                    //there is no pivot in this column
      std::list<int> D;
      D.clear();
      for(int i = 1; i <= std::min(k, number_pivots_ - 1); i ++){
        int pivot_col = indexPivot(A, i);
        if(A.coeffRef(i,k) != 0 && pivot_col < k){
          D.push_front(A.coeffRef(i, pivot_col));
        }
      }
      x.coeffRef(k) = makeDiv(D, x.coeffRef(k));        //fix free variables so they are in F_delta


    }else{                                              //compute now the implied variables

      std::list<std::pair<int, double>> S;
      S.clear();
      for(int i = k+1; i < cols; i++){                  //construct the list S to do the dot Product
        std::pair<int, double> tuple;
        tuple.first = A.coeffRef(pivot,i);
        double test =x.coeffRef(i) / A.coeffRef(pivot,k);
        if(x.coeffRef(i) != ( test * A.coeffRef(pivot,k)))
           std::cout << "WARNING: can't devide" << " in row : " << i << std::endl;
        tuple.second = x.coeffRef(i) / A.coeffRef(pivot,k);
        if(tuple.first != 0)
          S.push_front(tuple);

      }
      double divided_B = b.coeffRef(pivot);
      divided_B = F_delta(divided_B / A.coeffRef(pivot, k));
      if( divided_B * A.coeffRef(pivot, k) !=  static_cast<double>(b.coeffRef(pivot)))
        std::cout << "WARNING: Can't handle the right hand side perfectly" << std::endl;
      x.coeffRef(k) = divided_B - safeDot(S);
    }
  }
}

double ExactConstraintSatisfaction::makeDiv(const std::list<int>& D, double x)
{
  if(D.empty()){
    return F_delta(x);
  }
  int d = lcm_list(D);
  double result = F_delta(x/d) * d;
  return result;

}

double ExactConstraintSatisfaction::safeDot(const std::list<std::pair<int, double> >& S)
{
  if(S.empty()) //case all Cij are zero after the pivot
    return 0;
  int safebreak = 9999999;
  std::list<std::pair<int, double>> P, N;
  P.clear();
  N.clear();
  int k = 0;
  double r = 0;                                         //return value of the dot

  for(std::pair<int, double> element : S){
    if(element.first * element.second > 0)
    {
      element.first = std::abs(element.first);
      element.second = std::abs(element.second);
      P.push_front(element);
    }
    else if(element.first * element.second < 0)
    {
      element.first = std::abs(element.first);
      element.second = - std::abs(element.second);
      N.push_front(element);
    }
  }

  while((!P.empty() || !N.empty()) && safebreak > 0){

    safebreak--;                                        //break out of the while after an amount of time
    if(!P.empty() && (r < 0 || N.empty())){
      const std::pair<int, double> element = P.front();
      P.pop_front();
      k = std::min(element.first, static_cast<int>(std::floor((delta_ - r)/element.second)));
      if(k <=0)
        std::cout << "k == 0" << std::endl;
      r = r + k * element.second;
      if(k < element.first)
        P.push_front({element.first - k, element.second});

    }else{

      const std::pair<int, double> element = N.front();
      N.pop_front();
      k = std::min(element.first, static_cast<int>(std::floor((-delta_ - r)/element.second)));
      if(k <= 0)
        std::cout << "k == 0" << std::endl;
      r = r + k * element.second;
      if(k < element.first)
        N.push_front({element.first - k, element.second});
    }
  }
  if(safebreak == 0)
    std::cout << "WARNING:: the set number of Iterations in Method : safedot is exceeded" << std::endl;
  return r;
}

