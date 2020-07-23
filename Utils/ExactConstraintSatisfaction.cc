#include "ExactConstraintSatisfaction.hh"

#include <CoMISo/Config/config.hh>

#include <CoMISo/NSolver/NProblemInterface.hh>
#include <CoMISo/Utils/CoMISoError.hh>
#include <vector>

ExactConstraintSatisfaction::ExactConstraintSatisfaction()
{

}


// ------------------------Helpfull Methods----------------------------------------

//row1 belongs to vector b, row2 to a row in the matrix
void ExactConstraintSatisfaction::swap_rows(SP_Matrix_R& mat, int row1, int row2){

  Eigen::SparseVector<int> row_2 = mat.row(row2);
  Eigen::SparseVector<int> row_1 = mat.row(row1);
  mat.row(row2) = row_1;
  mat.row(row1) = row_2;

//  mat.prune(0.0, 0);
  mat.makeCompressed();
  mat.finalize();

}


//We want to eliminate row1 in mat with the row corresponding to row2 in mat
//the row_2 has a pivot in (row2, col_p)
void ExactConstraintSatisfaction::eliminate_row(SP_Matrix_R& mat, int row1, int row2, int pivot_column)
{

//  int pivot_column = -1;
  const sparseVec row_2 = mat.row(row2);

//  for(sparseVec::InnerIterator it(row_2); it; ++it){
//    if(it.value() != 0){
//      pivot_column = it.index();
//      break;
//    }
//  }

//  if(counter != pivot_column)
//    std::cout << "UNGLEICH!!!!!!!!!!!!!!!!!!!!! " << counter << ", " << pivot_column << std::endl;

  if(pivot_column == -1)
    std::cout << "Error in eliminate_row (ExactConstraintSatsfaction.cc) : expected a pivot but didn't find any." << std::endl;

  int pivot_row1 = mat.coeff(row1, pivot_column);     //the element under the pivot

  if(pivot_row1 == 0)
    return;

  //declination here, to reduce runtime
  const sparseVec row_1 = mat.row(row1);
  int pivot_row2 = mat.coeff(row2, pivot_column);     //the pivot


//  if(counter == 109)
//  std::cout << mat << std::endl << std::endl << row_1 << std::endl << row_2 << std::endl;


  for(sparseVec::InnerIterator it(row_1); it; ++it){
    int col = it.col();
    int index = it.index();
    mat.coeffRef(row1, it.index()) = (pivot_row2 * it.value());
  }

  mat.makeCompressed();
  mat.finalize();

  for(sparseVec::InnerIterator it(row_2); it; ++it){
    mat.coeffRef(row1, it.index()) = mat.coeff(row1, it.index()) - (pivot_row1 * it.value());
  }

  mat.prune(0.0, 0);
  mat.makeCompressed();
  mat.finalize();
}

int ExactConstraintSatisfaction::gcd(const int a, const int b){
  if(b == 0){
    return std::abs(a);
  }
  return gcd(std::abs(b) , std::abs(a) % std::abs(b));
}

int ExactConstraintSatisfaction::gcd_row(const Eigen::SparseVector<int> row, const int b){

  int gcdValue = b;
  bool first = true;
  bool is_negativ = 0;

  for(Eigen::SparseVector<int>::InnerIterator it(row); it; ++it){
    if(it.value() != 0 && first){
      first = false;
      is_negativ = it.value() < 0;
    }
    gcdValue = gcd(gcdValue, it.value());
  }
  if(gcdValue == 0)
    return 1;
  if(is_negativ)
    gcdValue = std::abs(gcdValue) * -1;
  return gcdValue;
}

void ExactConstraintSatisfaction::print_matrix(const Eigen::SparseMatrix<int, Eigen::RowMajor> A){

  for(int i = 0; i < A.rows(); i++){
    std::cout << "row: " << i << " ";
    for(SP_Matrix_R::InnerIterator it(A, i); it; ++it){
      std::cout << "(" << it.index() << ", " << it.value() << "), ";
    }

    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void ExactConstraintSatisfaction::print_vector(Eigen::VectorXi b)
{
  std::cout << "the vector contains the elements: ";
  for(int i = 0; i < b.size(); i++){
    std::cout << b.coeffRef(i) << " ";
  }
  std::cout << std::endl;
}

int ExactConstraintSatisfaction::largest_exponent(const Eigen::VectorXd& x)
{

  int expo = -65;
  for(int i = 0; i < x.size(); i++){
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

int ExactConstraintSatisfaction::lcm(const std::vector<int>& D)
{
  int lcm_D = 1;
  for(int d : D)
  {
    lcm_D = lcm(lcm_D, d);
  }
  return lcm_D;
}

int ExactConstraintSatisfaction::index_pivot(const sparseVec& row)
{

  for(sparseVec::InnerIterator it(row); it; ++it)
  {
    if(it.value() != 0)
      return it.index();
  }
  return -1;
}

double ExactConstraintSatisfaction::get_delta()
{
  return delta_;
}

// ----------------------Matrix transformation-----------------------------------
void ExactConstraintSatisfaction::IREF_Gaussian(SP_Matrix_R& A, Eigen::VectorXi& b, const Eigen::VectorXd x){

  number_pivots_ = 0;
  int rows = A.rows();        //number of rows
  int cols = A.cols();        //number of columns
  int col_index = 0;         //save the next column where we search for a pivot

  for (int k = 0; k < rows; k++)
  {                                        //order Matrix after pivot

    //needed for the first line
    if(k == 0){
      int gcdValue = gcd_row(A.row(k), b.coeffRef(k));                     //compute the gcd to make the values as small as possible
      if(gcdValue == 1)
        continue;

      for(Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(A, k); it; ++it)
      {
        it.valueRef() = (it.value() / gcdValue);
      }
      b.coeffRef(k) = b.coeffRef(k) / gcdValue;
    }

    if(A.coeff(k, col_index) == 0)
    {

      int pivot_row = -1;

      for(; col_index < cols; col_index++){

        Eigen::SparseVector<int> col = A.col(col_index);
        for(Eigen::SparseVector<int>::InnerIterator it(col); it; ++it)
        {
          if(it.value() != 0 && it.index() >= k)
          {
            pivot_row = it.index();
            break;
          }
        }
        if(pivot_row != -1)
        {
          if(k != pivot_row)
          {
            swap_rows(A, k, pivot_row);
            //std::cout << "error is in Gaus at " << counter << ": " << (A.cast<double>() * x).squaredNorm() << std::endl;
          }
          col_index++;
          break;
        }
      }

      if(col_index == cols)
        break;

      if(pivot_row == -1)
        continue;
    }
    else
    {
      col_index++;//col_index = k +1;
    }
    number_pivots_++;                                                                       //we have found a pivot in column: col_index - 1, and row: pivot_row

    int col_p = col_index -1;

    for(int i = k+1; i < rows; ++i)
    {
      if(A.coeff(i, col_p) == 0)
        continue;

      b.coeffRef(i) = (A.coeff(k, col_p) * b.coeff(i) - A.coeff(i, col_p) * b.coeff(k));    //so we don't delete entries in A for the computation of b
      eliminate_row(A, i, k, col_p);     //Robin  : maybe write the method directly here
    }

    for(int i = k; i < rows; i++)
    {

      int gcdValue = gcd_row(A.row(i), b.coeffRef(i));                                       //compute the gcd to make the values as small as possible

      if(gcdValue == 1)
        continue;

      for(Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(A, i); it; ++it)
      {
        it.valueRef() = (it.value() / gcdValue);
      }
      b.coeffRef(i) = b.coeffRef(i) / gcdValue;
    }
  }
}

void ExactConstraintSatisfaction::IRREF_Jordan(SP_Matrix_R& A, Eigen::VectorXi &b)
{
  for(int k = number_pivots_ - 1; k > 0; k--){
    int pivot_col = -1;
    for(Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(A, k); it; ++it){
      if(it.value() != 0){
        pivot_col = it.index();
        break;
      }
    }
    if(pivot_col == -1)
      std::cout << "Error in IRREF_Jordan(ExactConstraintSatisfaction.cc) : couldn't find a pivot_col." << std::endl;

    for(int i = k-1; i >= 0; i--){                                                                      //eliminate row i with row k

      if(A.coeff(i, pivot_col) == 0)
        continue;
      if(b.coeff(i) != 0 || b.coeff(k) != 0)
        b.coeffRef(i) = (A.coeff(k, pivot_col) * b.coeff(i) - A.coeff(i, pivot_col) * b.coeff(k));        //do it in this order, so we don't delete entries in A for the computation of b
      eliminate_row(A, i, k, pivot_col);

      int gcdValue = gcd_row(A.row(i), b.coeffRef(i));                                                   //compute the gcd to make the values as small as possible
      if(gcdValue == 1)
        continue;

      for(Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(A, i); it; ++it){
        it.valueRef() = (it.value() / gcdValue);
      }
      b.coeffRef(i) = b.coeffRef(i) / gcdValue;
    }
  }
  A.makeCompressed();
  A.finalize();
  A.prune(0.0, 0);
}

//-------------------------------------Evaluation--------------------------------------------------------------------

void ExactConstraintSatisfaction::evaluation(SP_Matrix_R& _A, Eigen::VectorXi& b, Eigen::VectorXd& x, const Eigen::VectorXd values)
{
  //debug
  double time_G = 0.0;
  double time_J = 0.0;
  double time_e = 0.0;
  double time_count;
  //debug

  time_count = clock();
  IREF_Gaussian(_A, b, x);
  time_G = clock() - time_count;
  time_G = time_G/CLOCKS_PER_SEC;
  _A.prune(0.0 , 0);
  _A.makeCompressed();
  _A.finalize();
  std::cout << "error is after Gaus: " << (_A.cast<double>() * x).squaredNorm() << std::endl;

  //debug
  time_count = clock();
  //debug

  IRREF_Jordan(_A, b);

  //debug
  time_J = clock() - time_count;
  time_J = time_J/CLOCKS_PER_SEC;
  //debug

  _A.prune(0.0 , 0);
  _A.makeCompressed();
  _A.finalize();
  std::cout << "error is after Jordan: " << (_A.cast<double>() * x).squaredNorm() << std::endl;

  //debug
  time_count = clock();
  //debug


  SP_Matrix_C A = _A;         //change the matrix type to allow easier iteration


  int cols = A.cols();
  std::cout << "largest Expo" << std::endl;
  largest_exponent(values);
  std::cout << "findPivo" << std::endl;
  for(int k = cols -1; k >= 0; k--)
  {
    auto pivot_row = get_pivot_row_new(A, _A, k);

    if(pivot_row == -1)
    {
      //there is no pivot in this column
      auto D = get_divisors_new(A, _A, k);
      x.coeffRef(k) = makeDiv(D, x.coeffRef(k));            //fix free variables so they are in F_delta
    }
    else
    {

      auto S = get_dot_product_elements_new(_A, x, pivot_row);

      double divided_B = b.coeffRef(pivot_row);
      divided_B = F_delta(divided_B / A.coeffRef(pivot_row, k));
      if( divided_B * A.coeffRef(pivot_row, k) !=  static_cast<double>(b.coeffRef(pivot_row)))
        std::cout << "WARNING: Can't handle the right hand side perfectly" << std::endl;
      x.coeffRef(k) = divided_B - safeDot(S);

    }
  }

  //debug
  time_e = clock() - time_count;
  time_e = time_e/CLOCKS_PER_SEC;

  std::cout.precision(64);
  std::cout << "time for IREF: " << time_G << " Time for IRREF: " << time_J<< " Time for evaluation: " << time_e << std::endl;
  //debug

}

double ExactConstraintSatisfaction::makeDiv(const std::vector<int>& D, double x)
{
  if(D.empty()){
    return F_delta(x);
  }
  int d = lcm(D);
  double result = F_delta(x/d) * d;
  return result;

}

double ExactConstraintSatisfaction::safeDot(const std::vector<std::pair<int, double> >& S)
{
  if (S.empty()) //case all Cij are zero after the pivot
    return 0;
  int safebreak = 9999999;
  std::vector<std::pair<int, double>> P;
  std::vector<std::pair<int, double>> N;

  int k = 0;
  double r = 0;                                         //return value of the dot

  for(auto element : S){
    if(element.first * element.second > 0)
    {
      element.first = std::abs(element.first);
      element.second = std::abs(element.second);
      P.push_back(element);
    }
    else if(element.first * element.second < 0)
    {
      element.first = std::abs(element.first);
      element.second = - std::abs(element.second);
      N.push_back(element);
    }
  }

  while((!P.empty() || !N.empty()) && safebreak > 0)
  {
    safebreak--;                                        //break out of the while after an amount of time
    if(!P.empty() && (r < 0 || N.empty()))
    {
      const std::pair<int, double> element = P.back();
      P.pop_back();
      double test_value = element.second;
      if(test_value < 0.00000001)
      { //to prevent overflow through the dividing
        k = element.first;
      }
      else
      {
        k = std::min(element.first, static_cast<int>(std::floor((delta_ - r)/element.second)));
      }
      r = r + k * element.second;
      if(k < element.first)
        P.push_back({element.first - k, element.second});
    }
    else
    {
      const std::pair<int, double> element = N.back();
      N.pop_back();
      double test_value = element.second;
      if(std::abs(test_value) < 0.00000001)
      {            //to prevent overflow through the dividing
        k = element.first;
      }
      else
      {
        k = std::min(element.first, static_cast<int>(std::floor((-delta_ - r)/element.second)));
      }
      r = r + k * element.second;
      if(k < element.first)
        N.push_back({element.first - k, element.second});
    }

    if(k == 0)
    {
      std::cout << "ERROR: The representable range of double values (delta) is to small (ExactConstraintSatisfaction.cc)" << std::endl;
      return -99999999;
    }
  }
  if(safebreak == 0)
    std::cout << "WARNING:: the set number of Iterations in Method : safedot is exceeded" << std::endl;

  return r;
}

int ExactConstraintSatisfaction::get_pivot_row_student(const SP_Matrix_C& A, int col)
{
  int pivot_row = -1;
  for(SP_Matrix_C::InnerIterator it(A, col); it; ++it)
  {
    if(it.value() != 0)
    {
      int index = index_pivot(A.row(it.index()));
      if(index == col)
      {
        pivot_row = it.index();
        break;
      }
    }
    else
    {
      COMISO_THROW_TODO("There should be no non zero values in the matrix");
    }
  }
  return pivot_row;
}

int ExactConstraintSatisfaction::get_pivot_row_new(const SP_Matrix_C& A, const SP_Matrix_R& _A, int col)
{
  auto collumn = A.col(col);
  if (collumn.nonZeros() != 1) // a pivot is allways the only entry in a column
    return -1;

  auto row = SP_Matrix_C::InnerIterator(A, col).index();

  // check if col is the first element in row
  auto first_index = SP_Matrix_R::InnerIterator(_A, row).index();

  if (first_index == col)
    return row;

  return -1;
}

std::vector<int> ExactConstraintSatisfaction::get_divisors_student(const ExactConstraintSatisfaction::SP_Matrix_C& A, int col)
{
  std::vector<int> D;
  for(SP_Matrix_C::InnerIterator it(A, col); it; ++it)
  {
    COMISO_THROW_TODO_if(it.value() == 0, "There should be no zeros left in the matrix");
    if(it.value() != 0 && it.index() <= col && it.index() < number_pivots_)
    {
      int pivot_col = index_pivot(A.row(it.index()));
      D.push_back(A.coeff(it.index(), pivot_col));
    }
  }

  return D;
}

std::vector<int> ExactConstraintSatisfaction::get_divisors_new(const SP_Matrix_C& A, const SP_Matrix_R& _A, int col)
{
  std::vector<int> D;
  for(SP_Matrix_C::InnerIterator it(A, col); it; ++it)
  {
    COMISO_THROW_TODO_if(it.value() == 0,              "There should be no zeros left in the matrix");
    if (it.index() >= number_pivots_)
      std::cout << A << std::endl;
    COMISO_THROW_TODO_if(it.index() >= number_pivots_, "The matrix should only contain number_pivots non empty rows");
    COMISO_THROW_TODO_if(it.index() > col,             "The matrix should not contain elements below the diagonal");

    D.push_back(SP_Matrix_R::InnerIterator(_A, it.index()).value());
  }
  return D;
}

std::vector<std::pair<int, double> > ExactConstraintSatisfaction::get_dot_product_elements_student(const SP_Matrix_C& A, const Eigen::VectorXd& x,  int k, int pivot_row)
{
  std::vector<std::pair<int, double>> S;
  int cols = A.cols();
  for(int i = k+1; i < cols; i++)
  {                      //construct the list S to do the dot Product
    std::pair<int, double> pair;
    pair.first = A.coeff(pivot_row,i);
    double test = x.coeff(i) / A.coeff(pivot_row,k);
    if(x.coeff(i) != ( test * A.coeff(pivot_row,k)))
      std::cout << "WARNING: can't devide" << " in row : " << i << std::endl;
    pair.second = x.coeff(i) / A.coeff(pivot_row,k);
    if(pair.first != 0)
      S.push_back(pair);
  }

  return S;
}

std::vector<std::pair<int, double> > ExactConstraintSatisfaction::get_dot_product_elements_new(const SP_Matrix_R& A, const Eigen::VectorXd& x, int pivot_row)
{
  std::vector<std::pair<int, double>> S;

  SP_Matrix_R::InnerIterator it(A, pivot_row);
  auto pivot_val = it.value();

  while (++it)
  {
    std::pair<int, double> pair;
    pair.first = it.value();
    auto tmp = x.coeff(it.index());
    COMISO_THROW_TODO_if((tmp / pivot_val) * pivot_val != tmp, "element in x is not divisible by pivot element");
    pair.second = tmp / pivot_val;
    S.push_back(pair);
  }
  return S;
}

