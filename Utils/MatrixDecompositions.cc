#include <CoMISo/Utils/MatrixDecompositions.hh>
#include <CoMISo/Utils/MatrixDecompositionsT_impl.cc>


namespace COMISO {

template std::unique_ptr<MatrixDecomposition<double>> make_decomposition(MatrixDecompositionAlgorithm);

} // namespace COMISO
