#include <CoMISo/Utils/MatrixDecomposition.hh>
#include <CoMISo/Utils/MatrixDecompositionT_impl.hh>


namespace COMISO {

template std::unique_ptr<MatrixDecomposition<double>> make_decomposition(MatrixDecompositionAlgorithm);

} // namespace COMISO
