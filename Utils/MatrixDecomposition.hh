#pragma once

#include <memory>
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace COMISO {

enum class MatrixDecompositionAlgorithm {
    // TODO: split ordering into separate parameter?
    Eigen_SimplicialLDLD_MetisOrdering,
    Eigen_SimplicialLDLT_NaturalOrdering,
    Eigen_SimplicialLDLT_AMDOrdering,
    Eigen_SimplicialCholesky,
    // TODO LLT with eigen?
    Cholmod_Supernodal,
    Cholmod_SimplicialLLT,
    UmfPack_LU,
};


template<typename Scalar=double>
class MatrixDecomposition {
public:
    using SMatrix = Eigen::SparseMatrix<Scalar>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    virtual ~MatrixDecomposition() = default;

    virtual void analyzePattern(const SMatrix& a) = 0;
    virtual void compute(const SMatrix& a) = 0;
    virtual void factorize(const SMatrix& a) = 0;

    virtual Eigen::ComputationInfo info() const = 0;
    virtual Vector solve(Vector const&) const = 0;
};

template<typename Scalar=double>
std::unique_ptr<MatrixDecomposition<Scalar>> make_decomposition(MatrixDecompositionAlgorithm);

extern template
std::unique_ptr<MatrixDecomposition<double>> make_decomposition(MatrixDecompositionAlgorithm);

} // namespace COMISO
