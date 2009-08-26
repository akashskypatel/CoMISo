//=============================================================================
//
//  CLASS HarmonicExampleT
//
//=============================================================================


#ifndef ACG_HARMONICEXAMPLET_HH
#define ACG_HARMONICEXAMPLET_HH


//== INCLUDES =================================================================
#include <gmm/gmm.h>

#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/System/omstream.hh>

#include <vector>
#include <string>
#include <iostream>
#include <cfloat> //for FLT_{MIN, MAX}
#include <CoMISo/Utils/StopWatch.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <CoMISo/Solver/GMM_Tools.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <CoMISo/Solver/ConstrainedSolver.hh>
#include <CoMISo/Utils/ColorCoder.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================


	      
/** \class HarmonicExampleT HarmonicExampleT.hh <ACG/.../HarmonicExampleT.hh>

    Brief Description.
  
    A more elaborate description follows.
*/

template <class MeshT>
class HarmonicExampleT
{
public:

  typedef typename MeshT::Scalar         Scalar;
	typedef typename MeshT::Point          Point;
	typedef typename OpenMesh::Vec2d      Vec2d;
	typedef typename MeshT::VertexHandle   VH;
	typedef typename MeshT::FaceHandle     FH;
	typedef typename MeshT::EdgeHandle     EH;
	typedef typename MeshT::HalfedgeHandle HEH;
	typedef typename MeshT::FaceIter				FIter;
	typedef typename MeshT::EdgeIter				EIter;
	typedef typename MeshT::VertexIter			VIter;
	typedef typename MeshT::HalfedgeIter		HIter;
	typedef typename MeshT::FaceFaceIter	FFIter;
	typedef typename MeshT::FaceVertexIter	FVIter;
	typedef typename MeshT::VertexFaceIter	VFIter;
	typedef typename MeshT::VertexEdgeIter	VEIter;
	typedef typename MeshT::FaceHalfedgeIter FHIter;
	typedef typename MeshT::VertexOHalfedgeIter VOHIter;
	typedef typename MeshT::VertexIHalfedgeIter VIHIter;
	typedef typename MeshT::TexCoord2D			TexCoord2D;
	typedef unsigned int uint;
	typedef typename gmm::linalg_traits< gmm::rsvector< Scalar > > GMMRSTRAIT;
	typedef typename gmm::linalg_traits< gmm::wsvector< Scalar > > GMMWSTRAIT;
	typedef typename GMMRSTRAIT::iterator GMMColIter;
	typedef typename GMMWSTRAIT::iterator wGMMColIter;

  typedef typename gmm::wsvector< double > SVectorT;
  typedef gmm::linalg_traits<SVectorT>::const_iterator SVIter;
  typedef typename gmm::col_matrix< gmm::wsvector< Scalar > > CMatrixT;
  typedef typename gmm::row_matrix< gmm::wsvector< Scalar > > RMatrixT;
  typedef typename std::vector< Scalar > VectorT;
  typedef typename std::vector< int >    VectorIT;
  typedef gmm::dense_matrix<Scalar>      DenseMatrix;
  typedef gmm::csc_matrix<double>        CSCMatrix;
  //  typedef ACG::MCSCMatrix                MCSCMatrix;

 

  /// Constructor
  HarmonicExampleT( MeshT& _mesh) : mesh_(_mesh)
  {}
 
  /// Destructor
  ~HarmonicExampleT() {}

  void compute();

private:

  // create a simple connectivity-based laplace matrix
  void init_laplace_matrix();

  int  add_laplace_row( VH _vh, int _row);
  void add_integer_constraints();

  /// adds some constraints to the laplace system.
  /// The vertices of each boundary are constrained to the same value. 
  /// The different boundaries are constrained to differ by e.g. 5
  int  add_boundary_constraints();
  
  void color_mesh();

  // laplace matrix 
  CMatrixT laplace_;

  // solution
  std::vector< double > x_;
  
  // right hand side
  std::vector< double > rhs_;

  // constraints
  RMatrixT constraints_;

  // integer constraints (indices of variables to be rounded)
  std::vector< int > ids_to_round_;

  MeshT& mesh_;
};


//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_HARMONICEXAMPLET_C)
#define ACG_HARMONICEXAMPLE_TEMPLATES
#include "HarmonicExampleT.cc"
#endif
//=============================================================================
#endif // ACG_HARMONICEXAMPLET_HH defined
//=============================================================================

