//=============================================================================
//
// CLASS HarmonicExamplePerObjectData
//
//=============================================================================


#ifndef HARMONICEXAMPLEPEROBJECTDATA_HH
#define HARMONICEXAMPLEPEROBJECTDATA_HH

//== INCLUDES =================================================================


#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenFlipper/common/perObjectData.hh>
#include "HarmonicExampleT.hh"


//== CLASS DEFINITION =========================================================


template <class MeshT>
class HarmonicExamplePerObjectDataT : public PerObjectData
{
public:

  HarmonicExamplePerObjectDataT( MeshT& _mesh) : harmonicexample_(_mesh) 
  {}
  
  virtual
  ~HarmonicExamplePerObjectDataT() 
  {}
  
  ACG::HarmonicExampleT<MeshT>& harmonicexample() { return harmonicexample_;}

private:
  // create an FeatureLine
  ACG::HarmonicExampleT<MeshT> harmonicexample_;
};


//=============================================================================
#endif // HARMONICEXAMPLEPEROBJECTDATA_HH defined
//=============================================================================

