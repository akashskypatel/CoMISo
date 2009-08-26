//=============================================================================
//
//  CLASS HarmonicExamplePlugin
//
//=============================================================================


#ifndef HARMONICEXAMPLEPLUGIN_HH
#define HARMONICEXAMPLEPLUGIN_HH


//== INCLUDES =================================================================

#include <QObject>
#include <QMenuBar>
#include <QSpinBox>

#include <OpenFlipper/common/Types.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/KeyInterface.hh>
#include <OpenFlipper/BasePlugin/MouseInterface.hh>
#include <OpenFlipper/BasePlugin/PickingInterface.hh>
#include <OpenFlipper/BasePlugin/ScriptInterface.hh>

#include <ACG/QtWidgets/QtExaminerViewer.hh>

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include "HarmonicExampleToolbar.hh"
#include "HarmonicExamplePerObjectDataT.hh"
#include "HarmonicExampleT.hh"


//== CLASS DEFINITION =========================================================


class HarmonicExamplePlugin : public QObject, BaseInterface, ToolboxInterface, KeyInterface, ScriptInterface, MouseInterface, PickingInterface
{
  Q_OBJECT
  Q_INTERFACES(BaseInterface)
  Q_INTERFACES(ToolboxInterface)
  Q_INTERFACES(KeyInterface)
  Q_INTERFACES(ScriptInterface)
  Q_INTERFACES(MouseInterface)
  Q_INTERFACES(PickingInterface)


  // typedef for easy access
  typedef ACG::HarmonicExampleT<TriMesh> HarmonicExample;
  typedef HarmonicExamplePerObjectDataT<TriMesh>   POD;

signals:
  void updateView();
  void updatedObject(int);


private slots:

  // initialization functions
  void initializePlugin();
  void pluginsInitialized();


  // compute
  void slotCompute();


public :

  ~HarmonicExamplePlugin() {};


  bool initializeToolbox(QWidget*& _widget);

  QString name() { return (QString("HarmonicExample")); };
  QString description( ) { return (QString("Computes the HarmonicExample of the the active Mesh")); }; 

private :

  // return name of per object data
  const char * pod_name() { return "HARMONICEXAMPLE_PER_OBJECT_DATA";}

  // get HarmonicExample object for a given object
  HarmonicExample* get_harmonicexample_object( BaseObjectData* _object );
  


private :
  /// Widget for Toolbox
  HarmonicExampleToolbar* tool_;
};


//=============================================================================
#endif // HARMONICEXAMPLEPLUGIN_HH defined
//=============================================================================

