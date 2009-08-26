#include "ui_HarmonicExampleToolbarBase.hh"
#include <QtGui>

class HarmonicExampleToolbar : public QWidget, public Ui::HarmonicExampleToolbarBase
{
  Q_OBJECT

public:
  HarmonicExampleToolbar(QWidget * parent = 0)
    : QWidget(parent)
  {
    setupUi(this);
  }
};
