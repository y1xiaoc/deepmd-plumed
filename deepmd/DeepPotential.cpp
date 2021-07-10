#include "core/PlumedMain.h"
#include "Common.h"
#include "DeepModelBase.h"
#include "deepmd/DeepPot.h"

namespace PLMD {
namespace dp_plmd { // to avoid conflicts with libdeepmd
//+PLUMEDOC COLVAR DEEPPOTENTIAL
/*
Calculate the potential energy for the system using a deep potential model.

A binary graph file of the model and a text file specify the type of each atom are needed.
The type file should be a list of integers with the i-th element correspond to the type 
in the deep tensor model of the i-th atom in the system (or specified by you).
The output is scaled by a factor given by UNIT_CVT, which defaults to 1.
By default will use periodic boundary conditions, which will be handled by 
the deep tensor model automatically. In case NOPBC flag is specified, the box
will be ignored and there will be no the pbc handling.

\par Examples

Here's a simple example showing how to use this CV (also the default values of the keywords):
\plumedfile
dp: DEEPPOTENTIAL MODEL=graph.pb ATYPE=type.raw UNIT_CVT=96.487
\endplumedfile
*/
//+ENDPLUMEDOC
class DeepPotential : public DeepModelBase<deepmd::DeepPot, 1> {
public:
  explicit DeepPotential(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(DeepPotential,"DEEPPOTENTIAL")

void DeepPotential::registerKeywords(Keywords& keys) {
  DeepModelBase<deepmd::DeepPot, 1>::registerKeywords(keys);
}

DeepPotential::DeepPotential(const ActionOptions&ao):
  Action(ao),
  DeepModelBase(ao)
{
  if (output_unit < 0) {
    output_unit = global_energy_unit / plumed.getAtoms().getUnits().getEnergy();
  }
  log.printf("  output unit conversion set to %f\n", output_unit);

  addValueWithDerivatives(); setNotPeriodic();
  
  requestAtoms(atoms); 
}


}
}
