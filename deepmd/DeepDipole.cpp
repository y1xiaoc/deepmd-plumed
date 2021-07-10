#include "core/PlumedMain.h"
#include "Common.h"
#include "DeepModelBase.h"
#include "deepmd/DeepTensor.h"

namespace PLMD {
namespace dp_plmd { // to avoid conflicts with libdeepmd

constexpr unsigned int odim = 3;
const std::array<std::string, odim> cpnts = {"x", "y", "z"};

//+PLUMEDOC COLVAR DEEPDIPOLE
/*
Calculate the dipole moment vector for the system using a deep tensor model.

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
dipole: DEEPDIPOLE MODEL=dipole.pb ATYPE=type.raw UNIT_CVT=1.0
\endplumedfile
*/
//+ENDPLUMEDOC
class DeepDipole : public DeepModelBase<deepmd::DeepTensor, odim> {
public:
  explicit DeepDipole(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(DeepDipole,"DEEPDIPOLE")

void DeepDipole::registerKeywords(Keywords& keys) {
  DeepModelBase<deepmd::DeepTensor, odim>::registerKeywords(keys);
  for (unsigned k = 0; k < odim; ++k) {
    keys.addOutputComponent(cpnts[k],"COMPONENTS","the "+cpnts[k]+"-component of the dipole");
  }
}

DeepDipole::DeepDipole(const ActionOptions&ao):
  Action(ao),
  DeepModelBase(ao)
{  
  if (output_unit < 0) {
    output_unit = global_dipole_unit;
  }
  log.printf("  output unit conversion set to %f\n", output_unit);

  for (unsigned k = 0; k < odim; ++k) {
    addComponentWithDerivatives(cpnts[k]); componentIsNotPeriodic(cpnts[k]);
  }
  
  if (dp.output_dim() != odim) { 
    throw std::runtime_error( "invalid graph file! the output dimension should be 3" ); 
  }

  requestAtoms(atoms); 
}


}
}
