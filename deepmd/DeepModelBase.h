#pragma once

#include "Common.h"
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"

namespace PLMD {
namespace dp_plmd { // to avoid conflicts with libdeepmd

template <class DP, unsigned int ODIM>
class DeepModelBase : public Colvar {
public:
  explicit DeepModelBase (const ActionOptions&);
  void calculate () override;
  static void registerKeywords (Keywords& keys);
protected:
  std::vector<AtomNumber> atoms; 
  bool nopbc;
  // the deepmd model class, can be DeepPot or DeepTensor 
  DP dp; 
  // the type list of each atoms
  std::vector<int> atype;
  // the conversion const for the model output
  double output_unit;
  // the conversion const for the model input (atom position)
  double length_unit;
private:
  void dp_compute (std::vector<FLOAT_PREC> & output,
                   std::vector<FLOAT_PREC> & force,
                   std::vector<FLOAT_PREC> & virial,
                   const std::vector<FLOAT_PREC> & coord,
                   const std::vector<int> & atype,
                   const std::vector<FLOAT_PREC> & box);
};

}
}