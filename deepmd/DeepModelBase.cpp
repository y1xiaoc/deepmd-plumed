#include "core/ActionAtomistic.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include "DeepModelBase.h"
#include "Common.h"

#include "deepmd/DeepTensor.h"
#include "deepmd/DeepPot.h"

namespace PLMD {
namespace dp_plmd { // to avoid conflicts with libdeepmd

template <class DP, unsigned int ODIM>
void
DeepModelBase<DP, ODIM>::
registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","the group of atoms we are calculating the CV for (defaults to the whole system)");
  keys.add("compulsory","MODEL","graph.pb","the deep model binary graph file");
  keys.add("compulsory","ATYPE","type.raw" ,"the file specify the type (in the model) of each atom");
  keys.add("optional","UNIT_CVT","the unit conversion constant of model output (will be multiplied to the graph output, default is converting to plumed energy unit for potential and 1.0 for tensor)");
}

template <class DP, unsigned int ODIM>
DeepModelBase<DP, ODIM>::
DeepModelBase(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nopbc(false)
{
  // make sure the length unit passed to graph is Angstrom
  length_unit = global_length_unit / plumed.getAtoms().getUnits().getLength();
  // try to use the default unit. but user can change this
  output_unit = -1;

  parseAtomList("ATOMS",atoms);
  parseFlag("NOPBC",nopbc);
  std::string graph_file;
  parse("MODEL", graph_file);
  std::string type_file;
  parse("ATYPE", type_file);
  parse("UNIT_CVT", output_unit);

  checkRead();
  
  // default use all atoms; otherwise warn the user
  unsigned nmdatoms = plumed.getAtoms().getNatoms();
  if (atoms.size() == 0) {
    atoms.resize(nmdatoms);
    for(unsigned i = 0; i < nmdatoms; ++i) {
      atoms[i].setIndex(i);
    }
    log.printf("  of all %u atoms in the system\n",static_cast<unsigned>(atoms.size()));
  } else {
    log.printf("  of %u atoms\n",static_cast<unsigned>(atoms.size()));
  }
  for(unsigned i = 0; i < atoms.size(); ++i) {
    log.printf("  %d", atoms[i].serial());
  }
  log.printf("  \n");
  if (atoms.size() != nmdatoms) { 
    log.printf("  # of atoms provided: %d != # of atoms in the system: %d\n", atoms.size(), nmdatoms); 
    log.printf("  Please make sure you know what you are doing!\n");
  }

  if(nopbc) { log.printf("  without periodic boundary conditions\n"); }
  else      { log.printf("  using periodic boundary conditions\n"); }

  // load deepmd model from graph file; check dimention is correct
  log.printf("  using graph file:  %s \n", graph_file.c_str());
  dp.init(graph_file);
  // dp.print_summary("  DeePMD: ");
  log.printf("  DeepTensor model initialized successfully\n");

  // load atom types that will be feed into the graph
  load_atype(atype, type_file);
  if (atype.size() != atoms.size()) { 
    throw std::runtime_error( "invalid atom type file! the size should be equal to the number of atoms" ); 
  }
  log.printf("  assign type to atoms\n");
  for(unsigned i = 0; i < atype.size(); ++i) {
    log.printf("  %d", atype[i]);
  }
  log.printf("  \n");
}

template <class DP, unsigned int ODIM>
void
DeepModelBase<DP, ODIM>::
calculate() {
  if (!nopbc) { makeWhole(); }
  unsigned N = getNumberOfAtoms();
  std::vector<FLOAT_PREC> _output(ODIM);
  std::vector<FLOAT_PREC> _force (ODIM * N * 3);
  std::vector<FLOAT_PREC> _virial(ODIM * 9);
  std::vector<FLOAT_PREC> _coord (N * 3);
  std::vector<FLOAT_PREC> _box   (9); 
  IndexConverter ic(ODIM, N, 3);
  
  // copy atom coords
  for (unsigned i = 0; i < N; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      _coord[ic.f(i,j)] = getPosition(i)[j] / length_unit;
    }
  }
  // copy box tensor
  Tensor box = getBox();
  if (nopbc) { 
    _box.clear(); // if size is not 9 then no pbc
  } else {
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        _box[i * 3 + j] = box(i, j) / length_unit;
      }
    }
  }

  dp_compute(_output, _force, _virial, _coord, atype, _box);

  // get back output
  for (unsigned k = 0; k < ODIM; ++k) {
    getPntrToComponent(k)->set(_output[k] * output_unit);
  }
  // get back force
  for (unsigned k = 0; k < ODIM; ++k) {
    for(unsigned i = 0; i < N; i++) {
      setAtomsDerivatives(
        getPntrToComponent(k), 
        i,
        - Vector(_force[ic.f(k, i, 0)],
                 _force[ic.f(k, i, 1)], 
                 _force[ic.f(k, i, 2)]) 
        * output_unit / length_unit
      );
    }
  }
  // get back virial
  for (unsigned k = 0; k < ODIM; ++k) {
    // setBoxDerivativesNoPbc(getPntrToComponent(k));
    // ^^^^^^^^ will ask plumed to calculate virial ^^^^^^^^
    setBoxDerivatives(
      getPntrToComponent(k),
      Tensor( // we manually transpose here; the dp convention is different
        _virial[k*9 + 0], _virial[k*9 + 3], _virial[k*9 + 6], 
        _virial[k*9 + 1], _virial[k*9 + 4], _virial[k*9 + 7], 
        _virial[k*9 + 2], _virial[k*9 + 5], _virial[k*9 + 8]) 
      * output_unit
    );
  }
}

template <class DP, unsigned int ODIM>
void
DeepModelBase<DP, ODIM>::
dp_compute (std::vector<FLOAT_PREC> & output,
            std::vector<FLOAT_PREC> & force,
            std::vector<FLOAT_PREC> & virial,
            const std::vector<FLOAT_PREC> & coord,
            const std::vector<int> & atype,
            const std::vector<FLOAT_PREC> & box) 
{
  dp.compute(output, force, virial, coord, atype, box);
}

template<>
void
DeepModelBase<deepmd::DeepPot, 1>::
dp_compute (std::vector<FLOAT_PREC> & output,
            std::vector<FLOAT_PREC> & force,
            std::vector<FLOAT_PREC> & virial,
            const std::vector<FLOAT_PREC> & coord,
            const std::vector<int> & atype,
            const std::vector<FLOAT_PREC> & box) 
{
  double _output;
  dp.compute(_output, force, virial, coord, atype, box);
  output.clear();
  output.push_back(_output);
}

template class DeepModelBase<deepmd::DeepPot, 1>;
template class DeepModelBase<deepmd::DeepTensor, 3>;
template class DeepModelBase<deepmd::DeepTensor, 9>;

}
}