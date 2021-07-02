/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)
   See http://www.plumed.org for more information.
   This file is part of plumed, version 2.
   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionAtomistic.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"

#include "Common.h"
#include "deepmd/DeepTensor.h"

namespace PLMD {
namespace dp_plmd { // to avoid conflicts with libdeepmd
//+PLUMEDOC COLVAR DEEPPOLAR
/*
Calculate the polarizability tensor for the system using a deep tensor model.

A binary graph file of the model and a text file specify the type of each atom are needed.
The type file should be a list of integers with the i-th element correspond to the type 
in the deep tensor model of the i-th atom in the system (or specified by you).
The output is scaled by a factor given by UNIT_CVT, which defaults to 1.
By default will use periodic boundary conditions, which will be handled by 
the deep tensor model automatically. In case NOPBC flag is specified, the box
will be enlarged to avoid the pbc handling. That should be used in a non pbc system.

\par Examples

Here's a simple example showing how to use this CV (also the default values of the keywords):
\plumedfile
polar: DEEPPOLAR MODEL=polar.pb ATYPE=type.raw UNIT_CVT=1.0
\endplumedfile
*/
//+ENDPLUMEDOC
class DeepPolar : public Colvar {
  std::vector<AtomNumber> atoms; 
  bool nopbc;
public:
  explicit DeepPolar(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
private:
  deepmd::DeepTensor dp; 
  std::vector<int> atype;
  double polar_unit;
  double length_unit;
  constexpr static int odim = 9;
  static const std::array<std::string, odim> cpnts;
};

PLUMED_REGISTER_ACTION(DeepPolar,"DEEPPOLAR")

void DeepPolar::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","the group of atoms we are calculating the polarizability for (defaults to the whole system)");
  keys.add("compulsory","MODEL","polar.pb","the DeepPolar model binary graph file");
  keys.add("compulsory","ATYPE","type.raw" ,"the file specify the type (in the model) of each atom");
  keys.add("optional","UNIT_CVT","the unit conversion constant of output polar tensor (will be multiplied to the graph output, default is 1.0)");
  keys.addOutputComponent("xx","COMPONENTS","the xx-component of the polarizability tensor");
  keys.addOutputComponent("xy","COMPONENTS","the xy-component of the polarizability tensor");
  keys.addOutputComponent("xz","COMPONENTS","the xz-component of the polarizability tensor");
  keys.addOutputComponent("yx","COMPONENTS","the yx-component of the polarizability tensor");
  keys.addOutputComponent("yy","COMPONENTS","the yy-component of the polarizability tensor");
  keys.addOutputComponent("yz","COMPONENTS","the yz-component of the polarizability tensor");
  keys.addOutputComponent("zx","COMPONENTS","the zx-component of the polarizability tensor");
  keys.addOutputComponent("zy","COMPONENTS","the zy-component of the polarizability tensor");
  keys.addOutputComponent("zz","COMPONENTS","the zz-component of the polarizability tensor");
}

const std::array<std::string, DeepPolar::odim> DeepPolar::cpnts = {
  "xx", "xy", "xz", 
  "yx", "yy", "yz", 
  "zx", "zy", "zz", 
};

DeepPolar::DeepPolar(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nopbc(false),
  polar_unit(global_polar_unit),
  length_unit(global_length_unit)
{
  parseAtomList("ATOMS",atoms);
  parseFlag("NOPBC",nopbc);
  std::string graph_file;
  parse("MODEL", graph_file);
  std::string type_file;
  parse("ATYPE", type_file);
  parse("UNIT_CVT", polar_unit);

  checkRead();
  
  for (unsigned k = 0; k < odim; ++k) {
    addComponentWithDerivatives(cpnts[k]); componentIsNotPeriodic(cpnts[k]);
  }

  // make sure the length unit passed to graph is Angstrom
  length_unit /= plumed.getAtoms().getUnits().getLength();

  // default use all atoms; otherwise warn the user
  if (atoms.size() == 0) {
    atoms.resize(getTotAtoms());
    for(unsigned i = 0; i < getTotAtoms(); ++i) {
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
  if (atoms.size() != getTotAtoms()) { 
    log.printf("  # of atoms provided: %d != # of atoms in the system: %d\n", atoms.size(), getTotAtoms()); 
    log.printf("  Please make sure you know what you are doing!\n");
  }

  if(nopbc) { log.printf("  without periodic boundary conditions\n"); }
  else      { log.printf("  using periodic boundary conditions\n"); }

  // load deepmd model from graph file; check dimention is correct
  log.printf("  using graph file:  %s \n", graph_file.c_str());
  dp.init(graph_file);
  if (dp.output_dim() != odim) { 
    throw std::runtime_error( "invalid graph file! the output dimension should be 3" ); 
  }
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

  // the output unit will be multiple to the graph output
  log.printf("  output unit conversion set to %f\n", polar_unit);

  requestAtoms(atoms); 
}

void DeepPolar::calculate()
{
  if (!nopbc) { makeWhole(); }
  unsigned N = getNumberOfAtoms();
  std::vector<FLOAT_PREC> _polar(odim);
  std::vector<FLOAT_PREC> _force (odim * N * 3);
  std::vector<FLOAT_PREC> _virial(odim * 9);
  std::vector<FLOAT_PREC> _coord (N * 3);
  std::vector<FLOAT_PREC> _box   (9); 
  IndexConverter ic(odim, N, 3);
  
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

  dp.compute(_polar, _force, _virial, _coord, atype, _box);

  // get back polar
  for (unsigned k = 0; k < odim; ++k) {
    getPntrToComponent(cpnts[k])->set(_polar[k] * polar_unit);
  }
  // get back force
  for (unsigned k = 0; k < odim; ++k) {
    for(unsigned i = 0; i < N; i++) {
      setAtomsDerivatives(
        getPntrToComponent(cpnts[k]), 
        i,
        - Vector(_force[ic.f(k, i, 0)],
                 _force[ic.f(k, i, 1)], 
                 _force[ic.f(k, i, 2)]) 
        * polar_unit / length_unit
      );
    }
  }
  // get back virial
  for (unsigned k = 0; k < odim; ++k) {
    // setBoxDerivativesNoPbc(getPntrToComponent(cpnts[k]));
    // ^^^^^^^^ will ask plumed to calculate virial ^^^^^^^^
    setBoxDerivatives(
      getPntrToComponent(cpnts[k]),
      Tensor( // we manually transpose here; the dp convention is different
        _virial[k*9 + 0], _virial[k*9 + 3], _virial[k*9 + 6], 
        _virial[k*9 + 1], _virial[k*9 + 4], _virial[k*9 + 7], 
        _virial[k*9 + 2], _virial[k*9 + 5], _virial[k*9 + 8]) 
    );
  }
}

}
}
