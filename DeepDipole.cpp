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
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"

#include <string>
#include <cmath>

#include "deepmd/DeepTensor.h"

using namespace std;
// using namespace tensorflow;

namespace PLMD {
namespace deepmd_plmd { // to avoid conflicts with libdeepmd
//+PLUMEDOC COLVAR DEEPDIPOLE
/*
Calculate the dipole moment for a group of atoms.
When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding the molecule with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.
In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.
*/
//+ENDPLUMEDOC
class DeepDipole : public Colvar {
  vector<AtomNumber> atoms; 
  bool nopbc;
public:
  explicit DeepDipole(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
private:
  deepmd::DeepTensor dp;    
};

PLUMED_REGISTER_ACTION(DeepDipole,"DEEPDIPOLE")

void DeepDipole::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","the group of atoms we are calculating the dipole moment for"); // TODO assert GROUP: @allatoms
  keys.add("compulsory","MODEL","dipole.pb","the DeepDipole model file");
  keys.addOutputComponent("x","COMPONENTS","the x-component of the dipole");
  keys.addOutputComponent("y","COMPONENTS","the y-component of the dipole");
  keys.addOutputComponent("z","COMPONENTS","the z-component of the dipole");
}

DeepDipole::DeepDipole(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  nopbc(false)
{
  parseAtomList("ATOMS",atoms); // atoms list
  parseFlag("NOPBC",nopbc);
  string graph_file;
  parse("MODEL", graph_file);

  checkRead();
  
  addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
  addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");

  log.printf("  of %u atoms\n",static_cast<unsigned>(atoms.size()));
  for(unsigned int i=0; i<atoms.size(); ++i) {
    log.printf("  %d", atoms[i].serial());
  }
  log.printf("  \n");
  if(nopbc) log.printf("  without periodic boundary conditions\n");
  else      log.printf("  using periodic boundary conditions\n");

  log.printf("  using graph file:  %s \n", graph_file.c_str());
  dp.init(graph_file);
  if (dp.output_dim() != 3) { throw runtime_error( "invalid graph file! the output dimension should be 3" ); }
  log.printf("  DeepTensor model initialized successfully");

  requestAtoms(atoms); 
}

void DeepDipole::calculate()
{
    if(!nopbc) makeWhole();
    unsigned N=getNumberOfAtoms();
    Vector tot_dipole;
    for(unsigned i=0; i<N; ++i) {
    tot_dipole += getPosition(i);
    }
    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");
    for(unsigned i=0; i<N; i++) {
      setAtomsDerivatives(valuex,i,Vector(1.0,0.0,0.0));
      setAtomsDerivatives(valuey,i,Vector(0.0,1.0,0.0));
      setAtomsDerivatives(valuez,i,Vector(0.0,0.0,1.0));
    }
    setBoxDerivativesNoPbc(valuex);
    setBoxDerivativesNoPbc(valuey);
    setBoxDerivativesNoPbc(valuez);
    valuex->set(tot_dipole[0]);
    valuey->set(tot_dipole[1]);
    valuez->set(tot_dipole[2]);
  }

// // DeepWannier, to be completed
// void DeepDipole::calculate()
// {
//     if(!nopbc) makeWhole();
//     unsigned N=getNumberOfAtoms();  // this function returns the number of needed atoms instead of all atoms
//     // get input tensor,  declare atomic dipole and tot_dipole as tensors
//     // get atomic_dipole = DW_model(input_tensor)
//     // atomic dipole is (ncenter,3) matrix, in our case, ncenter = N/5
//     for(unsigned i=0; i<ncenter; ++i) {
//     tot_dipole += atomic_dipole[i];
//     }
//     Value* valuex=getPntrToComponent("x");
//     Value* valuey=getPntrToComponent("y");
//     Value* valuez=getPntrToComponent("z");
//     // declare DP_grad_x/y/z, evaluate gradient of tot_dipole w.r.t atomic coordinates
//     for(unsigned i=0; i<N; i++) {
//       setAtomsDerivatives(valuex,i,DP_grad_x[i]);
//       setAtomsDerivatives(valuey,i,DP_grad_y[i]);
//       setAtomsDerivatives(valuez,i,DP_grad_z[i]);
//     }
//     // evaluate gradient of DP_grad_x/y/z w.r.t cell, const Tensor&d is (3,3) matrix.
//     setBoxDerivatives(valuex, const Tensor&d)
//     setBoxDerivatives(valuey, const Tensor&d)
//     setBoxDerivatives(valuez, const Tensor&d)
//     valuex->set(tot_dipole[0]);
//     valuey->set(tot_dipole[1]);
//     valuez->set(tot_dipole[2]);
//   }
}
}
