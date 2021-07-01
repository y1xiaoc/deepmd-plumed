#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>


namespace PLMD {
namespace dp_plmd {


#ifdef HIGH_PREC
typedef double FLOAT_PREC;
#else
typedef float FLOAT_PREC;
#endif


// model unit in plumed default unit
// 1 Angstrom (deep model) = 0.1 nm (plumed)
const double global_length_unit = 0.1;
// 1 eV (deep model) = 96.487 kJ/mol (plumed)
const double global_energy_unit = 96.487;
// do not convert dipole or polarizability by default 
const double global_dipole_unit = 1.0;
const double global_polar_unit  = 1.0;


// load vector of types from specified file
void load_atype(std::vector<int> & atype, const std::string& filename);


// ravel and unravel indices
struct IndexConverter
{
  const unsigned nout;
  const unsigned natom;
  const unsigned ndim;
public:
  IndexConverter() = delete;
  IndexConverter(const unsigned & natom, const unsigned & ndim = 3): 
    nout(1), natom(natom), ndim(ndim) 
    { if (natom <= 0 || ndim <= 0) throw std::runtime_error("Invalid dimension"); }
  IndexConverter(const unsigned & nout, const unsigned & natom, const unsigned & ndim): 
    nout(nout), natom(natom), ndim(ndim) 
    { if (nout <= 0 || natom <= 0 || ndim <= 0) throw std::runtime_error("Invalid dimension"); }
  // flat (unravel) 2d index 
  unsigned f(const int & iatom, const int & idim);
  // flat (unravel) 3d index
  unsigned f(const int & iout, const int & iatom, const int & idim);
  // ravel to 2d index
  void r(unsigned & iatom, unsigned & idim, const unsigned & multi_idx);
  // ravel to 3d index
  void r(unsigned & iout, unsigned & iatom, unsigned & idim, const unsigned & multi_idx);
};

}
}

