# include "Common.h"


void 
PLMD::
dp_plmd::
load_atype(std::vector<int> & atype, 
           const std::string & filename) 
{
  std::ifstream fin(filename);
  int element;
  atype.clear();
  while (fin >> element) {
    atype.push_back(element);
  }
}


unsigned 
PLMD::
dp_plmd::
IndexConverter::
f(const int & iatom, 
  const int & idim)
{
  return PLMD::dp_plmd::IndexConverter::f(0, iatom, idim);
}

unsigned 
PLMD::
dp_plmd::
IndexConverter::
f(const int & iout,
  const int & iatom, 
  const int & idim)
{
  int io = iout  >= 0 ? iout  : nout  + iout;
  int ia = iatom >= 0 ? iatom : natom + iatom;
  int id = idim  >= 0 ? idim  : ndim  + idim;
  if (io < 0 || (io >= nout && nout != 1) || 
      ia < 0 || ia >= natom ||
      id < 0 || id >= ndim) { 
    throw std::runtime_error ( "Index out of range!" ); 
  }
  return io * (natom * ndim) + ia * (ndim) + id;
}

void 
PLMD::
dp_plmd::
IndexConverter::
r(unsigned & iatom, 
  unsigned & idim,
  const unsigned & multi_idx)
{
  unsigned tmp_iout;
  PLMD::dp_plmd::IndexConverter::r(tmp_iout, iatom, idim, multi_idx);
  if (tmp_iout > 0) { throw std::runtime_error ( "Index out of range!" ); }
}

void 
PLMD::
dp_plmd::
IndexConverter::
r(unsigned & iout,
  unsigned & iatom, 
  unsigned & idim,
  const unsigned & multi_idx)
{
  iout  = (multi_idx / natom / ndim);
  iatom = (multi_idx / ndim) % natom;
  idim  = multi_idx % ndim;
  if (iout >= nout && nout != 1) { throw std::runtime_error ( "Index out of range!" ); }
}

