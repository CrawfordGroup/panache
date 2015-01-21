#ifndef PANACHE_COMMON_CPP_H
#define PANACHE_COMMON_CPP_H

#include "panache/BasisSet.h"
#include "panache/Molecule.h"
#include "panache/SimpleMatrix.h"


// Note - exceptions are turned on for the ifstream object
// so that any parsing errors just throw an exeption. Catch those,
// and throw an exception
panache::SharedBasisSet ReadBasisFile(panache::SharedMolecule mol, const std::string & filename);

panache::SharedMolecule ReadMoleculeFile(const std::string & filename);

std::shared_ptr<panache::SimpleMatrix> ReadCMatrixFile(const std::string & filename);

std::vector<double> ReadOrbEnFile(const std::string & filename);

int ReadNocc(const std::string & filename);

std::string GetNextArg(int & i, int argc, char ** argv);

int GetIArg(int & i, int argc, char ** argv);

#endif
