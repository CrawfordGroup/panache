/*! \file
 * \brief Generic, local 4-index tensor storage (header)
 * \ingroup fourindexgroup
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_LOCAL4INDEXTENSOR_H
#define PANACHE_LOCAL4INDEXTENSOR_H

#include "panache/tensorbackend/FourIndexTensor.h"

#include <vector>

namespace panache
{


class LocalQTensor;  // the component 3-index tensor


class Local4IndexTensor : public FourIndexTensor
{
public:
    Local4IndexTensor(LocalQTensor * left, LocalQTensor * right);

protected:
    virtual int GetNLocalIntegrals_(void) const;
    virtual FourIndexIntegral LocalIntegral_(int index);

private:
    LocalQTensor * left_;
    LocalQTensor * right_;

    std::vector<double> contscratch_;

};

} // close namespace panache

#endif
