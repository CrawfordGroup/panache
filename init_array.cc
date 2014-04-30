/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*
** \file
** \brief Initialize an array of doubles
** \ingroup CIOMR
*/

#include <cstdlib>
#include <cstring>
#include "Output.h"

namespace panache {

/*!
** init_array(): This function initializes an array of doubles of
** length 'size' and returns a pointer to the first element
**
** \param size = length of array (unsigned long to allow large arrays)
**
** Returns: pointer to new array
**
** \ingroup CIOMR
*/
double * init_array(unsigned long int size)
{
  double *array;

  if ((array = (double *) malloc(size*(unsigned long int)sizeof(double)))
    == NULL) {
    output::printf("init_array: trouble allocating memory \n");
    output::printf("size = %ld\n",size);
    exit(1);
  }
  bzero(array,size*(unsigned long int)sizeof(double));
  return(array);
}

}

