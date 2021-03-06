/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef QUADR_HPP
#define QUADR_HPP

#include <iostream>
#include "functions.hpp"
using namespace std;

/*
  1D,2D,3D Romberg non-singular integration a'la Numerical Recipes.
  Integration bounds must be constants.
  Header file is quadr.H.
  Test program is tstquadr.C.
  Same in Mathematica is tstquadr.ma.
*/

/*
  The iterations are stopped when the results changes by less than absacc.
*/



double Romberg(const T1DFunction& func, double a, double b, double absacc);
double Romberg(const T2DFunction& func, double a, double b, double c, double d, double absacc);
double Romberg(const T3DFunction& func, double a, double b, double c, double d, double e, double f, double absacc);

#endif
