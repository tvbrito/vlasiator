/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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
 */

#ifndef MAGNETOSPHERE_H
#define MAGNETOSPHERE_H

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {
   class Magnetosphere: public TriAxisSearch {
    public:
      Magnetosphere();
      virtual ~Magnetosphere();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell) const;
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,
                                         const int& popID
                                        ) const;
      
    protected:
      Real getDistribValue(
                           creal& x,creal& y, creal& z,
                           creal& vx, creal& vy, creal& vz,
                           creal& dvx, creal& dvy, creal& dvz
                          ) const;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual std::vector<std::array<Real, 3> > getV0(
                                                      creal x,
                                                      creal y,
                                                      creal z
                                                     ) const;
      
      Real T;
      Real tailRho;
      Real V0[3];
      Real ionosphereV0[3];
      Real constBgB[3];
      bool noDipoleInSW;
      Real ionosphereRho;
      Real ionosphereRadius;
      Real ionosphereTaperRadius;
      uint ionosphereGeometry;
      Real center[3];
      Real dipoleScalingFactor;
      Real dipoleMirrorLocationX;
      uint dipoleType;
      uint nSpaceSamples;
      uint nVelocitySamples;
   }; // class Magnetosphere
} // namespace projects

#endif

