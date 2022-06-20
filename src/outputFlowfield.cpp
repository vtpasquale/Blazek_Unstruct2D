/// @file outputFlowfield.cpp
///
/// Output of the flow field in Vis2D format.
///
//*****************************************************************************
//
//  (c) J. Blazek, CFD Consulting & Analysis, www.cfd-ca.de
//  Created February 15, 2014
//  Last modification: January 9, 2017
//
//=============================================================================
//
//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public License
//  as published by the Free Software Foundation; either version 2
//  of the License, or (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//
//*****************************************************************************

#include <ctime>
#include <sstream>
#include <stdexcept>
#include "output.h"

using namespace std;

/// Writes out selected quantities for the whole flow field in Vis2D format.
///
/// @param geometry    geometrical data
/// @param fluidProps  fluid properties
/// @param bndConds    boundary conditions
/// @param iter        current iteration
/// @exception         std::runtime_error  file cannot be opened
///
void Output::Flowfield( const Geometry &geometry, const FluidProps &fluidProps,
                        const BndConds &bndConds, int iter ) const
{
  int nn = geometry.nndInt;
  int  i, m;
  REAL rrho[nn], u[nn], v[nn], e[nn], press[nn], temp[nn], c[nn], mach[nn],
       ttot[nn], ptot[nn], machis[nn], ptloss[nn], pratio[nn], ptotinf,
       gam1[nn], ggm1[nn], visc[nn];

  ptotinf = 0.;
  if (external)
  {
    REAL gam1Temp    = fluidProps.gamma - 1.0;
    REAL ggm1Temp    = fluidProps.gamma/gam1Temp;
    ptotinf = bndConds.pinf*POW((1.0+0.5*gam1Temp*bndConds.machinf*bndConds.machinf),ggm1Temp);
  }

 // Compute Output data
  for (i=0; i<geometry.nndInt; i++)
  {
    rrho[i]  = 1.0/fluidProps.cv[i].dens;
    u[i]     = fluidProps.cv[i].xmom*rrho[i];
    v[i]     = fluidProps.cv[i].ymom*rrho[i];
    e[i]     = fluidProps.cv[i].ener*rrho[i];
    press[i] = fluidProps.dv[i].press;
    temp[i]  = fluidProps.dv[i].temp;
    c[i]     = fluidProps.dv[i].csoun;
    gam1[i]  = fluidProps.dv[i].gamma - 1.0;
    ggm1[i]  = fluidProps.dv[i].gamma/gam1[i];
    mach[i]  = SQRT(u[i]*u[i]+v[i]*v[i])/c[i];
    ttot[i]  = (e[i]+press[i]*rrho[i])/fluidProps.dv[i].cpgas;
    ptot[i]  = press[i]*POW(ttot[i]/temp[i],ggm1[i]);
    if (external)
    {
      ptloss[i] = 1.0 - ptot[i]/ptotinf;
      pratio[i] = ptotinf/press[i];
    }
    else
    {
      ptloss[i] = 1.0 - ptot[i]/bndConds.ptinl;
      pratio[i] = bndConds.ptinl/press[i];
    }
    machis[i] = (POW(pratio[i],1.0/ggm1[i])-1.0)*2.0/gam1[i];
    machis[i] = MAX(machis[i], 0.0);
    machis[i] = SQRT(machis[i]);

    if (fluidProps.equsType == Equations::NavierStokes)
      visc[i] = fluidProps.dvLam[i].mue;
    else
      visc[i] = 0.0;
  }

  // open file
  stringstream str;
  str << right;
  str.width(5);
  str.fill('0');
  str << iter;
  string fname = fnameFlow + str.str() + ".vtk";
  
  ofstream stream( fname,ofstream::out | ofstream::trunc );
  if (stream.fail()) throw runtime_error( "could not open plot file (flow field) for writing." );

  // Header
  stream << "# vtk DataFile Version 3.0" << endl;
  stream << "U2D Output" << endl;
  stream << "ASCII" << endl;

  // POINTS
  stream << endl << "DATASET UNSTRUCTURED_GRID" << endl;
  stream << "POINTS " << geometry.nndInt  << " float" << endl;
  for (int i = 0; i <  geometry.nndInt ; i++)
      stream << geometry.coords[i].x << " " << geometry.coords[i].y << " 0.0" << endl;

  // CELLS
  stream << "CELLS " << geometry.nTria << " " << 4*geometry.nTria << endl;
  for (int i = 0; i < geometry.nTria; i++)
      stream << "3 " << geometry.tria[i].node[0] << " " << geometry.tria[i].node[1] << " " << geometry.tria[i].node[2] << endl;

  // CELL TYPES
  stream << "CELL_TYPES " << geometry.nTria << endl;
  for (int i = 0; i < geometry.nTria; i++)
      stream << "5" << endl;

  // Write output data 
  stream << "POINT_DATA " << geometry.nndInt << endl;
  for (m=0; m<MXQFIELD; m++)
  {
    if (varOn[m])
    {
    switch (m)
      {
        case 0:   // density
          stream << "SCALARS " << "rho" << " float" << endl;
          stream << "LOOKUP_TABLE default" << endl;
          for (i=0; i<geometry.nndInt; i++)
            stream << fluidProps.cv[i].dens << endl;
          break;
        case 1:   // u-velocity
          stream << "SCALARS " << "u" << " float" << endl;
          stream << "LOOKUP_TABLE default" << endl;
          for (i=0; i<geometry.nndInt; i++)
            stream << u[i] << endl;
          break;
        case 2:   // v-velocity
            stream << "SCALARS " << "v" << " float" << endl;
            stream << "LOOKUP_TABLE default" << endl;
            for (i=0; i<geometry.nndInt; i++)
                stream << v[i] << endl;
          break;
        case 3:   // static pressure
          stream << "SCALARS " << "p" << " float" << endl;
          stream << "LOOKUP_TABLE default" << endl;
          for (i=0; i<geometry.nndInt; i++)
            stream << press[i] << endl;
          break;
        case 4:   // total pressure
          stream << "SCALARS " << "ptot" << " float" << endl;
          stream << "LOOKUP_TABLE default" << endl;
          for (i=0; i<geometry.nndInt; i++)
            stream << ptot[i] << endl;
          break;
        case 5:   // static temperature
          stream << "SCALARS " << "temp" << " float" << endl;
          stream << "LOOKUP_TABLE default" << endl;
          for (i=0; i<geometry.nndInt; i++)
            stream << temp[i] << endl;
          break;
        case 6:   // total temperature
          stream << "SCALARS " << "ttot" << " float" << endl;
          stream << "LOOKUP_TABLE default" << endl;
          for (i=0; i<geometry.nndInt; i++)
            stream << ttot[i] << endl;
          break;
        case 7:   // local Mach number
          stream << "SCALARS " << "mach" << " float" << endl;
          stream << "LOOKUP_TABLE default" << endl;
          for (i=0; i<geometry.nndInt; i++)
            stream << mach[i] << endl;
          break;
        case 8:   // isentropic Mach number
          stream << "SCALARS " << "machis" << " float" << endl;
          stream << "LOOKUP_TABLE default" << endl;
          for (i=0; i<geometry.nndInt; i++)
            stream << machis[i] << endl;
          break;
        case 9:   // total pressure loss
          stream << "SCALARS " << "ptloss" << " float" << endl;
          stream << "LOOKUP_TABLE default" << endl;
          for (i=0; i<geometry.nndInt; i++)
            stream << ptloss[i] << endl;
          break;
        case 10:  // laminar viscosity coefficient
          stream << "SCALARS " << "visc" << " float" << endl;
          stream << "LOOKUP_TABLE default" << endl;
          for (i=0; i<geometry.nndInt; i++)
            stream << visc[i] << endl;
          break;
      }
    }
  }

  // stream << scientific;
  // stream.precision(8);

  stream.close();
}
