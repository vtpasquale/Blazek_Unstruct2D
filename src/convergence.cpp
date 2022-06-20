/// @file convergence.cpp
///
/// Output of the convergence history.
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

#include "solver.h"

using namespace std;

/// Monitors the convergence, prints it out and stores it in a file. For external
/// flow, it also computes and prints out the lift, the drag and the moment coefficients.
/// For internal flow, it computes and prints out the mass flow and the mass flow ratio.
///
void Solver::Convergence( Output &output )
{
  int  idr;
  REAL drmax;

  // compute the 2-norm of the density change

  timeDiscr.DensityChange( geometry,fluidProps,drho,drmax,idr );

  if (iter == 1)
  {
    drho1 = drho;
    drho1 = MAX(drho1,1.0e-33);
    drho  = 1.0;
  }
  else
  {
    drho = drho/drho1;
    drho = MAX(drho,1.0e-33);
  }

  // compute forces & moments (external flow)

  if (bndConds.flowType == FlowType::External)
  {
    Forces();
    bndConds.cl = cl;  // update lift coefficient (used for vortex correction at far-field)

  // compute mass flow and mass flow ratio (internal flow)
  }
  else
  {
    MassFlow();
  }

  //  print out / store

  output.WriteConvergence( iter,drho,drmax,idr,cl,cd,cm,mflow,mfratio );

  cout << scientific << right;
  cout.precision(4);
  cout.fill(' ');

  if (bndConds.flowType == FlowType::External)
  {
    cout.width(5);
    cout << iter;
    cout.width(14);
    cout << LOG10(drho);
    cout.width(14);
    cout << drmax;
    cout.width(8);
    cout << idr;
    cout.width(14);
    cout << cl;
    cout.width(14);
    cout << cd;
    cout.width(14);
    cout << cm << endl;
  }
  else
  {
    cout.width(5);
    cout << iter;
    cout.width(14);
    cout << LOG10(drho);
    cout.width(14);
    cout << drmax;
    cout.width(8);
    cout << idr;
    cout.width(14);
    cout << mflow;
    cout.width(14);
    cout << mfratio << endl;
  }
}
