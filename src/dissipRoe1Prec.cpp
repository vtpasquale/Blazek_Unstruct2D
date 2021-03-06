/// @file dissipRoe1Prec.cpp
///
/// Computation of 1st-order upwind dissipation with preconditioning.
///
//*****************************************************************************
//
//  (c) J. Blazek, CFD Consulting & Analysis, www.cfd-ca.de
//  Created February 15, 2014
//  Last modification: September 20, 2014
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

#include "spaceDiscr.h"

/// Computes upwind dissipation according to 1st-order Roe's flux-difference
/// splitting scheme. The dissipation terms are multiplied by a preconditioning
/// matrix for low Mach numbers (see Eqs. (9.54), (9.82)).
///
/// @param geometry    geometrical data
/// @param fluidProps  fluid properties
/// @param precond     low Mach-number preconditioning
/// @param beta        dissipation-blending coefficient
///
void SpaceDiscr::DissipRoe1Prec( const Geometry &geometry, const FluidProps &fluidProps,
                                 const Precond &precond, REAL beta )
{
  int  i, j, ie;
  REAL beta5, ds, gam1, rrho, rl, ul, vl, pl, tl, hl, rr, ur, vr, pr, tr, hr,
       rav, dd, dd1, uav, vav, pav, tav, hav, cav, q2a, uv, h1, h2, h5, delta,
       eabs1, eabs2, eabs5, a1, a4, a5, cc2, ra1g, rhop, rhoT, hT, theta;
  REAL fd[5], nVec[3], gtpMat[5][5], tp1Mat[5][5];
  REAL wVec[5], wpVec[5], wrlVec[5], dumVec[5];

  CONSVARS       *cv   = fluidProps.cv;
  DEPVARS        *dv   = fluidProps.dv;
  NODE           *sij  = geometry.sij;
  Geometry::EDGE *edge = geometry.edge;

  beta5 = 0.5*beta;

  nVec[2]   = 0.0;  // no 3rd dimension here
  wrlVec[3] = 0.0;
  wVec[3]   = 0.0;
  wpVec[3]  = 0.0;

  for (ie=0; ie<geometry.nEdges; ie++)
  {
    i       = edge[ie].i;
    j       = edge[ie].j;
    ds      = SQRT(sij[ie].x*sij[ie].x+sij[ie].y*sij[ie].y);
    nVec[0] = sij[ie].x/ds;
    nVec[1] = sij[ie].y/ds;

    // left & right state

    rrho = 1.0/cv[i].dens;
    rl   = cv[i].dens;
    ul   = cv[i].xmom*rrho;
    vl   = cv[i].ymom*rrho;
    pl   = dv[i].press;
    tl   = dv[i].temp;
    hl   = (pl+cv[i].ener)*rrho;

    rrho = 1.0/cv[j].dens;
    rr   = cv[j].dens;
    ur   = cv[j].xmom*rrho;
    vr   = cv[j].ymom*rrho;
    pr   = dv[j].press;
    tr   = dv[j].temp;
    hr   = (pr+cv[j].ener)*rrho;

    // Roe's average

    rav  = SQRT(rl*rr);
    gam1 = 0.5*(dv[i].gamma+dv[j].gamma) - 1.0;
    dd   = rav/rl;
    dd1  = 1.0/(1.0+dd);
    uav  = (ul+dd*ur)*dd1;
    vav  = (vl+dd*vr)*dd1;
    pav  = (pl+dd*pr)*dd1;
    tav  = (tl+dd*tr)*dd1;
    hav  = (hl+dd*hr)*dd1;
    q2a  = uav*uav + vav*vav;
    cav  = SQRT(gam1*(hav-0.5*q2a));
    uv   = uav*nVec[0] + vav*nVec[1];

    // preconditioning

    wrlVec[0] = rr - rl;
    wrlVec[1] = cv[j].xmom - cv[i].xmom;
    wrlVec[2] = cv[j].ymom - cv[i].ymom;
    wrlVec[4] = cv[j].ener - cv[i].ener;

    wVec[0] = rav;
    wVec[1] = rav*uav;
    wVec[2] = rav*vav;
    wVec[4] = rav*hav - pav;

    wpVec[0] = pav;
    wpVec[1] = uav;
    wpVec[2] = vav;
    wpVec[4] = tav;

    rhop  =  rav/pav;
    rhoT  = -rav/tav;
    hT    = 0.5*(dv[i].cpgas+dv[j].cpgas);
    theta = precond.ComputeTheta( gam1+1.0,cav,q2a );

    delta = epsEntr*cav;
    h2    = ABS(uv);
    eabs2 = EntropyCorr( h2,delta );
    a1    = wVec[0]*rhop*hT + rhoT;
    ra1g  = 1.0/(wVec[0]*theta*hT + rhoT);
    a4    = a1*ra1g;
    a5    = wVec[0]*hT*ra1g;
    cc2   = 0.5*SQRT((uv*uv)*((a4-1.0)*(a4-1.0))+4.0*a5);
    h2    = 0.5*(a4+1.0)*uv;
    h1    = ABS(h2 + cc2);
    h5    = ABS(h2 - cc2);
    eabs1 = EntropyCorr( h1,delta );
    eabs5 = EntropyCorr( h5,delta );

    precond.RightEigenvec( wVec,wpVec,nVec,uv,hav,theta,rhop,rhoT,0.0,hT,gtpMat );
    precond.MatprodTp1_P1( wVec,wpVec,nVec,uv,hav,theta,rhop,rhoT,0.0,hT,q2a,tp1Mat );

    // upwind fluxes

    precond.MatVecProd5( tp1Mat,wrlVec,dumVec );
    dumVec[0] *= eabs2;
    dumVec[1] *= eabs2;
    dumVec[2] *= eabs2;
    dumVec[3] *= eabs1;
    dumVec[4] *= eabs5;
    precond.MatVecProd5( gtpMat,dumVec,fd );

    // edge contributions to dissipation

    ds           *= beta5;
    diss[i].dens += fd[0]*ds;
    diss[i].xmom += fd[1]*ds;
    diss[i].ymom += fd[2]*ds;
    diss[i].ener += fd[4]*ds;

    diss[j].dens -= fd[0]*ds;
    diss[j].xmom -= fd[1]*ds;
    diss[j].ymom -= fd[2]*ds;
    diss[j].ener -= fd[4]*ds;
  }
}
