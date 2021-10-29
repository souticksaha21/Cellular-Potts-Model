#include "stdio.h"
#include "math.h"
#include "mersenne_twister.hpp"

#include "cpm6g_param.cpp"
// includes alpha, lambda_A, epsilon, r, eta; scan over expt g values

int main(void)
{

  // instantiate random number generator
  MTRand rg(1);

  // for file writing
  FILE *fpparam, *fpX, *fpY, *fpNp, *fpt, *fpv, *fpCI, *fpCR;
  char filename[FILENAME_MAX];

  // allocate variables
  int i, j, k, l, z, ig, k1, k2;
  int xmin, xmax, ymin, ymax, x1, x2, y1, y2;
  int r1, r2, r3;
  double r4, r5;
  double t; // s, time
  int sigma1, sigma2; // pixel identities (0 = ECM, 1 = cell)
  int NT; // number of trials in a MC time step
  double A0 = 4*R*R; // cell area (um)
  int N = (int) A0/a/a; // average number of cell pixels
  int L = (int) 2*R/a; // initial cell length in pixels
  int kmax = 10*N; // max pixel number (for file writing)
  int X[kmax]; // cell pixel x locations
  int Y[kmax]; // cell pixel y locations
  int Np; // current number of cell pixels
  int dN; // perimeter change
  double du, w; // energy change and work
  int kn; // index of neighbor
  int nn; // number of neighbor's neighbors
  double tn; // next measurement time
  int b = 1; // buffer for sampling area
  int Nm = (int) T/deltat; // number of measurements (excluding t = 0)
  double xcm; // center of mass x position
  double ycm; // center of mass y position
  double xcmt; // center of mass x position, updated every tau
  double ycmt; // center of mass y position, updated every tau
  double xcms[Nm+1]; // vector of center of mass x positions
  double ycms[Nm+1]; // vector of center of mass y positions
  int icm; // index for cm vectors
  double dxcm, dycm, dxcmt, dycmt; // changes in CM components
  double cmnorm;
  double v; // speed
  double CI; // chemotactic index
  double CR; // chemotactic ratio
  double vsum, dist;
  int ns[kmax], n1[kmax], n2[kmax]; // molecule number samples
  double nbar; // mean of ns
  double lambda, p; // for Poisson sample generation
  double lambda1, p1; // for Poisson sample generation
  double lambda2, p2; // for Poisson sample generation
  double qx, qy; // components of sensing vector
  double px, py; // components of polarization vector
  double g; // scanned parameter
  double gs[4]; // 1/um^4, gradient
  double c0;
  gs[0] = 0; // values
  gs[1] = 1*.6/1e3;
  gs[2] = 5*.6/1e3;
  gs[3] = 50*.6/1e3;
  

  // write parameters
  sprintf(filename,"%s.param.dat",base);
  fpparam = fopen(filename,"w");
  fprintf(fpparam,"%f\n",R);
  fprintf(fpparam,"%f\n",a);
  fprintf(fpparam,"%f\n",tau);
  fprintf(fpparam,"%f\n",T);
  fprintf(fpparam,"%f\n",deltat);
  fclose(fpparam);

  // loop over parameters
  for (ig = 0; ig < 4; ig++){
    g = gs[ig];
    c0 = 0.5*g*1000;
   // chem(ig) = 0.5*g
    printf("g = %f\n",g);
    
    // open stat files
    sprintf(filename,"%s_g%d.v.dat",base,ig);
    fpv = fopen(filename,"w");
    sprintf(filename,"%s_g%d.CI.dat",base,ig);
    fpCI = fopen(filename,"w");
    sprintf(filename,"%s_g%d.CR.dat",base,ig);
    fpCR = fopen(filename,"w");
    
    // loop over trials
    for (z = 0; z < Z; z++){
      //printf("z = %d\n",z);
      
      // initialize polarization vector
      px = 0;
      py = 0;
      
      // initialize cell pixel locations (with buffer zeros at end)
      Np = N;
      for (i = 0; i < L; i++){
	for (j = 0; j < L; j++){
	  k = j + L*i;
	  X[k] = i;
	  Y[k] = j;
	}
      }
      for (k = Np; k < kmax; k++){
	X[k] = 0;
	Y[k] = 0;
      }
      
      // calculate center of mass
      xcm = 0;
      ycm = 0;
      for (k = 0; k < Np; k++){
	xcm = xcm + X[k];
	ycm = ycm + Y[k];
      }
      xcm = xcm/Np;
      ycm = ycm/Np;
      xcmt = xcm;
      ycmt = ycm;
      
      // initialize center of mass change
      dxcmt = 0;
      dycmt = 0;
      
      // simulate
      t = 0;
      tn = 0;
      icm = 0;
      while (t < T){
	
	// every deltat
	if (t >= tn){
	  
	  // record (x,y) for center of mass
	  xcms[icm] = xcm;
	  ycms[icm] = ycm;
	  tn = tn + deltat;
	  icm = icm + 1;
	}
	
	// get Poisson sample at each pixel
	// note that if lambda < 0, then n = 0;
	nbar = 0;
          for (k = 0; k < Np; k++){
            lambda1 = (c0 + g*X[k])*a*a*a;
            n1[k] = 0;
            i=n1[k];
            r5 = rg.randDblExc();
            p1 = -log(r5);
            while (p1 < lambda1){
              r5 = rg.randDblExc();
              p1 = p1 - log(r5);
              n1[k] = n1[k] + 1;
              i=n1[k];
            }
            lambda2 = (float)(Nrecep)/(float)(Np);
            n2[k] = 0;
            j=n2[k];
            r5 = rg.randDblExc();
            p2 = -log(r5);
            while (p2 < lambda2 && j < 1.2*i){
              r5 = rg.randDblExc();
              p2 = p2 - log(r5);
              n2[k] = n2[k] + 1;
              j=n2[k];
            }
            //jm=std::min(i,j);
            ns[k]=std::min(i,j);
            nbar = nbar + ns[k];
           // fprintf(mintest,"%d %d %d \n",i,j,jm);
          }
          nbar = nbar/Np;	
	
	// calculate sensing vector
	qx = 0;
	qy = 0;
	for (k = 0; k < Np; k++){
	  qx = qx + (ns[k]-nbar)*(X[k]-xcm)
	    /sqrt((X[k]-xcm)*(X[k]-xcm) + (Y[k]-ycm)*(Y[k]-ycm));
	  qy = qy + (ns[k]-nbar)*(Y[k]-ycm)
	    /sqrt((X[k]-xcm)*(X[k]-xcm) + (Y[k]-ycm)*(Y[k]-ycm));
	}
	qx = qx/Np;
	qy = qy/Np;
	
	// calculate polarization vector
	if ((dxcmt == 0) && (dycmt == 0)){
	  cmnorm = 1;
	}
	else{
	  cmnorm = sqrt(dxcmt*dxcmt+dycmt*dycmt);
	}
	px = px + r*tau*(-px + epsilon*qx + eta*dxcmt/cmnorm);
	py = py + r*tau*(-py + epsilon*qy + eta*dycmt/cmnorm);
	
	// find min/max x and y locations
	xmin = X[0];
	xmax = X[0];
	ymin = Y[0];
	ymax = Y[0];
	for (k = 0; k < Np; k++){
	  if (X[k] < xmin){ xmin = X[k]; }
	  if (X[k] > xmax){ xmax = X[k]; }
	  if (Y[k] < ymin){ ymin = Y[k]; }
	  if (Y[k] > ymax){ ymax = Y[k]; }
	}
	
	// calculate number of samples
	NT = (xmax-xmin+2*b+1)*(ymax-ymin+2*b+1);
	
	// loop over samples
	for (l = 0; l < NT; l++){
	  
	  // randomly pick a pixel
	  r1 = rg.randInt(xmax-xmin+2*b);
	  r2 = rg.randInt(ymax-ymin+2*b);
	  x1 = (int) (xmin-b + r1);
	  y1 = (int) (ymin-b + r2);
	  
	  // randomly pick a neighbor of that pixel
	  r3 = rg.randInt(3);
	  if (r3 == 0){ x2 = x1-1; y2 = y1; }
	  if (r3 == 1){ x2 = x1+1; y2 = y1; }
	  if (r3 == 2){ x2 = x1; y2 = y1-1; }
	  if (r3 == 3){ x2 = x1; y2 = y1+1; }
	  
	  // find identities of pixel and neighbor
	  sigma1 = 0;
	  sigma2 = 0;
	  for (k = 0; k < Np; k++){
	    if ((X[k] == x1) && (Y[k] == y1)){
	      sigma1 = 1;
	    }
	    if ((X[k] == x2) && (Y[k] == y2)){
	      sigma2 = 1;
	      // save neighbor's index in case it becomes ECM
	      kn = k;
	    }
	  }
	  
	  // if indentities different, try to copy pixel identity to neighbor
	  // case 1: pixel is cell and neighbor is ECM
	  if ((sigma1 == 1) && (sigma2 == 0)){
	    
	    // find perimeter change based on neighbor's neighbors
	    nn = 0;
	    for (k = 0; k < Np; k++){
	      if ((X[k] == x2+1) && (Y[k] == y2)){
		nn = nn + 1;
	      }
	      if ((X[k] == x2-1) && (Y[k] == y2)){
		nn = nn + 1;
	      }
	      if ((X[k] == x2) && (Y[k] == y2+1)){
		nn = nn + 1;
	      }
	      if ((X[k] == x2) && (Y[k] == y2-1)){
		nn = nn + 1;
	      }
	    }
	    if (nn == 1){ dN = 2; }
	    if (nn == 2){ dN = 0; }
	    if (nn == 3){ dN = -2; }
	    if (nn == 4){ dN = -4; }
	    if (nn > 4){ printf("nn = %d\n",nn); }
	    if (nn == 0){ printf("nn = %d\n",nn); }
	    
	    // calculate energy change
	    du = alpha*a*dN
	      + lambdaA*a*a*a*a*((Np+1-N)*(Np+1-N) - (Np-N)*(Np-N));
	    
	    // calculate center of mass change
	    dxcm = (x2-xcm)/(Np+1);
	    dycm = (y2-ycm)/(Np+1);
	    
	    // calculate work
	    w = dxcm*px + dycm*py;
	    
	    // attempt the move
	    r4 = rg.rand();
	    if (r4 < exp(-du+w)){
	      X[Np] = x2;
	      Y[Np] = y2;
	      xcm = xcm + dxcm;
	      ycm = ycm + dycm;
	      Np = Np + 1;
	    }
	  }
	  
	  // case 2: pixel is ECM and neighbor is cell
	  if ((sigma1 == 0) && (sigma2 == 1)){
	    
	    // find perimeter change based on neighbor's neighbors
	    nn = 0;
	    for (k = 0; k < Np; k++){
	      if ((X[k] == x2+1) && (Y[k] == y2)){
		nn = nn + 1;
	      }
	      if ((X[k] == x2-1) && (Y[k] == y2)){
		nn = nn + 1;
	      }
	      if ((X[k] == x2) && (Y[k] == y2+1)){
		nn = nn + 1;
	      }
	      if ((X[k] == x2) && (Y[k] == y2-1)){
		nn = nn + 1;
	      }
	    }
	    if (nn == 0){ dN = -4; }
	    if (nn == 1){ dN = -2; }
	    if (nn == 2){ dN = 0; }
	    if (nn == 3){ dN = 2; }
	    if (nn > 3){ printf("nn = %d\n",nn); }
	    
	    // calculate energy change
	    du = alpha*a*dN
	      + lambdaA*a*a*a*a*((Np-1-N)*(Np-1-N) - (Np-N)*(Np-N));
	    
	    // calculate center of mass change
	    dxcm = -(x2-xcm)/(Np-1);
	    dycm = -(y2-ycm)/(Np-1);
	    
	    // calculate work
	    w = dxcm*px + dycm*py;
	    
	    // attempt the move
	    r4 = rg.rand();
	    if (r4 < exp(-du+w)){
	      // remove the neighbor pixel from the cell list and shift
	      for (k = kn; k < Np-1; k++){
		X[k] = X[k+1];
		Y[k] = Y[k+1];
	      }
	      X[Np-1] = 0;
	      Y[Np-1] = 0;
	      xcm = xcm + dxcm;
	      ycm = ycm + dycm;
	      Np = Np - 1;
	    }
	  }
	}
	
	// calculate total center of mass change in tau
	dxcmt = xcm - xcmt;
	dycmt = ycm - ycmt;
	xcmt = xcm;
	ycmt = ycm;
	
	// update time
	t = t + tau; 
      }
      
      // compute statistics
      // speed
      vsum = 0;
      for (i = 0; i < Nm; i++){
	vsum = vsum
	  + a*sqrt((xcms[i+1]-xcms[i])*(xcms[i+1]-xcms[i]) 
		   + (ycms[i+1]-ycms[i])*(ycms[i+1]-ycms[i]))/deltat;
      }
      v = vsum/Nm;
      
      // CI
      CI = (xcms[Nm]-xcms[0])/
	sqrt((xcms[Nm]-xcms[0])*(xcms[Nm]-xcms[0])
	     + (ycms[Nm]-ycms[0])*(ycms[Nm]-ycms[0]));
      
      // CR
      dist = 0;
      for (i = 0; i < Nm; i++){
	dist = dist + sqrt((xcms[i+1]-xcms[i])*(xcms[i+1]-xcms[i])
			   + (ycms[i+1]-ycms[i])*(ycms[i+1]-ycms[i]));
      }
      CR = sqrt((xcms[Nm]-xcms[0])*(xcms[Nm]-xcms[0])
		+ (ycms[Nm]-ycms[0])*(ycms[Nm]-ycms[0]))/dist;
      
      // write stats
      fprintf(fpv,"%f\n",v);
      fprintf(fpCI,"%f\n",CI);
      fprintf(fpCR,"%f\n",CR);
      
    }
    
    // close stat files
    fclose(fpv);
    fclose(fpCI);
    fclose(fpCR);
  }

  return 0;
  
}
