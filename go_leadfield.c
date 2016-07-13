#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void 
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  mxArray *lf;
  double nr, na, dar, F, tmp2;
  double a[3], tmp1[3], delF[3], B[3], cqr[3], um[3], rm[3];
  double *lf_p, *rm_p, *um_p, *R, *Q;
  char str[256];

  if (nrhs<4)
    mexErrMsgTxt("Not enough input arguments");
  if (nrhs>4)
    mexErrMsgTxt("Too many input arguments");

  if (mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 1");
  if (mxGetN(prhs[1])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 2");
  if (mxGetN(prhs[2])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 3");
  if (mxGetN(prhs[3])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 4");
  
  lf = mxCreateDoubleMatrix(1, 1, mxREAL);
  lf_p = mxGetData(lf);

  R      = mxGetData(prhs[0]);
  Q      = mxGetData(prhs[1]);
  rm_p   = mxGetData(prhs[2]);
  um_p   = mxGetData(prhs[3]);

  /* unpack the pointers (in the faintest hope this works) */

  rm[0] = rm_p[0];
  rm[1] = rm_p[1];
  rm[2] = rm_p[2];

  um[0] = um_p[0];
  um[1] = um_p[1];
  um[2] = um_p[2];

  /* compute the difference between this channel and the dipole position */
  a[0] = rm[0] - R[0];
  a[1] = rm[1] - R[1];
  a[2] = rm[2] - R[2];

  nr = sqrt(rm[0]*rm[0] + rm[1]*rm[1] + rm[2]*rm[2]); /* norm(rm) */
  na = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]); /* norm(a) */
  
  tmp1[0] = a[0] / na;
  tmp1[1] = a[1] / na;          /* a/norm(a) */
  tmp1[2] = a[2] / na;

  dar = tmp1[0]*R[0] + tmp1[1]*R[1] + tmp1[2]*R[2]; /* dot(tmp1,R) */

  cqr[0] = Q[1]*R[2] - Q[2]*R[1];
  cqr[1] = Q[2]*R[0] - Q[0]*R[2];   /* cross(Q,R) */
  cqr[2] = Q[0]*R[1] - Q[1]*R[0];

  F = na*(nr*na + nr*nr - (R[0]*rm[0] + R[1]*rm[1] + R[2]*rm[2]));

  delF[0] = ((na*na)/nr + dar + 2*na + 2*nr)*rm[0] - (na + 2*nr + dar)*R[0];
  delF[1] = ((na*na)/nr + dar + 2*na + 2*nr)*rm[1] - (na + 2*nr + dar)*R[1];
  delF[2] = ((na*na)/nr + dar + 2*na + 2*nr)*rm[2] - (na + 2*nr + dar)*R[2];

  tmp2 = cqr[0]*rm[0] + cqr[1]*rm[1] + cqr[2]*rm[2];  /* dot( cqr , rm) */

  B[0] = (F*cqr[0]-tmp2*delF[0])*1e-7/(F*F);
  B[1] = (F*cqr[1]-tmp2*delF[1])*1e-7/(F*F);
  B[2] = (F*cqr[2]-tmp2*delF[2])*1e-7/(F*F);

  lf_p[0] = B[0]*um[0] + B[1]*um[1] + B[2]*um[2]; 

  plhs[0] = lf;

  return;
}