#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "/Applications/NV5/idl91/external/include/idl_export.h"

static void bailout(char *msg)
{
  IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP, msg);
}

static void info(char *msg)
{
  IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, msg);
}

void assert_numeric(const char *description, IDL_VPTR arg)
{
  char msg[256];
  if (arg->type != IDL_TYP_DOUBLE && arg->type != IDL_TYP_FLOAT && arg->type != IDL_TYP_INT &&
      arg->type != IDL_TYP_UINT && arg->type != IDL_TYP_LONG && arg->type != IDL_TYP_ULONG) {
    char msg[256];
    snprintf(msg, sizeof(msg), "%s must be numeric", description);
    bailout(msg);
  }
}

void check_params(int argc, IDL_VPTR Argv[])
{
  IDL_ENSURE_ARRAY(Argv[0]);
  IDL_ENSURE_ARRAY(Argv[1]);
  if (Argv[1]->value.arr->n_dim != 1) {
    bailout("A (coefficients) must be a 1-dimensional array");
  }
  assert_numeric("X (1st arg)", Argv[0]);
  assert_numeric("A (coefficients)", Argv[1]);
}

void make_permanent_arr_from_template(IDL_VPTR x_vptr, IDL_VPTR dest)
{
  IDL_VPTR tmp;
  // First temporary then make permanent with IDL_VarCopy
  IDL_VarMakeTempFromTemplate(x_vptr, x_vptr->type, 0 /*not struct*/, &tmp, FALSE);
  IDL_VarCopy(tmp, dest);
}

void make_pder_array(IDL_VPTR x_vptr, IDL_VPTR a_vptr, IDL_VPTR pder, IDL_MEMINT pder_dim[2])
{
  if (pder) {
    pder_dim[0] = x_vptr->value.arr->n_elts;
    pder_dim[1] = a_vptr->value.arr->n_elts;
    IDL_VPTR tmp;
    IDL_MakeTempArray(IDL_TYP_DOUBLE, 2, pder_dim, IDL_ARR_INI_NOP, &tmp);
    IDL_VarCopy(tmp, pder);
  }
}

// ****************************************************************************************************
// ******************************************* COMP_POLY **********************************************
// ****************************************************************************************************

/*; Use         : COMP_POLY,X,A,F [,PDER]*/
static void COMP_POLY(int argc, IDL_VPTR Argv[], char *argk)
{
  check_params(argc, Argv);

  IDL_VPTR x_vptr = IDL_CvtDbl(1, Argv);
  IDL_VPTR a_vptr = IDL_CvtDbl(1, Argv + 1);

  IDL_VPTR f_vptr = Argv[2]; /* Output array */
  IDL_VPTR pder_vptr = argc > 3 ? Argv[3] : NULL;

  make_permanent_arr_from_template(x_vptr, f_vptr);

  double *x = (void *) x_vptr->value.arr->data;
  double *a = (void *) a_vptr->value.arr->data;
  double *f = (void *) f_vptr->value.arr->data;
  double *pder = NULL;

  IDL_MEMINT pder_dim[2];
  if (pder_vptr) {
    make_pder_array(x_vptr, a_vptr, pder_vptr, pder_dim);
    pder = (void *) pder_vptr->value.arr->data;
  }

  for (IDL_MEMINT ix = 0; ix < x_vptr->value.arr->n_elts; ix++) {
    f[ix] = 0.0;
    for (IDL_MEMINT aix = 0; aix < a_vptr->value.arr->n_elts; aix++) {
      f[ix] += a[aix] * pow(x[ix], aix);
      if (pder_vptr) {
        // \frac{\partial}{\partial a[aix]} a[aix] * pow(...) = pow(...)
        pder[ix + aix * pder_dim[0]] = pow(x[ix], aix);
      }
    }
  }
  // These might be temp. due to type conversion
  IDL_DELTMP(x_vptr);
  IDL_DELTMP(a_vptr);
}

// ****************************************************************************************************
// ******************************************* COMP_GAUSS *********************************************
// ****************************************************************************************************

/*; Use         : COMP_GAUSS,X,A,F [,PDER]*/
static void COMP_GAUSS(int argc, IDL_VPTR Argv[], char *argk)
{
  char msg[256];
  info("COMP_GAUSS: Running!");
  check_params(argc, Argv);

  IDL_VPTR x_vptr = IDL_CvtDbl(1, Argv); /* Input array */
  IDL_VPTR a_vptr = IDL_CvtDbl(1, Argv + 1); /* Coefficients */
  IDL_VPTR f_vptr = Argv[2]; /* Output array */
  IDL_VPTR pder_vptr = argc > 3 ? Argv[3] : NULL; /* Partial derivatives, optional */

  make_permanent_arr_from_template(x_vptr, f_vptr);

  double *a = (void *) a_vptr->value.arr->data;
  double *f = (void *) f_vptr->value.arr->data;
  double *x = (void *) x_vptr->value.arr->data;
  double *pder = NULL;

  IDL_MEMINT pder_dim[2];
  if (pder_vptr) {
    make_pder_array(x_vptr, a_vptr, pder_vptr, pder_dim);
    pder = (void *) pder_vptr->value.arr->data;
  }

  for (IDL_MEMINT ix = 0; ix < x_vptr->value.arr->n_elts; ix++) {
    double z = (x[ix] - a[1]) / a[2];
    double z2 = z * z;
    double kern = exp(-z2 * 0.5);
    kern = (z2 < 1000.0) ? kern : 0.0; // Avoid exp overflow
    f[ix] = a[0] * kern;

    if (pder_vptr) {
      pder[ix + 0 * pder_dim[0]] = kern;
      pder[ix + 1 * pder_dim[0]] = f[ix] * z / a[2];
      pder[ix + 2 * pder_dim[0]] = pder[ix + 1 * pder_dim[0]] * z;
    }
  }
  // These might be temp. due to type conversion
  IDL_DELTMP(x_vptr);
  IDL_DELTMP(a_vptr);
}

/*
  z = (x-a(1))/a(2)
  z2 = z*z
  ix = where(z2 LT 1000.0) ;; Exp(-1000) == 0 unless quadruple precision
  kern = exp(-z2(ix)*0.5)

  f(ix) = a(0)*kern

  pder(ix,0) = kern
  pder(ix,1) = f(ix) * z(ix)/a(2)
  pder(ix,2) = pder(ix,1) * z(ix)
  help,pder
*/

/* Testing:
x = [500.000, 500.200, 500.400, 500.600, 500.800, 501.000, 501.200, 501.400, 501.600, 501.800, 502.000, 502.200, 502.400, 502.600, 502.800, 503.000, 503.200, 503.400, 503.600, 503.800] 
a = [45., 502., 0.4]
comp_gauss,x,a,f,pder
window,0
plot,x,f
oplot,x,pder[*,0]*10
oplot,x,pder[*,1]/10.
oplot,x,pder[*,2]/20.
*/

// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************

int IDL_Load(void)
{
  static IDL_SYSFUN_DEF2 pro_def[] = {{(IDL_SYSRTN_GENERIC) COMP_POLY, "COMP_POLY", 3, 4, 0, 0},
                                      {(IDL_SYSRTN_GENERIC) COMP_GAUSS, "COMP_GAUSS", 3, 4, 0, 0}};
  return IDL_SysRtnAdd(pro_def, FALSE, 2);
}
