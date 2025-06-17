#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "/Applications/NV5/idl91/external/include/idl_export.h"

#define MASK 1 /* Enable/disable keyword mask */
#define IKWOF(a) IDL_KW_OFFSETOF(a)

static void bailout(char *msg)
{
  IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP, msg);
}

static void info(char *msg)
{
  IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_INFO, msg);
}

#define NO_STRUCT_DEF 0

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
  IDL_ENSURE_ARRAY(Argv[0]); /* Ensure X is an array */
  IDL_ENSURE_ARRAY(Argv[1]); /* Ensure A is an array */
  if (Argv[1]->value.arr->n_dim != 1) {
    bailout("A (coefficients) must be a 1-dimensional array");
  }
  assert_numeric("X (1st arg)", Argv[0]);
  assert_numeric("A (coefficients)", Argv[1]);
}

void make_arr_from_template(IDL_VPTR x_vptr, IDL_VPTR dest)
{
  IDL_VPTR tmp;
  // Make F output array, first temporary then make permanent with IDL_VarCopy
  IDL_VarMakeTempFromTemplate(x_vptr, x_vptr->type, NO_STRUCT_DEF, &tmp, FALSE);
  IDL_VarCopy(tmp, dest);
}

void make_pder_array(IDL_VPTR x_vptr, IDL_VPTR a_vptr, IDL_VPTR pder, IDL_MEMINT pder_dim[2])
{
  if (pder) {
    pder_dim[0] = x_vptr->value.arr->n_elts; // Number of elements in x_vptr
    pder_dim[1] = a_vptr->value.arr->n_elts; // Number of coefficients in a_vptr
    IDL_VPTR tmp;
    IDL_MakeTempArray(IDL_TYP_FLOAT, 2, pder_dim, IDL_ARR_INI_NOP, &tmp);
    IDL_VarCopy(tmp, pder);
  }
}

/*; Use         : COMP_POLY,X,A,F [,PDER]*/
static void COMP_POLY(int argc, IDL_VPTR Argv[], char *argk)
{
  check_params(argc, Argv);

  IDL_VPTR x_vptr = IDL_CvtFlt(1, Argv); /* Input array */
  IDL_VPTR a_vptr = IDL_CvtFlt(1, Argv + 1); /* Coefficients */

  IDL_VPTR f_vptr = Argv[2]; /* Output array */
  IDL_VPTR pder_vptr = argc > 3 ? Argv[3] : NULL; /* Partial derivatives, optional */

  make_arr_from_template(x_vptr, f_vptr);

  float *x = (void *) x_vptr->value.arr->data;
  float *a = (void *) a_vptr->value.arr->data;
  float *f = (void *) f_vptr->value.arr->data;
  float *pder = NULL;

  IDL_MEMINT pder_dim[2];
  if (pder_vptr) {
    make_pder_array(x_vptr, a_vptr, pder_vptr, pder_dim);
    pder = (void *) pder_vptr->value.arr->data;
  }

  for (IDL_MEMINT xindex = 0; xindex < x_vptr->value.arr->n_elts; xindex++) {
    f[xindex] = 0.0;
    for (IDL_MEMINT aindex = 0; aindex < a_vptr->value.arr->n_elts; aindex++) {
      f[xindex] += a[aindex] * pow(x[xindex], aindex);
      if (pder_vptr) {
        // \frac{\dell }{\dell C} C * x^j  = x^j
        pder[xindex + aindex * pder_dim[0]] = pow(x[xindex], aindex);
      }
    }
  }
  // These might be temp. due to type conversion
  IDL_DELTMP(x_vptr);
  IDL_DELTMP(a_vptr);
}

// **************************************************************************

/*; Use         : COMP_GAUSS,X,A,F [,PDER]*/
static void COMP_GAUSS(int argc, IDL_VPTR Argv[], char *argk)
{
  char msg[256];
  info("COMP_GAUSS: Running!");
  check_params(argc, Argv);

  IDL_VPTR x_vptr = IDL_CvtFlt(1, Argv); /* Input array */
  IDL_VPTR a_vptr = IDL_CvtFlt(1, Argv + 1); /* Coefficients */
  IDL_VPTR f_vptr = Argv[2]; /* Output array */
  IDL_VPTR pder_vptr = argc > 3 ? Argv[3] : NULL; /* Partial derivatives, optional */

  make_arr_from_template(x_vptr, f_vptr);

  float *a = (void *) a_vptr->value.arr->data;
  float *f = (void *) f_vptr->value.arr->data;
  float *x_data = (void *) x_vptr->value.arr->data;
  float *pder = NULL;

  // snprintf(msg, sizeof(msg), "a = [%.2f, %.2f, %.2f]", a[0], a[1], a[2]);
  // info(msg);
  // snprintf(msg, sizeof(msg), "x_vptr->value.arr->n_elts = %lld", x_vptr->value.arr->n_elts);
  // info(msg);
  IDL_MEMINT pder_dim[2];
  if (pder_vptr) {
    make_pder_array(x_vptr, a_vptr, pder_vptr, pder_dim);
    pder = (void *) pder_vptr->value.arr->data;
  }

  for (IDL_MEMINT ix = 0; ix < x_vptr->value.arr->n_elts; ix++) {
    float x = x_data[ix];
    /* z = (x-a(1))/a(2) */
    float z = (x - a[1]) / a[2];
    /* z2 = z*z */
    float z2 = z * z;
    /* kern = exp(-z2(ix)*0.5) */
    float kern = exp(-z2 * 0.5);
    kern = (z2 < 1000.0) ? kern : 0.0; // Avoid exp overflow
    /* f(ix) = a(0)*exp(-z2(ix)*0.5) */
    f[ix] = a[0] * kern;
    // snprintf(msg, sizeof(msg), "COMP_GAUSS: ix = %lld, x = %.2f, z = %.2f, kern = %.6f", ix, x, z, kern);
    // info(msg);

    if (pder_vptr) {
      /* pder(ix,0) = kern  */
      pder[ix + 0 * pder_dim[0]] = kern;
      /* pder[ix,1] = f(ix) * z(ix)/a(2) */
      pder[ix + 1 * pder_dim[0]] = f[ix] * z / a[2];
      /* pder[ix,2] = pder(ix,1) * z(ix) */
      pder[ix + 2 * pder_dim[0]] = pder[ix + 1 * pder_dim[0]] * z;
    }
  }
  // These might be temp. due to type conversion
  IDL_DELTMP(x_vptr);
  IDL_DELTMP(a_vptr);
}

/*
PRO comp_gauss,x,a,f,pder
   message, /info, "COMP_GAUSS running"
  on_error,0

  nx = n_elements(x)

  z = (x-a(1))/a(2)
  z2 = z*z
  ix = where(z2 LT 1000.0) ;; Exp(-1000) == 0 unless quadruple precision


  f = make_array(size=size(z))

  IF ix(0) EQ -1 THEN BEGIN
     IF n_params() EQ 4 THEN BEGIN
        pder = fltarr(nx,3)
     END
     return
  END


  IF n_params() EQ 3 THEN BEGIN
     f(ix) = a(0)*exp(-z2(ix)*0.5)
     return
  END

  kern = exp(-z2(ix)*0.5)

  f(ix) = a(0)*kern

  pder = fltarr(nx,3)
  pder(ix,0) = kern         ;;
  pder(ix,1) = f(ix) * z(ix)/a(2)   ;; a(0)exp(-0.5*(x-a(1))^2/a(2)^2) * (x-a(1))/a(2)^2
  pder(ix,2) = pder(ix,1) * z(ix)   ;; a(0)exp(...) * (x-a(1))^2/a(2)^3
  help,pder

END
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

/*
Testing:
x=(findgen(21)-10)/5 + 500
comp_gauss,x,[1.01,500,1],f,pder
window,0
plot,x,f
plot,x,f,yrange=[-1.2,1.2]
oplot,x,pder[*,0]
plot,x,f,yrange=[-1.2,1.2],ystyle=1
oplot,x,pder[*,0]
oplot,x,pder[*,1]
oplot,x,pder[*,2]
*/
/*

 */

int IDL_Load(void)
{
  static IDL_SYSFUN_DEF2 pro_def[] = {{(IDL_SYSRTN_GENERIC) COMP_POLY, "COMP_POLY", 3, 4, 0, 0},
                                      {(IDL_SYSRTN_GENERIC) COMP_GAUSS, "COMP_GAUSS", 3, 4, 0, 0}};
  return IDL_SysRtnAdd(pro_def, FALSE, 2);
}
