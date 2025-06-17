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

void make_arr_from_template(IDL_VPTR x, IDL_VPTR dest)
{
  IDL_VPTR tmp;
  // Make F output array, first temporary then make permanent with IDL_VarCopy
  IDL_VarMakeTempFromTemplate(x, x->type, NO_STRUCT_DEF, &tmp, FALSE);
  IDL_VarCopy(tmp, dest);
}

void make_pder_array(IDL_VPTR x, IDL_VPTR a, IDL_VPTR pder, IDL_MEMINT pder_dim[2])
{
  if (pder) {
    pder_dim[0] = x->value.arr->n_elts; // Number of elements in x
    pder_dim[1] = a->value.arr->n_elts; // Number of coefficients in a
    IDL_VPTR tmp;
    IDL_MakeTempArray(IDL_TYP_DOUBLE, 2, pder_dim, IDL_ARR_INI_NOP, &tmp);
    IDL_VarCopy(tmp, pder);
  }
}

/*; Use         : COMP_POLY,X,A,F [,PDER]*/
static void COMP_POLY(int argc, IDL_VPTR Argv[], char *argk)
{
  check_params(argc, Argv);

  IDL_VPTR x = IDL_CvtDbl(1, Argv);     /* Input array */
  IDL_VPTR a = IDL_CvtDbl(1, Argv + 1); /* Coefficients */

  IDL_VPTR f = Argv[2];                      /* Output array */
  IDL_VPTR pder = argc > 3 ? Argv[3] : NULL; /* Partial derivatives, optional */

  make_arr_from_template(x, f);

  double *x_data = (void *) x->value.arr->data;
  double *a_data = (void *) a->value.arr->data;
  double *f_data = (void *) f->value.arr->data;
  double *pder_data = NULL;

  IDL_MEMINT pder_dim[2];
  if (pder) {
    make_pder_array(x, a, pder, pder_dim);
    pder_data = (void *) pder->value.arr->data;
  }

  for (IDL_MEMINT xindex = 0; xindex < x->value.arr->n_elts; xindex++) {
    f_data[xindex] = 0.0;
    for (IDL_MEMINT aindex = 0; aindex < a->value.arr->n_elts; aindex++) {
      f_data[xindex] += a_data[aindex] * pow(x_data[xindex], aindex);
      if (pder) {
        // \frac{\dell }{\dell C} C * x^j  = x^j
        pder_data[xindex + aindex * pder_dim[0]] = pow(x_data[xindex], aindex);
      }
    }
  }
  // These might be temp. due to type conversion
  IDL_DELTMP(x);
  IDL_DELTMP(a);
}

/*; Use         : COMP_GAUSS,X,A,F [,PDER]*/
static void COMP_GAUSS(int argc, IDL_VPTR Argv[], char *argk)
{
  info("COMP_GAUSS: Running!");
  check_params(argc, Argv);

  IDL_VPTR x = IDL_CvtDbl(1, Argv);     /* Input array */
  IDL_VPTR a = IDL_CvtDbl(1, Argv + 1); /* Coefficients */

  IDL_VPTR f = Argv[2];                      /* Output array */
  IDL_VPTR pder = argc > 3 ? Argv[3] : NULL; /* Partial derivatives, optional */

  make_arr_from_template(x, f);

  double *x_data = (void *) x->value.arr->data;
  double *a_data = (void *) a->value.arr->data;
  double *f_data = (void *) f->value.arr->data;
  double *pder_data = NULL;

  IDL_MEMINT pder_dim[2];
  if (pder) {
    make_pder_array(x, a, pder, pder_dim);
    pder_data = (void *) pder->value.arr->data;
  }

  double A = a_data[0];       // Coefficient a[0]
  double lambda0 = a_data[1]; // Center a[1]
  double w = a_data[2];       // Width a[2]

  double kern;       // e^((lambda-lambda0)^2/w^2*0.5)
  double lambdadiff; //
  for (IDL_MEMINT xindex = 0; xindex < x->value.arr->n_elts; xindex++) {
    double lambda = x_data[xindex];
    double z = (lambda - lambda0) / w;
    double z2 = z * z;
    double kern = exp(-z2 * 0.5);
    if (z2 > 1000.) {
      kern = 0.0; // Avoid overflow for large z^2
      // info("COMP_GAUSS: Using exp(-0.5 * z^2) for small z^2");
    }
    double F = A * kern;
    f_data[xindex] = F;
    for (IDL_MEMINT aindex = 0; aindex < a->value.arr->n_elts; aindex++) {
      if (pder) {
        if (aindex == 0) {
          pder_data[xindex + aindex * pder_dim[0]] = kern; // Derivative w.r.t. a[0]
        }
        if (aindex == 1) {
          // pder_data[xindex + aindex * pder_dim[0]] = F * z / (w * w);
          pder_data[xindex + aindex * pder_dim[0]] = F * z / w;
        }
        if (aindex == 2) {
          double pder1 = F * z / w;
          // pder_data[xindex + aindex * pder_dim[0]] = F * z2 / (w * w * w);
          pder_data[xindex + aindex * pder_dim[0]] = pder1 * z;
        }
      }
    }
  }
  // These might be temp. due to type conversion
  IDL_DELTMP(x);
  IDL_DELTMP(a);
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
  pder(ix,2) = pder(ix,1) * z(ix) ;; a(0)exp(...) * (x-a(1))^2/a(2)^3
  help,pder

END
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
