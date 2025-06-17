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
  info("COMP_POLY: Running!");
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

  double exp_term;   // e^((lambda-lambda0)^2/w^2*0.5)
  double lambdadiff; //
  for (IDL_MEMINT xindex = 0; xindex < x->value.arr->n_elts; xindex++) {
    double lambda = x_data[xindex];
    double z = (lambda - lambda0) / w; // Standardized variable
    double z2 = z * z;                 // z squared
    double exp_term = exp(-0.5 * z2);  // Exponential term
    double F = A * exp_term;           // Gaussian function value
    f_data[xindex] = F;                // Gaussian function value
    for (IDL_MEMINT aindex = 0; aindex < a->value.arr->n_elts; aindex++) {
      if (pder) {
        if (aindex == 0) {
          // Derivative w.r.t. a[0]
          // \frac{\dell f}{\dell a[0]} = e^{-\frac{(x-a[1])^2}{2a[2]^2}}
          pder_data[xindex + aindex * pder_dim[0]] = exp_term; // Derivative w.r.t. a[0]
        }
        if (aindex == 1) {
          // Derivative w.r.t. a[1]
          // \frac{\dell f}{\dell a[1]} = f * \frac{(x-a[1])}{a[2]^2}
          pder_data[xindex + aindex * pder_dim[0]] = F * z / w;
        }
        if (aindex == 2) {
          // Derivative w.r.t. a[2]
          // \frac{\dell f}{\dell a[2]} = f * \frac{(x-a[1])^2}{a[2]^3}
          pder_data[xindex + aindex * pder_dim[0]] = pder_data[xindex + 1 * pder_dim[0]] * z;
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
  z = ( x - a[1] ) / a(2)
  z2 = z*z

  f = a[0] * exp( -z2 * 0.5 )

  f = a[0] * exp( - (( x - a[1] ) / a[2]) ^2 * 0.5 )
  ;----

  kern = exp( -z2 * 0.5 )

  f = a(0)*kern

  pder[0] = kern         ;;
  pder[1] = f * z/a[2]   ;; a(0)exp(-0.5*(x-a(1))^2/a(2)^2) * (x-a(1))/a(2)^2
  pder[2] = pder[1] * z  ;; a(0)exp(...) * (x-a(1))^2/a(2)^3
 */

int IDL_Load(void)
{
  static IDL_SYSFUN_DEF2 pro_def[] = {{(IDL_SYSRTN_GENERIC) COMP_POLY, "COMP_POLY", 3, 4, 0, 0},
                                      {(IDL_SYSRTN_GENERIC) COMP_GAUSS, "COMP_GAUSS", 3, 4, 0, 0}};
  return IDL_SysRtnAdd(pro_def, FALSE, 2);
}
