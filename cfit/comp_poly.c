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

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define SWAP(a, b) (atemp = a, a = b, b = atemp)

#define COMP_POLY_CORE(TYPE, VAR_VALUE_FIELD, ALWAYS)                                                                  \
  {                                                                                                                    \
    TYPE *in = (TYPE *) a->value.arr->data;                                                                            \
    TYPE *out = (TYPE *) res->value.arr->data;                                                                         \
    TYPE *workspace = (TYPE *) Workspace->value.arr->data;                                                             \
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

/*; Use         : COMP_POLY,X,A,F [,PDER]*/
static void COMP_POLY(int argc, IDL_VPTR Argv[], char *argk)
{
  IDL_ENSURE_ARRAY(Argv[0]); /* Ensure X is an array */
  IDL_ENSURE_ARRAY(Argv[1]); /* Ensure A is an array */
  if (Argv[1]->value.arr->n_dim != 1) {
    bailout("A (coefficients) must be a 1-dimensional array");
  }
  assert_numeric("X (1st arg)", Argv[0]);
  assert_numeric("A (coefficients)", Argv[1]);

  IDL_VPTR x = IDL_CvtDbl(1, Argv);     /* Input array */
  IDL_VPTR a = IDL_CvtDbl(1, Argv + 1); /* Coefficients */

  IDL_VPTR f = Argv[2];                      /* Output array */
  IDL_VPTR pder = argc > 3 ? Argv[3] : NULL; /* Partial derivatives, optional */

  IDL_VPTR tmp;

  // Make F output array, first temporary then make permanent with IDL_VarCopy
  IDL_VarMakeTempFromTemplate(x, x->type, NO_STRUCT_DEF, &tmp, FALSE);
  IDL_VarCopy(tmp, f);

  IDL_MEMINT pder_dim[2];
  IDL_MEMINT x_elements = x->value.arr->n_elts;
  IDL_MEMINT a_coefficients = a->value.arr->n_elts;
  if (pder) {
    pder_dim[0] = x->value.arr->n_elts; // Number of elements in x
    pder_dim[1] = a->value.arr->n_elts; // Number of coefficients in a
    IDL_MakeTempArray(IDL_TYP_DOUBLE, 2, pder_dim, IDL_ARR_INI_NOP, &tmp);
    IDL_VarCopy(tmp, pder);
  }

  double *x_data = (void *) x->value.arr->data;
  double *a_data = (void *) a->value.arr->data;
  double *f_data = (void *) f->value.arr->data;
  double *pder_data = NULL;
  if (pder) {
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
      // These might be temp. due to type conversion
      IDL_DELTMP(x);
      IDL_DELTMP(a);
    }
  }
}

int IDL_Load(void)
{
  static IDL_SYSFUN_DEF2 pro_def[] = {{(IDL_SYSRTN_GENERIC) COMP_POLY, "COMP_POLY", 3, 4, 0, 0}};
  return IDL_SysRtnAdd(pro_def, FALSE, 1);
}

/*
; Explanation : Input coefficients A determine degree of polynomial, otherwise
;               this is straightforward - see e.g., CURVEFIT for explanations
;               about this type of function.
;
; Use         : COMP_POLY,X,A,F [,PDER]
;
; Inputs      : As all CURVEFIT evaluation functions
;
; Opt. Inputs : PDER
;
; Outputs     : F : Evaluated function
;
; Opt. Outputs: PDER : Partial derivatives.
;
; Keywords    : None.
;
; Calls       : None
;
; Common      : None.
;
; Restrictions: None.
;
; Side effects: None.
;
; Category    : Analysis
;
; Prev. Hist. : None.
;
; Written     : S.V.H.Haugan, UiO, 21 January 1997
;
; Modified    : Not yet
;
; Version     : 1, 21 January 1997
;-
PRO comp_poly,x,a,f,pder

  f = poly(x,a)

  IF N_params() EQ 4 THEN BEGIN
     nx = n_elements(x)
     nterms = N_elements(a)
     type = datatype(a,2)
     pder = make_array(nx,nterms,type=type,value=1.0)

     ;; Zero-order term

     ;; pder(*,0) = 1.0 ;;Already done..

     ;; First
     IF nterms GT 1 THEN pder(0,1) = x

     ;; Subsequent
     FOR i = 2,nterms-1 DO BEGIN
        pder(0,i) = x*pder(*,i-1)
     END
  END
END

*/