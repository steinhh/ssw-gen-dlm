#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <export.h>

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

/* Quicksort-based (partitioning) k'th element finder
   Based on http://www.mathcs.carleton.edu/courses/course_resources/
   cs227_w96/swanjorg/algorithms.html

   INPUTS:
   A - Array with data (N elements)
   K - The number of the element to be found
   N - The number of elements in A
*/

#define SWAP(a, b) (atemp = a, a = b, b = atemp)

#define COMP_POLY_CORE(TYPE, VAR_VALUE_FIELD, ALWAYS)     \
  {                                                       \
    TYPE *in = (TYPE *)a->value.arr->data;                \
    TYPE *out = (TYPE *)res->value.arr->data;             \
    TYPE *workspace = (TYPE *)Workspace->value.arr->data; \
  }

#define NO_STRUCT_DEF

/*; Use         : COMP_POLY,X,A,F [,PDER]*/
static void COMP_POLY(int argc, IDL_VPTR Argv[], char *argk)
{
  IDL_VPTR x = IDL_BasicTypeConversion(1, Argv, IDL_TYP_DOUBLE);     /* Input array */
  IDL_VPTR a = IDL_BasicTypeConversion(1, Argv + 1, IDL_TYP_DOUBLE); /* Coefficients */
  IDL_VPTR f = Argv[2];                                              /* Output array */
  IDL_VPTR pder = argc > 3 ? Argv[3] : NULL;                         /* Partial derivatives, optional */

  IDL_VPTR Workspace;

  /* TYPE CHECKING / ALLOCATION SECTION */
  if (!(x->flags & IDL_V_ARR))
  {
    bailout("X (1st arg) must be an array");
  }

  if (!(a->flags & IDL_V_ARR))
  {
    bailout("A (coefficients) must be an array");
  }
  if (a->value.arr->n_dim != 1)
  {
    bailout("A (coefficients) must be a 1-dimensional array");
  }

  IDL_VarCopy(x, f); /* Copy x to res */

  // x and a may be temp variables due to type conversion, delete if so
  // The IDL_DELTMP macro checks, does not hurt non-temp variables:
  IDL_DELTMP(x);
  IDL_DELTMP(a);

  // IDL_VarMakeTempFromTemplate(x, x->type, NO_STRUCT_DEF, &f, FALSE);

  /* Make workspace: */
  // IDL_MakeTempVector(a->type, b->value.l * c->value.l, IDL_ARR_INI_NOP, &Workspace);

  // if (a != Argv[0])
  //   IDL_DELTMP(a);
  // if (b != Argv[1])
  //   IDL_DELTMP(b);
  // if (c != Argv[2])
  //   IDL_DELTMP(c);
  // IDL_DELTMP(Workspace);

  // return res;
}

int IDL_Load(void)
{
  static IDL_SYSFUN_DEF2 pro_def[] =
      {{{COMP_POLY}, "COMP_POLY", 3, 4, 0, 0}};

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
; Calls       : None.
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