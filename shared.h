#pragma once
#ifndef EXTERN 
#define EXTERN extern
#endif

// Definition of global variables. This is simply done so one can 
// Avoid passing numbers around.

EXTERN int DEPTH;
EXTERN int MAXDEPTH;
EXTERN double CONV;
EXTERN double TFACTOR;
EXTERN Py_ssize_t DIM;
EXTERN Py_ssize_t LDIM; 
EXTERN Py_ssize_t DATALEN;
EXTERN PyObject* PLIST;
EXTERN PyObject* LATTICE;
EXTERN PyObject* THROTTLE = NULL;
EXTERN char* MODE = NULL;

// Now import itertools.product. We do it here so we only have to 
// import it once rather than once
EXTERN PyObject* ITERTOOLS; 
EXTERN PyObject* PRODUCT;
  
// declaration of function prototypes for ndfit
static inline PyObject* ndfit_getminimum(PyObject* list);
static inline PyObject* ndfit_callfunc(PyObject* func, PyObject* values, PyObject* params);
static PyObject* ndfit_maxdepth(PyObject* self, PyObject* args);
static PyObject* ndfit_dotproduct(PyObject* a, PyObject* b);
static PyObject* ndfit_dotdivide(PyObject* a, PyObject* b);
static PyObject* ndfit_dotadd(PyObject* a, PyObject* b);
static PyObject* ndfit_dotsubtract(PyObject* a, PyObject* b);
static double ndfit_entropy(PyObject* callfunc, PyObject* data, PyObject* params);
static PyObject* ndfit_permutatorshort(PyObject* step, double scale);
static PyObject* ndfit_permutatorfull(PyObject* step, double scale);
static PyObject* ndfit_next(PyObject* callfunc, PyObject* data,PyObject* params,PyObject* lattice);
PyObject* ndfit_recursive(PyObject* callfunc, PyObject* data, PyObject* params, PyObject* step, PyObject* lattice);
PyObject* ndfit_run(PyObject* self,PyObject *args, PyObject *kwds);

// Declaration of helper functions
PyObject* ndfit_product(PyObject* self, PyObject* args);
PyObject* ndfit_quotient(PyObject* self, PyObject* args);
PyObject* ndfit_sum(PyObject* self, PyObject* args);
PyObject* ndfit_difference(PyObject* self, PyObject* args);
PyObject* ndfit_functest(PyObject* self, PyObject* args);

// Declaration of initialization function
PyMODINIT_FUNC initndfit(void);

/////////////////////////////////////////////////////////
// ~~~~~~~~~~ DATA STRUCTURE DEFINITION  ~~~~~~~~~~~~~~//
/////////////////////////////////////////////////////////
// declaration of function prototypes for ndfitstruct

#ifndef NDSTRUCT
#define NDSTRUCT
typedef struct ndFit{
  PyObject_HEAD
  PyObject* data;
  PyObject* pList;
  PyObject* fitfunc;
  PyObject* errfunc;
  PyObject* lattice;
} ndFit;
#endif

PyTypeObject ndFitType;
PyMethodDef ndFit_methods;
PyMemberDef ndFit_members;

void ndFit_dealloc(ndFit* self);
PyObject* ndFit_new(PyTypeObject* type, PyObject* args, PyObject* kwds);




