//An N-dimensional curve fitting tool written in C Python
//Copyright (C) 2013  Michael Winters : micwinte@chalmers.se

//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; either version 2
//of the License, or (at your option) any later version.

//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

// Include Libraries
#include <Python.h>
#include <structmember.h>
#include <math.h>

#define EXTERN
#include "shared.h"

// Python exception obj
static PyObject* ndfitError;

//////////////////////
// Helper Functions //
//////////////////////
// Find minimum and calculate function methods
static inline PyObject* ndfit_getminimum(PyObject* list){
  PyList_Sort(list); return PyList_GetItem(list,0);}

static inline PyObject* ndfit_callfunc(PyObject* func, PyObject* values, PyObject* params){
  return PyObject_CallFunctionObjArgs(func,values,params,NULL);}

// Setters for maxdepth and convergence
static PyObject* ndfit_maxdepth(PyObject* self, PyObject* args){
  if(!PyArg_ParseTuple(args,"i",&MAXDEPTH)){return NULL;}
  Py_RETURN_NONE;
}
static PyObject* ndfit_convergence(PyObject* self, PyObject* args){
  if(!PyArg_ParseTuple(args,"d",&CONV)){return NULL;}
  Py_RETURN_NONE;
}
static PyObject* ndfit_throttle_factor(PyObject* self, PyObject* args){
  if(!PyArg_ParseTuple(args,"d",&TFACTOR)){return NULL;}
  Py_RETURN_NONE;
}

///////////////////////////////////////////
// List algebra methods needed for ndfit //
///////////////////////////////////////////
static PyObject* ndfit_dotproduct(PyObject* a, PyObject* b){
  double tmp;
  Py_ssize_t i; 
  Py_ssize_t size = PyList_Size(a);
  PyObject* product = PyList_New(size);

  for (i=0;i<size;i+=1){
      tmp = PyFloat_AsDouble(PyList_GetItem(a,i))*PyFloat_AsDouble(PyList_GetItem(b,i));
      PyList_SetItem(product,i,Py_BuildValue("d",tmp));
  }
  return product;
}

static PyObject* ndfit_dotadd(PyObject* a, PyObject* b){
  double tmp;
  Py_ssize_t i; 
  Py_ssize_t size = PyList_Size(a);
  PyObject* sum = PyList_New(size);

  for (i=0;i<size;i+=1){
      tmp = PyFloat_AsDouble(PyList_GetItem(a,i))+PyFloat_AsDouble(PyList_GetItem(b,i));
      PyList_SetItem(sum,i,Py_BuildValue("d",tmp));
  }
  return sum;
}

////////////////////////////////
// Entropy Calculation Method //
////////////////////////////////
static double ndfit_entropy(PyObject* callfunc, PyObject* data, PyObject* params){
  
  // Initialize counter
  Py_ssize_t i;

  // Initialize Residual List and tmp object (will not be returned)
  PyObject* values;
  PyObject* res_list;
  res_list = PyList_New(DATALEN);

  // Perform Loop to calculate residuals
  for(i=0;i<DATALEN;i+=1){
    Py_INCREF(data);
    values = ndfit_callfunc(callfunc,PyList_GetItem(data,i),params);
    PyList_SetItem(res_list,i,values);
    Py_DECREF(data);
  }
    
  // Sum over the residual list to get a final value
  double tmp;
  double sum = 0.0;
  for(i=0;i<DATALEN;i+=1){
    Py_INCREF(res_list);
    tmp = PyFloat_AsDouble(PyList_GetItem(res_list,i));
    sum+= tmp*tmp;
    Py_DECREF(res_list);
  }
  Py_DECREF(res_list);
  // Decref to things we are done with
  return (sqrt(sum)*(double)log((double)DATALEN))/((double)DATALEN);
};

//////////////////////////
// Build Lattice Method //
//////////////////////////
// A static method to build the permutator. This defines the fitting
// lattice. (e.g. the cube corners: +++,-++,+-+,++-,--+,-+-,+--,---)
// Note that the code only calls this method once.
static PyObject* ndfit_permutatorshort(PyObject* step, double scale){

  // Build [1,-1] for perumutation this is what 
  // will be fed to itertools
  PyObject* pm = Py_BuildValue("(d,d)",(1.0*scale)/sqrt(DIM),(-1.0*scale)/sqrt(DIM));

  // Pack these into the arglist. We need m of them where
  // m is the number of fitting parameters
  PyObject* args = PyTuple_New(DIM);

  Py_ssize_t i = 0;
  for(i=0; i<(Py_ssize_t)DIM;i+=1){
    PyTuple_SetItem(args,i,pm);
    Py_INCREF(pm);
  }
  PyObject* iterator = PyObject_CallObject(PRODUCT, args);
  
  // Result is an iterator which returns a tuple. We want to turn this 
  // into a static list of tuples which can be saved and used foever. 
  double tmp;
  Py_ssize_t j;
     
  PyObject* item;
  PyObject* lattice = PyList_New((Py_ssize_t)pow(2,DIM));
  // Get the (2^D) corners for the lattice
  for (i=0; i<(Py_ssize_t)pow(2,DIM); i+=1){
    PyObject* corner = PyList_New((Py_ssize_t)DIM);
    item = PyIter_Next(iterator);
    for (j=0; j<(Py_ssize_t)DIM; j+=1){
      tmp = PyFloat_AsDouble(PyTuple_GetItem(item,j));
      PyList_SetItem(corner,j,Py_BuildValue("d",tmp));
      Py_INCREF(corner);
    }
    PyList_SetItem(lattice,i,corner);
    Py_DECREF(item);
    Py_DECREF(corner);
  }
  Py_DECREF(iterator);
  
  for (i=0; i<(Py_ssize_t)pow(2,DIM); i+=1){
    PyList_SetItem(lattice,i,ndfit_dotproduct(step,PyList_GetItem(lattice,i)));
  } 

  // Clean up
  Py_DECREF(pm);
  Py_DECREF(args);
  LDIM = (Py_ssize_t)pow(2,DIM);
  return lattice;
}

static PyObject* ndfit_permutatorfull(PyObject* step, double scale){

    // Build [1,-1] for perumutation this is what 
  // will be fed to itertools
  PyObject* pm = Py_BuildValue("(d,d)",(1.0*scale)/sqrt(DIM),(-1.0*scale)/sqrt(DIM));

  // Pack these into the arglist. We need m of them where
  // m is the number of fitting parameters
  PyObject* args = PyTuple_New(DIM);

  Py_ssize_t i = 0;
  for(i=0; i<(Py_ssize_t)DIM;i+=1){
    PyTuple_SetItem(args,i,pm);
    Py_INCREF(pm);
  }
  PyObject* iterator = PyObject_CallObject(PRODUCT, args);
  
  // Result is an iterator which returns a tuple. We want to turn this 
  // into a static list of tuples which can be saved and used foever. 
  double tmp;
  Py_ssize_t j;
     
  PyObject* item;
  PyObject* lattice = PyList_New((Py_ssize_t)pow(2,DIM)+2*DIM);
  // Get the (2^D) corners for the lattice
  for (i=0; i<(Py_ssize_t)pow(2,DIM); i+=1){
    PyObject* corner = PyList_New((Py_ssize_t)DIM);
    item = PyIter_Next(iterator);
    for (j=0; j<(Py_ssize_t)DIM; j+=1){
      tmp = PyFloat_AsDouble(PyTuple_GetItem(item,j));
      PyList_SetItem(corner,j,Py_BuildValue("d",tmp));
      Py_INCREF(corner);
    }
    PyList_SetItem(lattice,i,corner);
    Py_DECREF(item);
    Py_DECREF(corner);
  }
  Py_DECREF(iterator);
  
  for (i=0; i<(Py_ssize_t)pow(2,DIM); i+=1){
    PyList_SetItem(lattice,i,ndfit_dotproduct(step,PyList_GetItem(lattice,i)));
  } 

  // Get the (2*DIM) edges for the lattice not that the numer of corners. 
  for (i=0; i<(Py_ssize_t)DIM; i+=1){
    PyObject* args2 = PyList_New(DIM);
    for (j=0; j<DIM; j+=1){PyList_SetItem(args2,j,Py_BuildValue("d",(double)0));} 
    PyList_SetItem(args2,i,Py_BuildValue("d",(double)scale));
    PyList_SetItem(lattice,i+pow(2,DIM),args2);
  }
  // Another round for the negative sides
  for (i=0; i<(Py_ssize_t)DIM; i+=1){
    PyObject* args2 = PyList_New(DIM);
    for (j=0; j<DIM; j+=1){PyList_SetItem(args2,j,Py_BuildValue("d",(double)0));} 
    PyList_SetItem(args2,i,Py_BuildValue("d",(double)(-1*scale)));
    PyList_SetItem(lattice,i+pow(2,DIM)+DIM,args2);
  }

  LDIM = (Py_ssize_t)(pow(2,DIM)+(2*DIM));
  return lattice;
}

//////////////////////////////////////////////////
// A method to calculate the recursive step one //
//////////////////////////////////////////////////
static PyObject* ndfit_next(PyObject* callfunc,PyObject* data,PyObject* params,PyObject* lattice){
  
  Py_INCREF(callfunc); 
  Py_INCREF(data);
  Py_INCREF(params);
  Py_INCREF(lattice);

  Py_ssize_t i = 0;
  Py_ssize_t lsize = PyList_Size(lattice); 
  PyObject* calc = PyList_New(lsize);

  PyObject* point;
  PyObject* entropy;
    
  for (i=0;i<lsize;i+=1){
    PyObject* tmp = PyTuple_New((Py_ssize_t)2);
 
    point = ndfit_dotadd(PyList_GetItem(lattice,i),params);
    entropy = Py_BuildValue("d",ndfit_entropy(callfunc, data, point));

    Py_INCREF(entropy);
    Py_INCREF(point);
               
    PyTuple_SetItem(tmp,0,entropy);
    PyTuple_SetItem(tmp,1,point);
    Py_INCREF(tmp);
    Py_INCREF(tmp);
     
    PyList_SetItem(calc,i,tmp);
    Py_INCREF(tmp);
  }

  Py_DECREF(callfunc); 
  Py_DECREF(data);
  Py_DECREF(params);
  Py_DECREF(lattice);

  return ndfit_getminimum(calc);
}

PyObject* ndfit_recursive(PyObject* callfunc, PyObject* data, PyObject* params, 
			  PyObject* step,PyObject* lattice){

  PyObject* next = ndfit_next(callfunc,data,params,lattice);
  PyList_SetItem(PLIST,DEPTH,next);
  DEPTH+=1;
  
  PyObject* tmp = PyList_GetItem(PLIST,DEPTH-2);
  Py_INCREF(PLIST);
  double entropy = PyFloat_AsDouble(PyTuple_GetItem(next,0));
  double check = PyFloat_AsDouble(PyTuple_GetItem(tmp,0));
  
  // Scale the lattice appropriately if throrrling is on
  if (entropy > CONV && (int)PyObject_IsTrue(THROTTLE)){ 
    if(!strcmp(MODE,"short")){
      LATTICE = ndfit_permutatorshort(step,((TFACTOR*entropy)+1.0));
    }
    else if(!strcmp(MODE,"full")){
      LATTICE = ndfit_permutatorfull(step,((TFACTOR*entropy)+1.0));
    }
    else{
      LATTICE = ndfit_permutatorshort(step,(TFACTOR*entropy)+1.0);
    }
  }
  // It we are throttling, then we would like to turn it off below 
  // the convergence to prevent overscaling
  else if ((int)PyObject_IsTrue(THROTTLE)){
    if(!strcmp(MODE,"short")){
      LATTICE = ndfit_permutatorshort(step,1.0);
    }
    else if(!strcmp(MODE,"full")){
      LATTICE = ndfit_permutatorfull(step,1.0);
    }
    else{
      LATTICE = ndfit_permutatorshort(step,1.0);
    }
  }
  // Stop Case 1: We arrived at the desired value 
  if (entropy<CONV && entropy>check){
    printf("Recursion Depth: %d\n",DEPTH);
    printf("Fit Entropy %f\n",check);
    return params;  
  } 
  // Stop Case 2: We have hit the maximim recursion depth
  else if(DEPTH==MAXDEPTH){
    printf("Exceeded Maximum Number of Recusive Steps %d\n",MAXDEPTH);
    printf("Fit Entropy: %f\n", check);
    return params;  
  }
  // Otherwise make the tail recursive call
  else{
    return ndfit_recursive(callfunc,data,PyTuple_GetItem(next,1),step,LATTICE);
  }
}

///////////////////////////////////////////////
// Main Method Runs Fit and Optimizes Params //
///////////////////////////////////////////////
PyObject* ndfit_run(PyObject* self,PyObject *args, PyObject *kwds){

  // Input Data we are reading from 
  PyObject* fitfunc;
  PyObject* callfunc; 
  PyObject* data;
  PyObject* params;
  PyObject* step;

  static char *kwlist[] = {"fitfunc","errfunc","data","params","step","mode","throttle",NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOOO|sO", kwlist, 
				   &fitfunc,&callfunc,&data,&params,
				   &step,&MODE,&THROTTLE)){return NULL;}

  // Read in the data and errr check the input
  if(!PyList_Check(data)){
    PyErr_SetString(ndfitError,"Data is not a list");
    return NULL;
  }
  if(!PyList_Check(params)){
    PyErr_SetString(ndfitError,"Params is not a list");
    return NULL;
  }
  if(!PyList_Check(step)){
    PyErr_SetString(ndfitError,"Step is not a list");
    return NULL;
  }
 // Check for lambdas
  if(!PyCallable_Check(fitfunc)){
    PyErr_SetString(ndfitError,"Invalid Fit Function");
    return NULL;   
  }
  if(!PyCallable_Check(callfunc)){
    PyErr_SetString(ndfitError,"Invalid Error Function");
    return NULL;   
  }

  // Increase ref counts if we didnt bail
  Py_INCREF(fitfunc);
  Py_INCREF(callfunc);
  Py_INCREF(data);
  Py_INCREF(params);


  // Test to see if the error function is even callable with the 
  // data provided. This prevents a segmentation fault with bad functions
  if(!ndfit_callfunc(callfunc,PyList_GetItem(data,0),params)){
     PyErr_SetString(ndfitError,"Unable to call error function. Check that input matches data");
     return NULL;
  }

  // Check if the throttling parameter has been set. 
  // If not, then set it to FALSE
  if (THROTTLE==NULL){
    THROTTLE = Py_False;
  }
  Py_INCREF(THROTTLE);

  // Set max depth and convergence to default values if not set already
  if(!CONV){CONV = 0.1;}
  if(!MAXDEPTH){MAXDEPTH = 1000;}
  if(!TFACTOR){TFACTOR = 1.0;}

  // Initialize the necessary gobal parameters based on data sets
  DEPTH = 0;
  DIM = PyList_Size(params);
  DATALEN = PyList_Size(data); 

  // Build the lattice and call the recursive code
  PLIST = PyList_New(MAXDEPTH);

  // Build the lattice from step. First we need to 
  // get itertools.product.
  ITERTOOLS = PyImport_ImportModule((char*)"itertools"); 
  PRODUCT = PyObject_GetAttrString(ITERTOOLS,(char*)"product");
  
  if(!strcmp(MODE,"short")){
    LATTICE = ndfit_permutatorshort(step,1.0);
  }
  else if(!strcmp(MODE,"full")){
    LATTICE = ndfit_permutatorfull(step,1.0);
  }
  else {
    LATTICE = ndfit_permutatorshort(step,1.0);
  }

  // Need to give the plist an initial value
  PyObject* first = ndfit_next(callfunc,data,params,LATTICE);
  PyList_SetItem(PLIST,DEPTH,first);
  Py_INCREF(first);
  DEPTH+=1;

  PyObject* result = ndfit_recursive(callfunc,data,params,step,LATTICE);

  // Build the final values
  PLIST = PyList_GetSlice(PLIST,0,DEPTH-1);
  PyObject* argList = Py_BuildValue("OOOOO",data,PLIST,fitfunc,callfunc,LATTICE);
  PyObject* ndfobj = PyObject_CallObject((PyObject*)&ndFitType,argList);
  Py_DECREF(argList);

  // Clean up
  Py_DECREF(fitfunc);
  Py_DECREF(callfunc);
  Py_DECREF(data);
  Py_DECREF(params);
  Py_DECREF(result);
  Py_DECREF(LATTICE);
  Py_DECREF(ITERTOOLS);
  Py_DECREF(PRODUCT);

  return ndfobj;
};

///////////////////////////
// Misc Useful Functions //
///////////////////////////

PyObject* ndfit_product(PyObject* self, PyObject* args){

  PyObject* a;
  PyObject* b;

  // Import lambda object
  if(!PyArg_ParseTuple(args,"OO",&a,&b)){return NULL;}
  if(!PyList_Check(a)){
    PyErr_SetString(ndfitError,"Arguments must be lists");
    return NULL;
  }
  if(!PyList_Check(b)){
    PyErr_SetString(ndfitError,"Arguments must be lists");
    return NULL;
  }
  
  Py_ssize_t sizea = PyList_Size(a);
  Py_ssize_t sizeb = PyList_Size(b);

  if((int)sizea != (int)sizeb){
     PyErr_SetString(ndfitError,"Lists must be the same size");
     return NULL;
  }

  double tmp;
  Py_ssize_t i; 
  PyObject* result = PyList_New(sizea);

  for (i=0;i<sizea;i+=1){
      tmp = PyFloat_AsDouble(PyList_GetItem(a,i))*PyFloat_AsDouble(PyList_GetItem(b,i));
      PyList_SetItem(result,i,Py_BuildValue("d",tmp));
  }
  return result;
}

PyObject* ndfit_quotient(PyObject* self, PyObject* args){

  PyObject* a;
  PyObject* b;

  // Import lambda object
  if(!PyArg_ParseTuple(args,"OO",&a,&b)){return NULL;}
  if(!PyList_Check(a)){
    PyErr_SetString(ndfitError,"Arguments must be lists");
    return NULL;
  }
  if(!PyList_Check(b)){
    PyErr_SetString(ndfitError,"Arguments must be lists");
    return NULL;
  }
  
  Py_ssize_t sizea = PyList_Size(a);
  Py_ssize_t sizeb = PyList_Size(b);

  if((int)sizea != (int)sizeb){
     PyErr_SetString(ndfitError,"Lists must be the same size");
     return NULL;
  }

  double tmp;
  Py_ssize_t i; 
  PyObject* result = PyList_New(sizea);

  for (i=0;i<sizea;i+=1){
      tmp = PyFloat_AsDouble(PyList_GetItem(a,i))/PyFloat_AsDouble(PyList_GetItem(b,i));
      PyList_SetItem(result,i,Py_BuildValue("d",tmp));
  }
  return result;
}

PyObject* ndfit_sum(PyObject* self, PyObject* args){

  PyObject* a;
  PyObject* b;

  // Import lambda object
  if(!PyArg_ParseTuple(args,"OO",&a,&b)){return NULL;}
  if(!PyList_Check(a)){
    PyErr_SetString(ndfitError,"Arguments must be lists");
    return NULL;
  }
  if(!PyList_Check(b)){
    PyErr_SetString(ndfitError,"Arguments must be lists");
    return NULL;
  }
  
  Py_ssize_t sizea = PyList_Size(a);
  Py_ssize_t sizeb = PyList_Size(b);

  if((int)sizea != (int)sizeb){
     PyErr_SetString(ndfitError,"Lists must be the same size");
     return NULL;
  }

  double tmp;
  Py_ssize_t i; 
  PyObject* result = PyList_New(sizea);

  for (i=0;i<sizea;i+=1){
      tmp = PyFloat_AsDouble(PyList_GetItem(a,i))+PyFloat_AsDouble(PyList_GetItem(b,i));
      PyList_SetItem(result,i,Py_BuildValue("d",tmp));
  }
  return result;
}

PyObject* ndfit_difference(PyObject* self, PyObject* args){

  PyObject* a;
  PyObject* b;

  // Import lambda object
  if(!PyArg_ParseTuple(args,"OO",&a,&b)){return NULL;}
  if(!PyList_Check(a)){
    PyErr_SetString(ndfitError,"Arguments must be lists");
    return NULL;
  }
  if(!PyList_Check(b)){
    PyErr_SetString(ndfitError,"Arguments must be lists");
    return NULL;
  }
  
  Py_ssize_t sizea = PyList_Size(a);
  Py_ssize_t sizeb = PyList_Size(b);

  if((int)sizea != (int)sizeb){
     PyErr_SetString(ndfitError,"Lists must be the same size");
     return NULL;
  }
   double tmp;
  Py_ssize_t i; 
  PyObject* result = PyList_New(sizea);

  for (i=0;i<sizea;i+=1){
      tmp = PyFloat_AsDouble(PyList_GetItem(a,i))-PyFloat_AsDouble(PyList_GetItem(b,i));
      PyList_SetItem(result,i,Py_BuildValue("d",tmp));
  }
  return result;
}

PyObject* ndfit_derivative(PyObject* self, PyObject* args){

  PyObject* a;
  PyObject* b;

  // Import lambda object
  if(!PyArg_ParseTuple(args,"OO",&a,&b)){return NULL;}
  if(!PyList_Check(a)){
    PyErr_SetString(ndfitError,"Arguments must be lists");
    return NULL;
  }
  if(!PyList_Check(b)){
    PyErr_SetString(ndfitError,"Arguments must be lists");
    return NULL;
  }
  
  Py_ssize_t sizea = PyList_Size(a);
  Py_ssize_t sizeb = PyList_Size(b);

  Py_INCREF(a);
  Py_INCREF(b);

  if((int)sizea != (int)sizeb){
     PyErr_SetString(ndfitError,"Lists must be the same size");
     return NULL;
  }
  double deltay;
  double deltax;
  Py_ssize_t i;
  PyObject* result = PyList_New(sizea);
  for (i=0;i<sizea;i+=1){
    // Take care of special cases where we are at the beginning 
    // and end of the data set
    if (i == 0){
      deltay = PyFloat_AsDouble(PyList_GetItem(a,i+1))-PyFloat_AsDouble(PyList_GetItem(a,i));
      deltay+=deltay;
      deltax = PyFloat_AsDouble(PyList_GetItem(b,i+1))-PyFloat_AsDouble(PyList_GetItem(b,i));
      deltax+=deltax;
      PyList_SetItem(result,i,Py_BuildValue("d",(deltay)/(deltax)));
    }
    // Note these must be else if to avoid passing into the else
    // block with i==0,1
    else if (i == 1){
      deltay = PyFloat_AsDouble(PyList_GetItem(a,i+1))-PyFloat_AsDouble(PyList_GetItem(a,i));
      deltay+= PyFloat_AsDouble(PyList_GetItem(a,i))-PyFloat_AsDouble(PyList_GetItem(a,i-1));
      deltax = PyFloat_AsDouble(PyList_GetItem(b,i+1))-PyFloat_AsDouble(PyList_GetItem(b,i));
      deltax+= PyFloat_AsDouble(PyList_GetItem(b,i))-PyFloat_AsDouble(PyList_GetItem(b,i-1));
      PyList_SetItem(result,i,Py_BuildValue("d",(deltay)/(deltax)));
    }
    else if  (i == (sizea-2)){
      deltay = PyFloat_AsDouble(PyList_GetItem(a,i+1))-PyFloat_AsDouble(PyList_GetItem(a,i));
      deltay+= PyFloat_AsDouble(PyList_GetItem(a,i))-PyFloat_AsDouble(PyList_GetItem(a,i-1));
      deltax = PyFloat_AsDouble(PyList_GetItem(b,i+1))-PyFloat_AsDouble(PyList_GetItem(b,i));
      deltax+= PyFloat_AsDouble(PyList_GetItem(b,i))-PyFloat_AsDouble(PyList_GetItem(b,i-1));
      PyList_SetItem(result,i,Py_BuildValue("d",(deltay)/(deltax)));
    }
    else if (i == (sizea-1)){
      deltay = PyFloat_AsDouble(PyList_GetItem(a,i))-PyFloat_AsDouble(PyList_GetItem(a,i-1));
      deltay+=deltay;
      deltax = PyFloat_AsDouble(PyList_GetItem(b,i))-PyFloat_AsDouble(PyList_GetItem(b,i-1));
      deltax+=deltax;
      PyList_SetItem(result,i,Py_BuildValue("d",(deltay)/(deltax)));
    }
    // deltax is always positive so we do not need to worry about 
    // dividing by zero.
    else {
      deltay = PyFloat_AsDouble(PyList_GetItem(a,i+2))-PyFloat_AsDouble(PyList_GetItem(a,i+1));
      deltay+= PyFloat_AsDouble(PyList_GetItem(a,i+1))-PyFloat_AsDouble(PyList_GetItem(a,i));
      deltay+= PyFloat_AsDouble(PyList_GetItem(a,i))-PyFloat_AsDouble(PyList_GetItem(a,i-1));
      deltay+= PyFloat_AsDouble(PyList_GetItem(a,i-1))-PyFloat_AsDouble(PyList_GetItem(a,i-2));

      deltax = PyFloat_AsDouble(PyList_GetItem(b,i+2))-PyFloat_AsDouble(PyList_GetItem(b,i+1));
      deltax+= PyFloat_AsDouble(PyList_GetItem(b,i+1))-PyFloat_AsDouble(PyList_GetItem(b,i));
      deltax+= PyFloat_AsDouble(PyList_GetItem(b,i))-PyFloat_AsDouble(PyList_GetItem(b,i-1));
      deltax+= PyFloat_AsDouble(PyList_GetItem(b,i-1))-PyFloat_AsDouble(PyList_GetItem(b,i-2));

      PyList_SetItem(result,i,Py_BuildValue("d",(deltay)/(deltax)));
    }
  }
  return result;
}

//Function to test lambda (not static = callable from elsewhere)
PyObject* ndfit_functest(PyObject* self, PyObject* args){

  PyObject* func;
  PyObject* point;
  PyObject* params;
  PyObject* result;
  
  // Import lambda object
  if(!PyArg_ParseTuple(args,"OOO",&func,&point,&params)){return NULL;}

  if(!PyCallable_Check(func)){
    PyErr_SetString(ndfitError,"Invalid Function");
    return NULL;   
  }
  else{result = PyObject_CallFunctionObjArgs(func,point,params,NULL);}
  
  if(!result){
    PyErr_SetString(ndfitError,"Unable to call function: check input parameters");
    return NULL;
  }

  Py_INCREF(result);
  Py_DECREF(func);
  Py_DECREF(point);
  Py_DECREF(params);
  return result;
};

////////////////////
// Function Table //
////////////////////

// Here is where we define the methods to appear in the module 
// Internal methods should not be listed here (like sum). This
// way they will not be available in the interpreter.
static PyMethodDef ndfit_methods[]={
  {"maxdepth", ndfit_maxdepth,METH_VARARGS,"set max recursion depth"},
  {"convergence", ndfit_convergence,METH_VARARGS,"set entropy convergence"},
  {"throttle_factor",ndfit_throttle_factor, METH_VARARGS,"set throttle factor"},
  {"run",ndfit_run,METH_VARARGS | METH_KEYWORDS,"main method"},
  {"evaluate_function",ndfit_functest, METH_VARARGS, "external method to check the function"},
  {"product",ndfit_product, METH_VARARGS, "external method to get elementwise product"},
  {"quotient",ndfit_quotient, METH_VARARGS, "external method to get elementwise quotient"},
  {"sum",ndfit_sum, METH_VARARGS, "external method to get elementwise sum"},
  {"difference",ndfit_difference, METH_VARARGS, "external method to get elementwise difference"},
  {"derivative",ndfit_derivative, METH_VARARGS,"external method to calculate the discrete derivitive"},
  {NULL,NULL,0,NULL}, /* Sentinel */
};

/////////////////////////////
// Initialization function //
/////////////////////////////
PyMODINIT_FUNC initndfit(void){
  PyObject *m;
  m = Py_InitModule("ndfit",ndfit_methods);
  if(m==NULL) return;

  // see first line
  ndfitError = PyErr_NewException("ndfit.error",NULL,NULL);
  Py_INCREF(ndfitError);
  PyModule_AddObject(m,"error",ndfitError);

  if (PyType_Ready(&ndFitType) < 0){return;}
  Py_INCREF(&ndFitType);
  PyModule_AddObject(m,"ndFit",(PyObject*)&ndFitType);
};
