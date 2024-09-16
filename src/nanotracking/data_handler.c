#define PY_SSIZE_T_CLEAN
#include <Python.h>
// #include "numpy/arrayobject.h"

// static PyObject *parse_data(PyObject *self, PyObject *args, PyObject *kwargs) {
static PyObject *data_handler_parse_data(PyObject *self, PyObject *args, PyObject *kwargs) {
// static PyCFunction parse_data(PyObject *self, PyObject *args, PyObject *kwargs) {
    // static char *keywords[] = {"bins", "sizes", "new_file_path", "num_data_points"};
    // PyObject *bins_arg = NULL, *bins = NULL;
    // PyObject *sizes_arg = NULL, *sizes = NULL;
    // const char *new_file_path;
    // int num_data_points;
    // if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!si", keywords, &PyArray_Type, &bins_arg, &PyArray_Type, &sizes_arg, &new_file_path, &num_data_points))
    //     return NULL;
    
    // bins = PyArray_FROM_OTF(bins_arg, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    // if (bins == NULL) goto fail;
    // sizes = PyArray_FROM_OTF(sizes_arg, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    // if (sizes == NULL) goto fail;
    static char *keywords[] = {"something"};
    int something;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i", keywords, &something))
        return NULL;
    
    // FILE *filePtr;
    // FILE *newFilePtr;
    // filePtr = fopen(dat_path, "r");
    // newFilePtr = fopen(new_file_path, "w");
    // for (int i = 0; i< num_data_points; i++) {
        
    //     fgets()
    // }
    // fclose(filePtr);
    
    // assert(PyArray_API);

    Py_INCREF(Py_None);
    return Py_None;
    
    // fail:
    //     Py_XDECREF(bins);
    //     Py_XDECREF(sizes);
    //     return NULL;
}

static PyMethodDef methods[] = {
    // {"parse_data",  parse_data, METH_VARARGS, ""},
    // {"parse_data",  data_handler_parse_data, METH_VARARGS, ""},
    {"parse_data",  data_handler_parse_data, METH_VARARGS | METH_KEYWORDS, ""},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "data_handler",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    methods
};

PyMODINIT_FUNC PyInit_data_handler(void) {
    // Initialise Numpy
    // assert(!PyErr_Occurred());
    // import_array();
    // if (PyErr_Occurred()) {
    //     printf("Failed to import numpy Python module(s).");
    //     return NULL;
    // }
    return PyModule_Create(&module);
}