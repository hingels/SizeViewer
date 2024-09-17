#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "numpy/arrayobject.h"

static PyObject *parse_data(PyObject *self, PyObject *args, PyObject *kwargs) {
    assert(PyArray_API);
    
    // static char *keywords[] = {"bins", "sizes", "paths", "num_data_points", NULL};
    static char *keywords[] = {"bins", "sizes", "sample_filename", "outputs_path", "num_data_points", NULL};
    // PyObject *bins_arg = NULL;
    // PyObject *sizes_arg = NULL;
    PyArrayObject *bins_arrObj = NULL;
    PyArrayObject *sizes_arrObj = NULL;
    // PyObject *paths_arg = NULL;//, *paths_arrObj = NULL;
    char *sample_filename = NULL;
    char *outputs_path = NULL;
    // int num_data_points = 0;
    int num_data_points = 0;
    // if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!i", keywords, &PyArray_Type, &bins_arg, &PyArray_Type, &sizes_arg, &PyDict_Type, &paths_arg, &num_data_points))
    // if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!ssi", keywords, &PyArray_Type, &bins_arg, &PyArray_Type, &sizes_arg, &sample_filename, &outputs_path, &num_data_points))
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!ssi", keywords, &PyArray_Type, &bins_arrObj, &PyArray_Type, &sizes_arrObj, &sample_filename, &outputs_path, &num_data_points))
        return NULL;
    
    if (!PyArray_Check(bins_arrObj) || !PyArray_ISNUMBER((PyArrayObject *) bins_arrObj)) {       // https://stackoverflow.com/a/72932084
        PyErr_SetString(PyExc_TypeError, "Bins argument must be a NumPy array with numeric dtype.");
        return NULL;
    }
    if (!PyArray_Check(sizes_arrObj) || !PyArray_ISNUMBER((PyArrayObject *) sizes_arrObj)) {       // https://stackoverflow.com/a/72932084
        PyErr_SetString(PyExc_TypeError, "Sizes argument must be a NumPy array with numeric dtype.");
        return NULL;
    }
    
    // npy_double bins[num_data_points];
    // npy_double sizes[num_data_points];
    // size_t data_size = sizeof(npy_double);
    double bins[num_data_points];
    double sizes[num_data_points];
    size_t data_size = sizeof(double);
    // npy_double *bins, *sizes;
    PySys_WriteStdout("Size %zu", data_size);
    // bins = malloc(num_data_points*data_size);
    // sizes = malloc(num_data_points*data_size);
    
    // PyObject *bins_arrObj = PyArray_FROM_OTF(bins_arg, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (bins_arrObj == NULL) goto fail;
    // PyObject *sizes_arrObj = PyArray_FROM_OTF(sizes_arg, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (sizes_arrObj == NULL) goto fail;

    // npy_double* bins[num_data_points];
    // npy_double* sizes[num_data_points];
    // npy_double* bins;
    // npy_double* sizes;
    // double* bins;
    // double* sizes;

    for (int i = 0; i < num_data_points; i++) {
        PySys_WriteStdout("Current bin %f value", (double) bins[i]);
    }
    PySys_WriteStdout("\n");
    bins[0] = 12345.0;
    
    if (PyArray_NDIM(bins_arrObj) != 1 || PyArray_NDIM(sizes_arrObj) != 1) goto fail;
    if ((int) *(PyArray_DIMS(bins_arrObj)) != num_data_points || (int) *(PyArray_DIMS(sizes_arrObj)) != num_data_points) goto fail;
    // if (PyArray_STRIDE(bins_arrObj) != …………… ||  ………) goto fail;
    // npy_double *bins = (npy_double *) PyArray_DATA(bins_arrObj);
    // npy_double *sizes = (npy_double *) PyArray_DATA(sizes_arrObj);
    // double *bins = (double *) PyArray_DATA(bins_arrObj);
    // double *sizes = (double *) PyArray_DATA(sizes_arrObj);
    // char *bins = PyArray_BYTES(bins_arrObj);
    // char *sizes = PyArray_BYTES(sizes_arrObj);
    // npy_intp dims[1] = {(npy_intp) num_data_points};
    npy_intp dims[1];
    // int ndim;
    PyArray_Descr *dtype = PyArray_DescrFromType(NPY_DOUBLE);
    Py_INCREF(bins_arrObj); // PyArray_AsCArray steals reference
    Py_INCREF(sizes_arrObj); // PyArray_AsCArray steals reference
    // if (PyArray_AsCArray((PyObject **) &bins_arrObj, (void *) &bins, dims, 1, dtype) != 0) {
    if (PyArray_AsCArray(&bins_arrObj, &bins, dims, 1, dtype) != 0) {
        PyErr_SetString(PyExc_TypeError, "Error converting to C array.");
        goto fail;
    }
    // if (PyArray_AsCArray((PyObject **) &sizes_arrObj, (void *) &sizes, dims, 1, dtype) != 0) {
    if (PyArray_AsCArray(&sizes_arrObj, &sizes, dims, 1, dtype) != 0) {
        PyErr_SetString(PyExc_TypeError, "Error converting to C array.");
        goto fail;
    }

    // FILE *bins_file = fopen((char *) PyDict_GetItemString(paths_arg, "bins"), "wb");
    // PyObject *bins_path = PyObject_Str(PyDict_GetItemString(paths_arg, "FAKE"));
    // char bins_path[100]; snprintf(bins_path, 100, "%s/%s", outputs_path, "bins");
    char filepath[100]; snprintf(filepath, 100, "%s/%s", outputs_path, sample_filename);
    FILE *file = fopen(filepath, "wb");
    if (file == NULL) {
        fclose(file);
        goto fail;
    }
    // size_t num_written = fwrite(bins, sizeof(npy_double), num_data_points, file);
    
    fwrite(bins, data_size, num_data_points, file);
    fwrite(sizes, data_size, num_data_points, file);
    // fwrite(bins, sizeof(double), num_data_points, file);
    // fwrite(sizes, sizeof(double), num_data_points, file);
    for (int i = 0; i < num_data_points; i++) {
        PySys_WriteStdout("Writing bin %f to file", (double) bins[i]);
    }
    fclose(file);
    PyArray_Free(bins_arrObj, bins);
    PyArray_Free(sizes_arrObj, sizes);
    // free(bins);
    // free(sizes);
    Py_RETURN_NONE;
    
    fail:
        Py_XDECREF(bins_arrObj);
        Py_XDECREF(sizes_arrObj);
        PyArray_Free(bins_arrObj, bins);
        PyArray_Free(sizes_arrObj, sizes);
        // free(bins);
        // free(sizes);
        return NULL;
}
static PyObject *read_data(PyObject *self, PyObject *args, PyObject *kwargs) {
    assert(PyArray_API);
    
    static char *keywords[] = {"sample_filename", "outputs_path", "num_data_points", NULL};
    char *sample_filename = NULL;
    char *outputs_path = NULL;
    int num_data_points = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ssi", keywords, &sample_filename, &outputs_path, &num_data_points))
        return NULL;

    char filepath[100]; snprintf(filepath, 100, "%s/%s", outputs_path, sample_filename);
    npy_double bins[num_data_points];
    npy_double sizes[num_data_points];
    FILE *file = fopen(filepath, "rb");
    if (file == NULL) {
        fclose(file);
        goto fail;
    }
    fread(&bins, sizeof(npy_double), num_data_points, file);
    fread(&sizes, sizeof(npy_double), num_data_points, file);
    fclose(file);
    npy_intp dims[1] = {(npy_intp) num_data_points};
    PyObject *bin_arr = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, bins);
    PyObject *size_arr = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, sizes);
    Py_INCREF(bin_arr);
    Py_INCREF(size_arr);
    PyObject *tuple = PyTuple_Pack(2, bin_arr, size_arr);
    Py_INCREF(tuple);
    return tuple;
    
    fail:
        return NULL;
}

static PyMethodDef methods[] = {
    // {"parse_data",  parse_data, METH_VARARGS | METH_KEYWORDS, ""},
    {"parse_data",  (PyCFunction)(void(*)(void)) parse_data, METH_VARARGS | METH_KEYWORDS, ""},
    {"read_data",  (PyCFunction)(void(*)(void)) read_data, METH_VARARGS | METH_KEYWORDS, ""},
    {NULL, NULL, 0, NULL}           // Sentinel
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "data_handler",     // Module name
    NULL,               // Module documentation
    -1,       /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    methods
};

PyMODINIT_FUNC PyInit_data_handler(void) {
    assert(!PyErr_Occurred());
    import_array(); // Initialize NumPy
    if (PyErr_Occurred()) {
        printf("Failed to import NumPy.");
        return NULL;
    }

    return PyModule_Create(&module);
}