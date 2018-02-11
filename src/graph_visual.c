#include </usr/include/python2.7/Python.h>

int graph_visual(int argc, char *argv[])
{
    PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pListArgs, *pValue;
    int i;
    int array[10] = {1, 3, 5}; int len = 3;

    if (argc < 3) {
        fprintf(stderr,"Usage: call pythonfile funcname [args]\n");
        return 1;
    }

    Py_Initialize();
    PyRun_SimpleString("import sys\nsys.path.append('.')\n");
    pName = PyString_FromString(argv[1]);
    /* Error checking of pName left out */

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, argv[2]);
        /* pFunc is a new reference */
        if (pFunc && PyCallable_Check(pFunc)) {
            { // set up a list argument
                pArgs = PyTuple_New(1);
                pListArgs = PyList_New(len);
                for (i = 0; i < len; ++i) {
                    pValue = PyInt_FromLong(array[i]);
                    PyList_SetItem(pListArgs, i, pValue);
                }
                PyTuple_SetItem(pArgs, 0, pListArgs);
            }
            pValue = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            if (pValue != NULL) {
                printf("Result of call: %ld\n", PyInt_AsLong(pValue));
                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", argv[2]);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", argv[1]);
        return 1;
    }
    Py_Finalize();
    return 0;
}
