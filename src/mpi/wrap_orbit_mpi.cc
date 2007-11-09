#include "Python.h"
#include "orbit_mpi.hh"

#include <cstring>
#include <iostream>

#include "wrap_orbit_mpi.hh"

//wrappers of mpi objects
#include "wrap_mpi_comm.hh"
#include "wrap_mpi_group.hh"
#include "wrap_mpi_status.hh"
#include "wrap_mpi_request.hh"
#include "wrap_mpi_datatype.hh"
#include "wrap_mpi_op.hh"


namespace wrap_orbit_mpi{
	
  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }
	
	#ifdef __cplusplus
	extern "C" {
		#endif
		
    static PyObject* mpi_get_processor_name(PyObject *self, PyObject *args) {
      char* name = new char[MPI_MAX_PROCESSOR_NAME];
      int len;
      ORBIT_MPI_Get_processor_name(name,&len);
      //PyObject* nm = PyString_FromStringAndSize(name, len);
			PyObject* nm = Py_BuildValue("s#",name, len);
      delete [] name;
      return nm;
    }
		
    static PyObject* mpi_comm_size(PyObject *self, PyObject *args) {
			if(PyTuple_Size(args) != 1){
				error("MPI_Comm_size(MPI_Comm) needs MPI_Comm as a parameter.");
			}			
			pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
      int result = 0;
      ORBIT_MPI_Comm_size(pyComm->comm,&result);
			return Py_BuildValue("i",result);
    }
		
    static PyObject* mpi_comm_rank(PyObject *self, PyObject *args) {
			if(PyTuple_Size(args) != 1){
				error("MPI_Comm_rank(MPI_Comm) needs MPI_Comm as a parameter.");
			}			
			pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
      int result = 0;
      ORBIT_MPI_Comm_rank(pyComm->comm,&result);
			return Py_BuildValue("i",result);
    }
		
    static PyObject* mpi_wtime(PyObject *self, PyObject *args) {
      double result ;
      result = (double) ORBIT_MPI_Wtime();
      return Py_BuildValue("d",result);
    }
		
    static PyObject* mpi_comm_group(PyObject *self, PyObject *args) {
			if(PyTuple_Size(args) != 2){
				error("MPI_Comm_group(MPI_Comm, MPI_Group) needs two parameters.");
			}
			pyORBIT_MPI_Comm* pyComm = (pyORBIT_MPI_Comm*) PyTuple_GetItem(args,0);
			pyORBIT_MPI_Group* pyGroup = (pyORBIT_MPI_Group*) PyTuple_GetItem(args,1);
			ORBIT_MPI_Comm_group(pyComm->comm,&pyGroup->group);
			return Py_BuildValue("i",1);
		}
		
		//Finalizes the execution of program
		//  the action is depended on the number of arguments
		//  () - no message
		//  (message) - will print message
		//this is a wrapper of
		// ORBIT_MPI_Finalize(const char* message)
		static PyObject* finalize(PyObject *self, PyObject *args){
			//if nVars == 0 no message
			//if nVars == 1 stop with message
			int nVars = PyTuple_Size(args);
			
			const char* message = NULL;
			
			if(nVars == 0 ||  nVars == 1){
				
				if(nVars == 1){
					//NO NEW OBJECT CREATED BY PyArg_ParseTuple!
					//NO NEED OF Py_DECREF()
					if(!PyArg_ParseTuple(	args,"s:finalize",&message)){
						ORBIT_MPI_Finalize("orbit_mpi - something wrong with error message.");
					}
					
					ORBIT_MPI_Finalize(message);
				}
				
			}
			else{
				ORBIT_MPI_Finalize("orbit_mpi. You should call finalize() or finalize(message)");
			}
			
			Py_INCREF(Py_None);
			return Py_None;
		}
		
		static PyMethodDef orbit_mpiMethods[] = {
			{ (char *)"MPI_Get_processor_name",  mpi_get_processor_name,   METH_VARARGS },
			{ (char *)"MPI_Comm_size",           mpi_comm_size,            METH_VARARGS },
			{ (char *)"MPI_Comm_rank",           mpi_comm_rank,            METH_VARARGS },
			{ (char *)"MPI_Wtime",               mpi_wtime,                METH_VARARGS },
			{ (char *)"MPI_Comm_group",          mpi_comm_group,           METH_VARARGS },
			{ (char *)"finalize",                finalize,                 METH_VARARGS },
			{ NULL, NULL }
		};
		
		
		void initorbit_mpi(void) {
			PyObject *m, *d;
			m = Py_InitModule((char*)"orbit_mpi",orbit_mpiMethods);
			d = PyModule_GetDict(m);
			
			//add MPI_Comm class and fields
			wrap_orbit_mpi_comm::init_orbit_mpi_comm(m);
			
			//add MPI_Group class and fields
			wrap_orbit_mpi_group::init_orbit_mpi_group(m);
			
			//add MPI_Status class and fields
			wrap_orbit_mpi_status::init_orbit_mpi_status(m);
			
			//add MPI_Request class and fields
			wrap_orbit_mpi_request::init_orbit_mpi_request(m);
			
			//add MPI_Datatype class and fields
			wrap_orbit_mpi_datatype::init_orbit_mpi_datatype(m);
			
			//add MPI_Op class and fields
			wrap_orbit_mpi_op::init_orbit_mpi_op(m);
		}
		
		
		#ifdef __cplusplus
	}
	#endif
	
	
	//end of namespace
}
