include ${PETSC_DIR}/lib/petsc/conf/variables
LAPACKE_LIB=-llapacke -llapack -lblas
LIBS = ${PETSC_KSP_LIB} ${LAPACKE_LIB}
PETSC_INC=-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
INC = ${PETSC_INC}

all: AD2D

clean:
	rm advdiff.o VMS.o IGA.o utilities.o AD2D.o AD2D
		
AD2D: AD2D.o utilities.o IGA.o VMS.o advdiff.o
	g++ -g -o AD2D AD2D.o utilities.o IGA.o VMS.o advdiff.o ${LIBS} ${INC}
		
AD2D.o: AD2D.cpp utilities.h IGA.h VMS.h advdiff.h
	g++ -c -g AD2D.cpp ${LIBS} ${INC}
	
utilities.o: utilities.cpp utilities.h IGA.h
	g++ -c -g utilities.cpp ${LIBS} ${INC}
	
IGA.o: IGA.cpp IGA.h
	g++ -c -g IGA.cpp ${LIBS} ${INC}
	
VMS.o: VMS.cpp VMS.h IGA.h
	g++ -c -g VMS.cpp ${LIBS} ${INC}
	
advdiff.o: advdiff.cpp advdiff.h IGA.h VMS.h
	g++ -c -g advdiff.cpp ${LIBS} ${INC}