include ${PETSC_DIR}/lib/petsc/conf/variables
LAPACKE_LIB=-llapacke -llapack -lblas
LIBS = ${PETSC_KSP_LIB} ${LAPACKE_LIB}
PETSC_INC=-I${PETSC_DIR}include -I${PETSC_DIR}${PETSC_ARCH}/include
INC = ${PETSC_INC}

all: AD2D.exe

clean:
		rm AD2D.o AD2D.exe
		
AD2D.exe: AD2D.o
		g++ -g -o AD2D.exe AD2D.o ${LIBS} ${INC}
		
AD2D.o: AD2D.cpp
		g++ -c -g AD2D.cpp -lpetscksp ${LIBS} ${INC}