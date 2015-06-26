all: IGA2D.exe

clean:
		rm IGA2D.o IGA2D.exe
		
IGA2D.exe: IGA2D.o
		g++ -g -o IGA2D.exe IGA2D.o
		
IGA2D.o: IGA2D.cpp
		g++ -c -g IGA2D.cpp