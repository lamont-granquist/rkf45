all: Rkf45.o
	nvcc -arch=sm_20 -rdc true Launcher.cu Rkf45.o 

Rkf45.o:
	nvcc -arch=sm_20 -dc Rkf45.cu LinkTest.cu
