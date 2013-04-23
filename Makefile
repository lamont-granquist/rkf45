all:
	nvcc cuda/Rkf45.cu cuda/Launcher.cu cuda/clireSimulatedIRPaths.cu cuda/clireDeviceLibrary.cu cuda/clireDevice.cu
