FC = gfortran
FFLAGS = -O3
LIBS = -llapack -lblas

TARGET = H2_podvr
SRC = podvr.f90 main.f90

$(TARGET): $(SRC)
		$(FC) $(FFLAGS) $(SRC) $(LIBS) -o $(TARGET)

clean:
		rm -f $(TARGET) *.o *.mod