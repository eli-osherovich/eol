FC = gfortran
FFLAGS = -g -shared -fPIC
LDFLAGS = -lblas

slibs := $(patsubst %.f,%.so,$(wildcard *.f))


all: $(slibs)

%.so: %.f
	$(FC) $(FFLAGS) $< $(LDFLAGS) -o $@



.PHONY : clean
clean :
	-rm *.so
