FC = gfortran
FFLAGS = -g -shared -fPIC -finit-integer=1
LDFLAGS = -lblas

slibs := $(patsubst %.f,%.so,$(wildcard *.f))


all: $(slibs)

%.so: %.f
	$(FC) $(FFLAGS) $< $(LDFLAGS) -o $@



.PHONY : clean
clean :
	-rm *.so
