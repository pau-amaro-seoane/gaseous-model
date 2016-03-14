#FFLAGS = -O   -tp p6 -pc 64
#FC = pgf77
#FFLAGS = -O6
FC = g77
SOURCE = getdata.f \
excomp.f pstart.f decode.f xplot.f dcode.f ncode.f 

OBJECTS = $(SOURCE:.f=.o)

star:	$(OBJECTS)
	$(FC) -fbounds-check $(FFLAGS) $(OBJECTS) -o getdata
