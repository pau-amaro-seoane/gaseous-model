INC1 = compar.f
INC2 = equiv.f
INC3 = params.f
#FFLAGS = -g -Mbounds -Mextend -Ktrap=fp
#FFLAGS = -fno-silent -ffixed-line-length-none -Wall -Wno-unused -march=i686 -g -fbounds-check  # g77 omptions
#FFLAGS = -fast -cg92
#FFLAGS = -fast -tp athlon -pc 64 -Mbounds
#FFLAGS = -Mextend -Bstatic -fast -fastsse -tp p6 -pc 64
#FFLAGS = -Bstatic -fast -fastsse -tp p7 -pc 64 # -Mbounds
#FFLAGS = -O6
#FFLAGS = -O2 -xtarget=ultra
#FFLAGS = -p -fast -cg92  # use prof -V../star afterwards
#FFLAGS = -g -C # -fast -xcg92
#FC = pgf77
FC = gfortran
SOURCE = star.f stevol.f dawson.f \
avint.f choose.f dens.f slope.f readpm.f readmo.f rgrid.f \
prints.f printo.f print.f i00eq.f i10eq.f init.f initst.f henyey.f \
i02eq.f i2ieq.f inform.f gidini.f rxterm.f i20eq.f inner.f outer.f \
ginit.f imreq.f girl.f info.f i30eq1.f i30eq.f i12eq.f formb2.f formb3.f \
formbs.f ediff.f lossc.f ahole.f lossco.f erf.f \
binsto.f relax.f escape.f orbit.f \
trmain.f derqp3.f stablz.f triple.f \
difsy3.f qpmod3.f tperi.f efac2.f efac3.f \
erel3.f stabl3.f trans3.f peri.f tides.f \
qumain.f quad.f newreg.f rchain.f endreg.f status.f \
trans4.f erel4.f derqp4.f difsy4.f rsort.f ichain.f

OBJECTS = $(SOURCE:.f=.o)

star:	$(OBJECTS) ranf.o
	$(FC) $(FFLAGS) $(OBJECTS) ranf.o -o ./spedi

starbig:	$(OBJECTS) ranf.o
		$(FC) $(FFLAGS) $(OBJECTS) ranf.o -o starbig

startest:	$(OBJECTS) ranf.o
		$(FC) $(FFLAGS) $(OBJECTS) ranf.o -o startest

starcheck:	$(OBJECTS) ranf.o
		$(FC) $(FFLAGS) $(OBJECTS) ranf.o -o starcheck

startrip:	$(OBJECTS) ranf.o
		$(FC) $(FFLAGS) $(OBJECTS) ranf.o -o startrip

$(OBJECTS): $(INC1) $(INC3)

print.o: $(INC2) $(INC3)

readmo.o: $(INC2) $(INC3)

trmain.o: $(INC3)

ranf.o:	ranf.f
	$(FC) -c ranf.f

clean:
	rm *.o
	rm spedi
