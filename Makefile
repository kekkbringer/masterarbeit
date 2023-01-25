IDIR = ../include
CC = g++
CFLAGS = -I$(IDIR) -std=c++20 -O3 -flto=auto -static #-mavx2 -mfma #-ffast-math #-march=native #-fopenmp #-Wpedantic
LFLAGS = -O3 -flto=auto -static #-mavx2 -mfma #-ffast-math #-march=nativ #-fopenmp
SANITIZE = #-fsanitize=address,undefined

ODIR = obj

all: ../bin/polly
rhf: ../bin/rhf

debug: CFLAGS += -gdwarf-4 -gstrict-dwarf
debug: LFLAGS += -gdwarf-4 -gstrict-dwarf
debug: ../bin/polly

../bin/polly: $(ODIR)/main.o $(ODIR)/cphf.o $(ODIR)/read_fourcenter.o $(ODIR)/read_herm.o $(ODIR)/read_spinor.o $(ODIR)/calc_stabmat.o $(ODIR)/fci_grad.o $(ODIR)/misc.o $(ODIR)/berry_rhs.o
	$(CC) -o ../bin/polly $(ODIR)/main.o $(ODIR)/cphf.o $(ODIR)/read_fourcenter.o $(ODIR)/read_herm.o $(ODIR)/read_spinor.o $(ODIR)/calc_stabmat.o $(ODIR)/fci_grad.o $(ODIR)/misc.o $(ODIR)/berry_rhs.o $(LFLAGS)

../bin/rhf: $(ODIR)/rhf.o $(ODIR)/cphf.o $(ODIR)/read_fourcenter.o $(ODIR)/read_herm.o $(ODIR)/read_spinor.o $(ODIR)/calc_stabmat.o $(ODIR)/fci_grad.o $(ODIR)/misc.o $(ODIR)/berry_rhs.o
	$(CC) -o ../bin/rhf $(ODIR)/rhf.o $(ODIR)/cphf.o $(ODIR)/read_fourcenter.o $(ODIR)/read_herm.o $(ODIR)/read_spinor.o $(ODIR)/calc_stabmat.o $(ODIR)/fci_grad.o $(ODIR)/misc.o $(ODIR)/berry_rhs.o $(LFLAGS)

$(ODIR)/calc_stabmat.o: calc_stabmat.cpp calc_stabmat.hpp
	$(CC) -c calc_stabmat.cpp -o $(ODIR)/calc_stabmat.o $(CFLAGS) $(SANITIZE)

$(ODIR)/misc.o: misc.cpp misc.hpp
	$(CC) -c misc.cpp -o $(ODIR)/misc.o $(CFLAGS) $(SANITIZE)

$(ODIR)/berry_rhs.o: berry_rhs.cpp berry_rhs.hpp
	$(CC) -c berry_rhs.cpp -o $(ODIR)/berry_rhs.o $(CFLAGS) $(SANITIZE)

$(ODIR)/cphf.o: cphf.cpp cphf.hpp
	$(CC) -c cphf.cpp -o $(ODIR)/cphf.o $(CFLAGS) $(SANITIZE)

$(ODIR)/read_fourcenter.o: read_fourcenter.cpp read_fourcenter.hpp
	$(CC) -c read_fourcenter.cpp -o $(ODIR)/read_fourcenter.o $(CFLAGS) $(SANITIZE)

$(ODIR)/read_herm.o: read_herm.cpp read_herm.hpp
	$(CC) -c read_herm.cpp -o $(ODIR)/read_herm.o $(CFLAGS) $(SANITIZE)

$(ODIR)/read_spinor.o: read_spinor.cpp read_spinor.hpp
	$(CC) -c read_spinor.cpp -o $(ODIR)/read_spinor.o $(CFLAGS) $(SANITIZE)

$(ODIR)/fci_grad.o: fci_grad.cpp fci_grad.hpp
	$(CC) -c fci_grad.cpp -o $(ODIR)/fci_grad.o $(CFLAGS) $(SANITIZE)

$(ODIR)/main.o: main.cpp $(ODIR)
	$(CC) -c main.cpp -o $(ODIR)/main.o $(CFLAGS) $(SANITIZE)

$(ODIR)/rhf.o: rhf.cpp $(ODIR)
	$(CC) -c rhf.cpp -o $(ODIR)/rhf.o $(CFLAGS) $(SANITIZE)

$(ODIR):
	[ -d $(ODIR) ] || mkdir -p $(ODIR)

.PHONY: static clean

static: $(ODIR)/main.o $(ODIR)/cphf.o $(ODIR)/read_fourcenter.o $(ODIR)/read_herm.o $(ODIR)/read_spinor.o $(ODIR)/calc_stabmat.o $(ODIR)/fci_grad.o $(ODIR)/misc.o $(ODIR)/berry_rhs.o
	$(CC) -o ../bin/polly_sl $(ODIR)/main.o $(ODIR)/cphf.o $(ODIR)/read_fourcenter.o $(ODIR)/read_herm.o $(ODIR)/read_spinor.o $(ODIR)/calc_stabmat.o $(ODIR)/fci_grad.o $(ODIR)/misc.o $(ODIR)/berry_rhs.o -static $(LFLAGS)


clean:
	rm -f ../bin/polly ../bin/polly_sl $(ODIR)/*.o
