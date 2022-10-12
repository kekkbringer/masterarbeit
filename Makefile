IDIR = ../include
CC = g++
CFLAGS = -I$(IDIR) -std=c++20 -O3 -flto=auto -ffast-math #-march=native #-mavx2 -mfma #-fopenmp #-Wpedantic
LFLAGS = -O3 -flto=auto -ffast-math #-march=native #-mavx2 -mfma #-fopenmp
SANITIZE = #-fsanitize=address,undefined

ODIR = obj

all: ../bin/polly

../bin/polly: $(ODIR)/main.o $(ODIR)/cphf.o $(ODIR)/read_fourcenter.o $(ODIR)/read_herm.o $(ODIR)/read_spinor.o $(ODIR)/calc_stabmat.o
	$(CC) -o ../bin/polly $(ODIR)/main.o $(ODIR)/cphf.o $(ODIR)/read_fourcenter.o $(ODIR)/read_herm.o $(ODIR)/read_spinor.o $(ODIR)/calc_stabmat.o $(LFLAGS)

$(ODIR)/calc_stabmat.o: calc_stabmat.cpp calc_stabmat.hpp
	$(CC) -c calc_stabmat.cpp -o $(ODIR)/calc_stabmat.o $(CFLAGS) $(SANITIZE)

$(ODIR)/cphf.o: cphf.cpp cphf.hpp
	$(CC) -c cphf.cpp -o $(ODIR)/cphf.o $(CFLAGS) $(SANITIZE)

$(ODIR)/read_fourcenter.o: read_fourcenter.cpp read_fourcenter.hpp
	$(CC) -c read_fourcenter.cpp -o $(ODIR)/read_fourcenter.o $(CFLAGS) $(SANITIZE)

$(ODIR)/read_herm.o: read_herm.cpp read_herm.hpp
	$(CC) -c read_herm.cpp -o $(ODIR)/read_herm.o $(CFLAGS) $(SANITIZE)

$(ODIR)/read_spinor.o: read_spinor.cpp read_spinor.hpp
	$(CC) -c read_spinor.cpp -o $(ODIR)/read_spinor.o $(CFLAGS) $(SANITIZE)

$(ODIR)/main.o: main.cpp
	$(CC) -c main.cpp -o $(ODIR)/main.o $(CFLAGS) $(SANITIZE)

.PHONY: static clean

static: $(ODIR)/main.o $(ODIR)/cphf.o $(ODIR)/read_fourcenter.o $(ODIR)/read_herm.o $(ODIR)/read_spinor.o $(ODIR)/calc_stabmat.o
	$(CC) -o ../bin/polly_sl $(ODIR)/main.o $(ODIR)/cphf.o $(ODIR)/read_fourcenter.o $(ODIR)/read_herm.o $(ODIR)/read_spinor.o $(ODIR)/calc_stabmat.o $(LFLAGS) -static


clean:
	rm -f ../bin/polly ../bin/polly_sl $(ODIR)/*.o
