GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Izstr/src -Iparallel-hashmap/parallel_hashmap/ -Wno-unused-parameter -Icxxopts/include -Iconcurrentqueue

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=-lz

_DEPS = fastqloader.h CommonUtils.h MBGCommon.h VectorWithDirection.h FastHasher.h SparseEdgeContainer.h HashList.h UnitigGraph.h BluntGraph.h ReadHelper.h HPCConsensus.h ErrorMaskHelper.h CompressedSequence.h ConsensusMaker.h StringIndex.h LittleBigVector.h MostlySparse2DHashmap.h RankBitvector.h TwobitLittleBigVector.h UnitigResolver.h CumulativeVector.h UnitigHelper.h BigVectorSet.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = MBG.o fastqloader.o CommonUtils.o main.o MBGCommon.o FastHasher.o SparseEdgeContainer.o HashList.o UnitigGraph.o BluntGraph.o HPCConsensus.o ErrorMaskHelper.o CompressedSequence.o ConsensusMaker.o StringIndex.o RankBitvector.o UnitigResolver.o UnitigHelper.o BigVectorSet.o ReadHelper.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(BINDIR)/MBG: $(OBJ)
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/main.o: $(SRCDIR)/main.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

all: $(BINDIR)/MBG

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
