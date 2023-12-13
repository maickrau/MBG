PLATFORM=$(shell uname -s)

GPP=$(CXX)

#  FreeBSD absolutely must have __XSI_VISIBLE=1 and _POSIX_C_SOURCE=200112L defined.
ifeq ($(PLATFORM), FreeBSD)
CPPFLAGS  = -D__XSI_VISIBLE=1
CPPFLAGS += -D_POSIX_C_SOURCE=200112L
endif

CPPFLAGS += -Wall -Wextra
CPPFLAGS += -Wno-unused-parameter
CPPFLAGS += -Wno-deprecated-declarations
CPPFLAGS += -Wno-unused-lambda-capture
CPPFLAGS += -Wno-unused-but-set-variable
CPPFLAGS += -std=c++17 -pthread
CPPFLAGS += -O3 -g
CPPFLAGS += -Izstr/src
CPPFLAGS += -Iparallel-hashmap/parallel_hashmap
CPPFLAGS += -Icxxopts/include
CPPFLAGS += -Iconcurrentqueue

CPPFLAGS += $(shell pkg-config --cflags zlib)
LIBS     += $(shell pkg-config --libs zlib)

ODIR=obj
BINDIR=bin
SRCDIR=src
LIBDIR=lib

_DEPS = fastqloader.h CommonUtils.h MBGCommon.h VectorWithDirection.h FastHasher.h SparseEdgeContainer.h HashList.h UnitigGraph.h BluntGraph.h ReadHelper.h HPCConsensus.h ErrorMaskHelper.h CompressedSequence.h ConsensusMaker.h StringIndex.h LittleBigVector.h MostlySparse2DHashmap.h RankBitvector.h TwobitLittleBigVector.h UnitigResolver.h CumulativeVector.h UnitigHelper.h BigVectorSet.h Serializer.h DumbSelect.h MsatValueVector.h Node.h KmerMatcher.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = MBG.o fastqloader.o CommonUtils.o MBGCommon.o FastHasher.o SparseEdgeContainer.o HashList.o UnitigGraph.o BluntGraph.o HPCConsensus.o ErrorMaskHelper.o CompressedSequence.o ConsensusMaker.o StringIndex.o RankBitvector.o UnitigResolver.o UnitigHelper.o BigVectorSet.o ReadHelper.o Serializer.o DumbSelect.o MsatValueVector.o Node.o KmerMatcher.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

#  MacOS isn't happy with static/dynamic flags.
ifeq ($(PLATFORM),Darwin)
   LINKFLAGS = $(LIBS) -pthread -static-libstdc++
else
   LINKFLAGS = -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -pthread -static-libstdc++
endif

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)
$(shell mkdir -p lib)

lib: $(LIBDIR)/mbg.a

$(LIBDIR)/mbg.a: $(OBJ) $(DEPS)
	ar rvs $@ $(OBJ)

$(BINDIR)/MBG: $(OBJ) $(ODIR)/main.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/main.o: $(SRCDIR)/main.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

all: $(BINDIR)/MBG

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
	rm -f $(LIBDIR)/*
