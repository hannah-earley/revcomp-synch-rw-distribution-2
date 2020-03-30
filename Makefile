all: expt

include ~/.files/util/proj-term.mk
include ~/.files/util/proj-view.mk
include ~/.files/util/proj-text.mk

CXX=clang++
CXXFLAGS=-Xclang -fopenmp -std=c++11 -O3 -Wall -Iinc/
LINKFLAGS=-lomp

SRCDIR=src
OBJDIR=obj

_DEPS=
DEPS=$(patsubst %,$(SRCDIR)/%,$(_DEPS))
_OBJ=expt.o
OBJ=$(patsubst %,$(OBJDIR)/%,$(_OBJ))

expt: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LINKFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

run: expt
	./expt

clean:
	-rm -f $(OBJDIR)/*.o
	-rmdir $(OBJDIR)

.PHONY: all run clean

