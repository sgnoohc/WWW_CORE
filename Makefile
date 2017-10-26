include Makefile.arch

#
# stuff to make
#
SOURCES=$(wildcard *.cc)
OBJECTS=$(SOURCES:.cc=.o)
LIB=WWW_CORE.so

#
# how to make it 
#

$(LIB): $(OBJECTS) 
	$(LD) $(LDFLAGS) $(SOFLAGS) $(OBJECTS) $(ROOTLIBS) -lTMVA -lEG -lGenVector -lXMLIO -lMLP -lTreePlayer -o $@

%.o:	%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

#
# target to build
# likelihood id library
#

rename:
	$(CXX) $(CXXFLAGS) RenameMapping.C -o $@

all: $(LIB) 
clean:
	rm -f *.o \
	rm -f *.d \
	rm -f *.so \
