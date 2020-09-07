# NOTE make needs tabs in indented lines
# Need to install g++, make ccpunit-devel
appname := ocusf_tests.exe

CXX := g++
CXXFLAGS := 
LDFLAGS := 
LDLIBS := -lcppunit
srcfiles := main.cpp ../utils.cpp
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

debug: CXXFLAGS += -DDEBUG -g
debug: CCFLAGS += -DDEBUG -g
debug: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS)

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(objects)

dist-clean: clean
	rm -f *~ .depend

include .depend
