# NOTE make needs tabs in indented lines
appname := ocusf.exe

CXX := g++
CXXFLAGS := 
LDFLAGS := 
LDLIBS := 
#Need to explictly name files because don't want to include test files
#srcfiles := $(shell find . -name "*.cpp")
srcfiles := OCUSF.cpp utils.cpp
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
