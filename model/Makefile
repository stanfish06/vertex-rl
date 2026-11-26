CXX ?= g++
CXXFLAGS ?= -O3
OPENMP ?= -fopenmp
JSON_LIBS ?= -ljson-c
READLINE_LIBS ?= -lreadline
GEOMETRY = three_d/mesh.cpp three_d/contact.cpp three_d/solver.cpp
DRIVER = $(GEOMETRY) three_d/config.cpp three_d/export.cpp three_d/driver.cpp
HEADERS3D = $(wildcard three_d/*.h)
HEADERS2D = $(wildcard *.h)

.PHONY: all
all: vertex vertex3d libvertex2d.so

vertex3d: three_d/main.cpp $(DRIVER) $(HEADERS3D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(OPENMP) three_d/main.cpp $(DRIVER) $(LDFLAGS) $(JSON_LIBS) -o $@

vertex: main.c $(DRIVER) $(HEADERS2D) $(HEADERS3D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(OPENMP) -DVERTEX_WITH_3D main.c $(DRIVER) $(LDFLAGS) $(READLINE_LIBS) $(JSON_LIBS) -o $@

libvertex2d.so: worker_api.cpp $(HEADERS2D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(OPENMP) -fPIC -shared worker_api.cpp $(LDFLAGS) -o $@
