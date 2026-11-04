BUILT_OBJECTS += obj/swiftcov_fast.o

# SwiftCov fast path: compiled at -O3 with C++17; upstream objects keep the
# release flags passed on the make command line.
obj/swiftcov_fast.o: override CXXFLAGS := -O3 -g0 -std=c++17

obj/swiftcov_fast.o: src/swiftcov/swiftcov_fast.cpp obj/swiftcov_fast.d
	$(CXX_COMPILE)
