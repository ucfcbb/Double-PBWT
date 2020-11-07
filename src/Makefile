# COMPILER
CXX = g++

# COMPILER FLAGS
CXXFLAGS = -std=c++17 -Wall -pedantic -O3

# path to header files
INCLUDES = -I../include

#LFLAGS = -L/home/pramesh/lib

SRCS = main.cpp doublePBWT.cpp helperFunctions.cpp

OBJS = $(SRCS:.cpp=.o)

# executable file
MAIN = d2pbwt_overlap_genetic

#### TARGETS
all: $(MAIN)

$(MAIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^
# DO NOT DELETE

main.o: ../include/doublePBWT.h ../include/helperFunctions.h
doublePBWT.o: ../include/doublePBWT.h ../include/helperFunctions.h
helperFunctions.o: ../include/helperFunctions.h
