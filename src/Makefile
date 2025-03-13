CXX           = icpc
CC            = icpc
CXXFLAGS      = -O3 \
                -DNDEBUG 
DEST          = ./
OBJS          = Ellipsoids.o \
                Cuboids.o \
                main.o \
                Lake.o \
                MeshGenerationDHexa.o \
                ObservationLine.o \
                ObservationPoint.o \
                ObservingSiteList.o \
                TopographyData.o \
                Util.o
PROGRAM       = makeDHexaMesh

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM)

clean:;		rm -f *.o *~ $(PROGRAM)
