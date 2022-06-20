# (c) J. Blazek, CFD Consulting & Analysis, support@cfd-ca.de

EXEC = ../run/Unstruct2D

CC = g++
LD = g++
OBJ = o
LIBS = -lm

CFLAGS = -O3 -m64 -Wall -std=c++11 -c
LDFLAGS = -O3 -m64 -s

SOURCES = bcondFarfield.cpp bcondInflow.cpp bcondOutflow.cpp bcondWallVisc.cpp \
          bndConds.cpp boundaryConditions.cpp checkMetrics.cpp computeMetrics.cpp \
          deleteTmpElist.cpp densityChange.cpp dependentVars.cpp dissipInit.cpp \
          dissipRoe1.cpp dissipRoe1Prec.cpp dissipRoe2.cpp dissipRoe2Prec.cpp \
          dummyNodes.cpp faceVectorsSymm.cpp faceVectorsVolumes.cpp \
          faceVectorsVolumesBound.cpp fluidProps.cpp fluxRoe1.cpp fluxRoe2.cpp \
          fluxViscous.cpp fluxWalls.cpp forces.cpp generateEdgelist.cpp \
          geometry.cpp getNumberBoundNodes.cpp gradients.cpp gradientsVisc.cpp \
          initConstants.cpp initSolution.cpp irsmoo.cpp limiter.cpp \
          limiterInit.cpp limiterRefVals.cpp massFlow.cpp outputConvergence.cpp \
          outputFlowfield.cpp outputSurfaces.cpp periodic.cpp readGrid.cpp \
          readSolution.cpp solve.cpp timeStep.cpp userInput.cpp main.cpp \
          output.cpp solver.cpp spaceDiscr.cpp timeDiscr.cpp \
          volumeProjections.cpp convergence.cpp writeSolution.cpp \
          zeroResiduals.cpp

OBJECTS = ${SOURCES:.cpp=.$(OBJ)}
.SUFFIXES: .cpp .$(OBJ)

.cpp.$(OBJ):
	$(CC) $(CFLAGS) $*.cpp

$(EXEC): $(OBJECTS)
	$(LD) $(LDFLAGS) -o $(EXEC) $(OBJECTS) $(LIBS)

clean:
	$(RM) *.$(OBJ)
