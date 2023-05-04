MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Config
FC := ifort
LD := $(FC)
RM := rm -f

# Source files
SRCS := ConvergenceResults.f90 \
	ConvEquiResults.f90 \
	DetailedAnalysis.f90 \
	EquilibriumCheck.f90 \
	generic_routines.f90 \
	globals.f90 \
	ImpulseResponse.f90 \
	LearningSimulation.f90 \
	LearningTrajectory.f90 \
	main.f90 \
	PI_routines.f90 \
	QGapToMaximum.f90 \
	QL_routines.f90
PROG := baseline
OBJS := $(subst .f90,.o,$(SRCS))

# Build commands
.PHONY: all clean
all: $(PROG)

$(PROG): $(OBJS)
	$(LD) -o $@ $+ -qopenmp

%.o %.mod: %.f90
	$(FC) -c -o $*.o $< -qopenmp
	@touch $@

# $(OBJS): %.o: %
# 	$(FC) -c -o $@ $< -qopenmp

# define dependencies between object files
ConvergenceResults.o: globals.o globals.mod QL_routines.o QL_routines.mod
ConvEquiResults.o: globals.o globals.mod QL_routines.o QL_routines.mod \
	generic_routines.o generic_routines.mod EquilibriumCheck.o EquilibriumCheck.mod
DetailedAnalysis.o: globals.o globals.mod QL_routines.o QL_routines.mod \
	ImpulseResponse.o ImpulseResponse.mod QGapToMaximum.o QGapToMaximum.mod
EquilibriumCheck.o: globals.o globals.mod QL_routines.o QL_routines.mod \
	generic_routines.o generic_routines.mod
globals.o: generic_routines.o generic_routines.mod
ImpulseResponse.o: globals.o globals.mod QL_routines.o QL_routines.mod \
	EquilibriumCheck.o EquilibriumCheck.mod
LearningSimulation.o: globals.o globals.mod QL_routines.o QL_routines.mod
LearningTrajectory.o: globals.o globals.mod QL_routines.o QL_routines.mod \
	generic_routines.o generic_routines.mod LearningSimulation.o LearningSimulation.mod \
	ImpulseResponse.o ImpulseResponse.mod
PI_routines.o: globals.o globals.mod QL_routines.o QL_routines.mod \
	generic_routines.o generic_routines.mod 
QGapToMaximum.o: globals.o globals.mod QL_routines.o QL_routines.mod \
	generic_routines.o generic_routines.mod EquilibriumCheck.o EquilibriumCheck.mod
QL_routines.o: globals.o globals.mod

main.o: globals.o globals.mod \
	LearningSimulation.o LearningSimulation.mod \
	ConvergenceResults.o ConvergenceResults.mod \
	ImpulseResponse.o ImpulseResponse.mod \
	EquilibriumCheck.o EquilibriumCheck.mod \
	QGapToMaximum.o QGapToMaximum.mod \
	LearningTrajectory.o LearningTrajectory.mod \
	DetailedAnalysis.o DetailedAnalysis.mod \
	QL_routines.o QL_routines.mod \
	PI_routines.o PI_routines.mod \
	generic_routines.o generic_routines.mod \
	ConvEquiResults.o ConvEquiResults.mod 

# rebuild all object files in case this Makefile changes
$(OBJS): $(MAKEFILE_LIST)

clean:
	$(RM) $(filter %.o, $(OBJS)) $(wildcard *.mod) 
