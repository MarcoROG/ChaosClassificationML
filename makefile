HEADERS = Headers/problem.h Headers/ProblemLibrary.h Headers/ExplicitEuler.h Headers/Parameters.h Headers/IO.h


genorb: GenOrbits.cxx
	g++ ${INCLUDE} GenOrbits.cxx ${HEADERS} ${MIZINT} -lm -o GenOrbits.o


genorb2: GenOrbits2.cxx
	g++ GenOrbits2.cxx ${HEADERS} -O3 -lm -o GenOrbits2.o

harmonic: HarmonicOsc.cxx
	g++ ${INCLUDE} HarmonicOsc.cxx ${HEADERS} ${MIZINT} -lm -o HarmonicOsc.o


harmoniceoc: HarmonicOscillator.cxx
	g++ ${INCLUDE} HarmonicOscillator.cxx ${HEADERS} ${MIZINT} -lm -o HarmonicEoc.o

all: harmonic henon

clearall:
	rm GenOrbits.o
