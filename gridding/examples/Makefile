wsgridderexample:	wspredictionexample.cpp ../wstackinggridder.cpp ../../fftwmultithreadenabler.cpp
	g++ -Wall -o wspredictionexample -I ../../external/aocommon/include/ -std=c++17 -DAVOID_CASACORE wspredictionexample.cpp ../wstackinggridder.cpp ../../fftwmultithreadenabler.cpp -lfftw3 -lfftw3_threads -lboost_date_time -lboost_thread -lboost_system

channelfitexample:	channelfitexample.cpp
	g++ -Wall -o channelfitexample -I ../../external/aocommon/include/ -std=c++17 channelfitexample.cpp ../../polynomialchannelfitter.cpp ../../polynomialfitter.cpp -lgsl -lgslcblas
