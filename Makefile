CPPFLAGS=-std=c++11 -Wall -O2 -l pthread

all: maxcut revmax preproc ba

maxcut: src/main.cpp src/mygraph.cpp src/algs.cpp
	g++ src/main.cpp -o maxcut  ${CPPFLAGS}
revmax: src/main.cpp src/mygraph.cpp src/algs.cpp
	g++ src/main.cpp -o revmax  ${CPPFLAGS} -DREVMAX
debug: src/main.cpp src/mygraph.cpp src/algs.cpp
	g++ src/main.cpp -o maxcut_debug -std=c++11 -Wall -Og -g 
preproc: src/preprocess.cpp src/mygraph.cpp
	g++ src/preprocess.cpp -o preproc  ${CPPFLAGS}
er: src/gen_er.cpp src/mygraph.cpp
	g++ -std=c++11 src/gen_er.cpp -o er
ba: src/gen_ba.cpp src/mygraph.cpp
	g++ -std=c++11 src/gen_ba.cpp -o ba
clean:
	rm ba maxcut preproc revmax 
