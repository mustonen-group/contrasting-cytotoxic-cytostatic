objt = traits.cpp population.cpp patch.cpp simulation.cpp main.cpp
cc = g++
cflags = -Wall -Wextra -Wno-unused-but-set-parameter -std=c++17 -march=native -fopenmp -O2
cflags_debug = -Wall -Wextra -std=c++17 -march=native -g
cflags_callgrind = -Wall -Wextra -Wno-unused-but-set-parameter -std=c++17 -march=native -g -O2
lflags = -ljsoncpp -I/home/jva/include
name = sim
name_debug = sim_dbg
name_callgrind = sim_cg

all : $(objt)
	$(cc) $(objt) $(cflags) $(lflags) -o $(name)

debug : $(objt)
	$(cc) $(objt) $(cflags_debug) $(lflags) -o $(name_debug)

callgrind : $(objt)
	$(cc) $(objt) $(cflags_callgrind) $(lflags) -o $(name_callgrind)
