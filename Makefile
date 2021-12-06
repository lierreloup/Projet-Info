CC = g++
CCO = $(CC) -c

all : bin/price

obj/fdm.o : src/fdm.cpp include/fdm.h include/pde.h
	$(CCO) $< -o $@

obj/option.o : src/option.cpp include/option.h 
	$(CCO) $< -o $@

obj/payoff.o : src/payoff.cpp include/payoff.h 
	$(CCO) $< -o $@

obj/pde.o : src/pde.cpp include/pde.h include/option.h
	$(CCO) $< -o $@

bin/price : obj/pde.o obj/payoff.o obj/option.o obj/fdm.o src/main.cpp
	$(CC) $^ -o bin/price