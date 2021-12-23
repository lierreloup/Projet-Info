CC = g++
CCO = $(CC) -c

all : bin/price

obj/fdm.o : src/fdm.cpp include/fdm.h include/pde.h
	$(CCO) $< -o $@

obj/option.o : src/option.cpp include/option.h 
	$(CCO) $< -o $@

obj/pde.o : src/pde.cpp include/pde.h include/option.h
	$(CCO) $< -o $@

obj/pricers.o : src/pricers.cpp
	$(CCO) $< -o $@

bin/price : obj/pde.o obj/option.o obj/fdm.o obj/pricers.o src/main.cpp 
	$(CC) $^ -o bin/price