CC = g++
CCO = $(CC) -c

all : bin/price

obj/fdm.o : src/fdm.cpp include/fdm.h include/pde.h
	$(CCO) $< -o $@

obj/option.o : src/option.cpp include/option.h 
	$(CCO) $< -o $@

obj/pde.o : src/pde.cpp include/pde.h include/option.h
	$(CCO) $< -o $@

obj/pricers.o : src/pricers.cpp include/pricers.h
	$(CCO) $< -o $@

obj/greeks.o : src/greeks.cpp include/greeks.h
	$(CCO) $< -o $@

obj/graph_builder.o : src/graph_builder.cpp include/graph_builder.h
	$(CCO) $< -o $@

obj/vba_interface.o : src/vba_interface.cpp include/vba_interface.h
	$(CCO) $< -o $@

bin/price : obj/pde.o obj/option.o obj/fdm.o obj/pricers.o obj/greeks.o obj/graph_builder.o obj/vba_interface.o src/main.cpp 
	$(CC) $^ -o bin/price
