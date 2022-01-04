CC = g++
CCO = $(CC) -c

all : bin/price

#obj/fdm_static.o :

obj/fdm.o : src/fdm.cpp include/fdm.h include/pde.h src/fdm_static.cpp
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

bin/test_fdm : src/test_fdm.cpp src/fdm_static.cpp obj/pde.o obj/option.o
	$(CC) $< obj/pde.o obj/option.o -o bin/test_fdm -lcriterion

bin/test_pricers : src/test_pricers.cpp obj/pde.o obj/option.o obj/fdm.o obj/pricers.o
	$(CC) $^ -o bin/test_pricers -lcriterion

tests : bin/test_fdm bin/test_pricers
