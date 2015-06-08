all:	ls

ls:	ls.o
	g++ -I/opt/ibm/cplex/include -I/opt/ibm/concert/include ls.o -L/opt/ibm/cplex/lib/x86-64_sles10_4.1/static_pic -L/opt/ibm/concert/lib/x86-64_sles10_4.1/static_pic  -O3 -mtune=native -fomit-frame-pointer -lpthread -lconcert -lilocplex -lcplex -lm -pthread  -lm -m64  -lz -lreadline -lncurses -o LS

ls.o:	ls.cpp
	g++ -I/opt/ibm/cplex/include -I/opt/ibm/concert/include  -O3 -Wall -m64 -std=c++0x -DIL_STD -c ls.cpp -o ls.o
	
