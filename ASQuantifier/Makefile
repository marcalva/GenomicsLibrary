CC = g++
FLAGS = -W -Wall -std=c++11 -stdlib=libc++
LD_LIBS = -L src/ -L /usr/local/lib/ -I include -I /usr/local/include -lbamtools -lz

ASQuantifier: src/main.cpp
	$(CC) $(FLAGS) src/*cpp -o $@ $(LD_LIBS) 
