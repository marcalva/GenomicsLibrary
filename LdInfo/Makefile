CC = g++
CFLAGS = -g -Wall -Wextra -std=c++11
IFLAGS = -I tabixpp/ -I tabixpp/htslib/ -I ../ParseText/
LFLAGS = -L tabixpp/htslib -lhts -lz

LdInfo : main.cpp LdInfo.o VcfRecord.o tabixpp/tabix.o
	$(CC) $(CFLAGS) $(IFLAGS) $(LFLAGS) -o LdInfo main.cpp LdInfo.o VcfRecord.o tabixpp/tabix.o

tabixpp/tabix.o : tabixpp/tabix.hpp tabixpp/tabix.cpp
	cd tabixpp
	make
	cd ..

LdInfo.o :  LdInfo.cpp LdInfo.h VcfRecord.h ../ParseText/ParseText.h
	$(CC) $(CFLAGS) $(IFLAGS) $(LFLAGS) -c LdInfo.cpp VcfRecord.cpp

VcfRecord.o : VcfRecord.cpp VcfRecord.h ../ParseText/ParseText.h
	$(CC) $(CFLAGS) $(IFLAGS) $(LFLAGS) -c VcfRecord.cpp
