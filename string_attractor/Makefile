SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -O3 -DNDEBUG -march=native -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -std=c++0x

all: text_to_st_att test

text_to_st_att:
	$(CC) $(CFLAGS) -o text_to_st_att ./src/main.cpp ./src/utils.cpp

test_functionality:
	$(CC) $(CFLAGS) -o test_functionality ./test/main.cpp ./src/utils.cpp
clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f text_to_st_att *.o
