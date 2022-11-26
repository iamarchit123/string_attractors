SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -O3 -DNDEBUG -std=c++0x
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CFLAGS += -D LINUX
endif
ifeq ($(UNAME_S),Darwin)
    CFLAGS += -D OSX
endif
all: text_to_st_att test

text_to_st_att:
	$(CC) $(CFLAGS) -o text_to_st_att ./src/main.cpp ./src/utils.cpp

test_functionality:
	$(CC) $(CFLAGS) -o test_functionality ./test/main.cpp ./src/utils.cpp

time_to_construct:
	$(CC) $(CFLAGS) -o test_functionality ./test/construction.cpp ./src/utils.cpp

time_to_query:
	$(CC) $(CFLAGS) -o test_functionality ./test/index_query.cpp ./src/utils.cpp

test_english:
	$(CC) $(CFLAGS) -o test_english ./test/main_english.cpp ./src/utils.cpp
clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f text_to_st_att *.o
