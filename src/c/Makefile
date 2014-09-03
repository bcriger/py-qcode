NAME = qcode_dist

#Sources and objects for making the library/test program:
LIB_SRCS = ${NAME}.c
LIB_OBJS = ${NAME}.o
TEST_SRCS = ${NAME}_test.c
TEST_OBJS = ${NAME}_test.o

#Compiler and flags
CC = gcc 
#CFLAGS := -g -fPIC -D_DEBUG -D_DEBUG_NEIGHBOURS #For debugging
CFLAGS := -O3 -fPIC #For optimization

#Source to include when compiling
INCLUDES :=

#Directory to search for shared libraries when linking/loading
LIBDIR := -L$(shell pwd) -Wl,-rpath=$(shell pwd)

#Shared Libraries
LIBS = -l${NAME} -lrt

.cpp.o:
	$(CC) $(CFLAGS) ${INCLUDES} $< -c -o $@

shared_library: ${LIB_SRCS}
	$(CC) $(CFLAGS) -shared -o lib${NAME}.so ${LIB_SRCS}

test_program: ${TEST_OBJS}
	$(CC) ${LIBDIR} -o $@ ${TEST_OBJS} ${LIBS}

clean:
	rm *.o *.so test_program

all: shared_library test_program