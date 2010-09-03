CC=gcc
CFLAGS=-O2 -std=gnu99 -Wall -g -ggdb

SRC_DIR=src
BIN_DIR=bin

all:	clean

tests:	test_square_matrix.o test_magnetic_monopole.o

square_matrix.o:
	${CC} ${CFLAGS} -c -fopenmp \
	-o ${BIN_DIR}/square_matrix.o \
	${SRC_DIR}/square_matrix.c

test_square_matrix.o: square_matrix.o
	${CC} ${CFLAGS} -fopenmp \
	-lm \
	-o ${BIN_DIR}/test_square_matrix.o \
	${BIN_DIR}/square_matrix.o \
	${SRC_DIR}/test_square_matrix.c

magnetic_monopole.o: square_matrix.o vector.o
	${CC} ${CFLAGS} -c -fopenmp \
	-o ${BIN_DIR}/magnetic_monopole.o \
	${SRC_DIR}/magnetic_monopole.c

test_magnetic_monopole.o: magnetic_monopole.o square_matrix.o vector.o
	${CC} ${CFLAGS} -fopenmp \
	-lm \
	-o ${BIN_DIR}/test_magnetic_monopole.o \
	${BIN_DIR}/magnetic_monopole.o \
	${BIN_DIR}/square_matrix.o \
	${BIN_DIR}/vector.o \
	${SRC_DIR}/test_magnetic_monopole.c

vector.o: 
	${CC} ${CFLAGS} -c -fopenmp \
	-o ${BIN_DIR}/vector.o \
	${SRC_DIR}/vector.c

clean:
	rm -f ${BIN_DIR}/*

run:	qrfact
	${BIN_DIR}/qrfact

.PHONY: clean
