CC=gcc
CFLAGS=-O2 -std=gnu99 -Wall -g -ggdb

SRC_DIR=src
BIN_DIR=bin

all:	clean qrfact

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

vector.o: 
	${CC} ${CFLAGS} -c -fopenmp \
	-o ${BIN_DIR}/vector.o \
	${SRC_DIR}/vector.c

qrfact.o:	
	${CC} ${CFLAGS} -c -fopenmp \
	-o ${BIN_DIR}/qrfact.o \
	${SRC_DIR}/qrfact.c

qrfact:	qrfact.o
	${CC} ${CFLAGS} ${BIN_DIR}/qrfact.o -lm -fopenmp -o ${BIN_DIR}/qrfact

clean:
	rm -f ${BIN_DIR}/*

run:	qrfact
	${BIN_DIR}/qrfact

.PHONY: clean
