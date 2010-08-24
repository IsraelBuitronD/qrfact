CC=gcc
CFLAGS=-O2 -std=gnu99 -Wall -g -ggdb

SRC_DIR=src
BIN_DIR=bin

all:	clean qrfact

tests:	test_matrix

qrfact.o:	
	${CC} ${CFLAGS} -c -fopenmp -o ${BIN_DIR}/qrfact.o ${SRC_DIR}/qrfact.c

matrix.o:
	${CC} ${CFLAGS} -c -o ${BIN_DIR}/matrix.o ${SRC_DIR}/matrix.c

test_matrix: matrix.o
	${CC} ${CFLAGS} -lm -o ${BIN_DIR}/test_matrix ${BIN_DIR}/matrix.o ${SRC_DIR}/test_matrix.c

qrfact:	matrix.o qrfact.o
	${CC} ${CFLAGS} -lm -fopenmp -o ${BIN_DIR}/qrfact ${BIN_DIR}/qrfact.o ${BIN_DIR}/matrix.o 

clean:
	rm -f ${BIN_DIR}/*

run:	qrfact
	${BIN_DIR}/qrfact

.PHONY: clean
