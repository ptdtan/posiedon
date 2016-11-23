all: ./src/rdis
	gcc -I./miniilib -lz -o a.out ./src/rdis.o ./minilib/bam.o

./src/rdis.o: ./src/rdis.c
	gcc -I./minilib -lz -o ./src/rdis.o -c ./src/rdis.c 

