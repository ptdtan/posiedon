all: ./src/rdis.o ./minilib/bam.o ./minilib/bgzf.o
	gcc -I./minilib -L./zlib -o rdis ./src/rdis.o ./minilib/bam.o ./minilib/bgzf.o -lz

./minilib/bam.o: ./minilib/bam.c
	gcc -I./minilib -L./zlib -lz -o ./minilib/bam.o -c ./minilib/bam.c
./minilib/bgzf.o: ./minilib/bgzf.c
	gcc -I./minilib -L./zlib -lz -o ./minilib/bgzf.o -c ./minilib/bgzf.c
./src/rdis.o: ./src/rdis.c
	gcc -I./minilib -L./zlib -lz -o ./src/rdis.o -c ./src/rdis.c
