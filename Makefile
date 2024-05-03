
all: mandelseries

mandelseries: mandelseries.o bitmap.o
	gcc -std=c99 mandelseries.o bitmap.o -o mandelseries -lpthread

mandelseries.o: mandelseries.c
	gcc -std=c99 -Wall -g -c mandelseries.c -o mandelseries.o

bitmap.o: bitmap.c
	gcc -std=c99 -Wall -g -c bitmap.c -o bitmap.o

clean:
	rm -f mandelseries.o bitmap.o mandelseries
