
all: mandel series

mandel: mandel.o bitmap.o
	gcc mandel.o bitmap.o -o mandel -lpthread

series: series.o bitmap.o
	gcc series.o bitmap.o -o series -lpthread

mandel.o: mandel.c
	gcc -Wall -g -c mandel.c -o mandel.o

series.o: series.c
	gcc -Wall -g -c series.c -o series.o

bitmap.o: bitmap.c
	gcc -Wall -g -c bitmap.c -o bitmap.o

clean:
	rm -f mandel.o series.o bitmap.o mandel series

 
