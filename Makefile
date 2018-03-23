
all: mandel series master_mandel

mandel: mandel.o bitmap.o
	gcc mandel.o bitmap.o -o mandel -lpthread -v -da -Q

series: series.o bitmap.o
	gcc series.o bitmap.o -o series -lpthread

master_mandel: master_mandel.o bitmap.o
	gcc master_mandel.o bitmap.o -o master_mandel -lpthread

mandel.o: mandel.c
	gcc -Wall -g -c mandel.c -o mandel.o

series.o: series.c
	gcc -Wall -g -c series.c -o series.o

master_mandel.o: master_mandel.c
	gcc -Wall -g -c master_mandel.c -o master_mandel.o

bitmap.o: bitmap.c
	gcc -Wall -g -c bitmap.c -o bitmap.o

clean:
	rm -f mandel.o series.o bitmap.o mandel series

 
