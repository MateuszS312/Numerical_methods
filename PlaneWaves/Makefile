CC = g++
LIBS =
#-lboost_thread 'allegro-config --libs'
CFLAGS = -Wall -pedantic -O2 -std=c++17
program01: main.o billing.o histogram.o
		echo Buduję program...
		$(CC) -o program01 main.o billing.o histogram.o $(CFLAGS) $(LIBS)
main.o:	main.cpp billing.hpp histogram.hpp
		echo Kompiluję plik main.cpp...
		$(CC) -o main.o -c main.cpp $(CFLAGS)
billing.o:	billing.cpp billing.hpp
		echo Kompiluję plik 1.cpp...
		$(CC) -o billing.o -c billing.cpp $(CFLAGS)
histogram.o:	histogram.cpp histogram.hpp
		echo Kompiluję plik 2.cpp...
		$(CC) -o histogram.o -c histogram.cpp $(CFLAGS)
clean:
		rm *.o
