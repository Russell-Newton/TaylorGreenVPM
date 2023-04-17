CC = nvc++
CFLAGS = -Wall -O2 -acc -Minfo=accel -std=c++17 
  
main: main.o vpm.o
	$(CC) $(CFLAGS) -o main main.o vpm.o 

main.o: main.cpp vpm.h 
	$(CC) $(CFLAGS) -c main.cpp
 
vpm.o: vpm.h
	$(CC) $(CFLAGS) -c vpm.cpp 

clean:
	rm -rf *.o 
