CC=g++
CFLAGS= -Wall -Ofast -std=c++11  -flto -pipe  -fopenmp -lz  -flto
LDFLAGS=-flto -lpthread -fopenmp -lz   -flto
EXEC=Test_indexes


all: $(EXEC)


Test_indexes:   Test_indexes.o NaiveIndex.o
	$(CC) -o $@ $^ $(LDFLAGS)

NaiveIndex.o: NaiveIndex.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

Test_indexes.o: Test_indexes.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
