CC=g++
CFLAGS= -Wall -Ofast -std=c++11  -flto -pipe  -fopenmp -lz -Weffc++ -pedantic
LDFLAGS=-flto -lpthread -fopenmp -lz
EXEC=Test_indexes


all: $(EXEC)


Test_indexes:   Test_indexes.o NaiveIndex.o TestTable.o EvaluationMaker.o
	$(CC) -o $@ $^ $(LDFLAGS)

NaiveIndex.o: NaiveIndex.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

Test_indexes.o: Test_indexes.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

TestTable.o: TestTable.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

FastaParser.o: FastaParser.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

EvaluationMaker.o: EvaluationMaker.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
