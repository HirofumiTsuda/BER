objs = main.o correlation.o code.o replace.o
objs2 = db.o correlation.o code.o replace.o
objs3 = db_dat.o correlation.o code.o replace.o
objs4 = db_QAM.o correlation.o code.o replace.o

main: $(objs)
	g++ -fopenmp -O3 -o main $(objs) -lm

db: $(objs2)
	g++ -fopenmp -O3 -o db $(objs2) -lm

db_dat: $(objs3)
	g++ -fopenmp -O3 -o db_dat $(objs3) -lm

db_QAM: $(objs4)
	g++ -fopenmp -O3 -o db_QAM $(objs4) -lm

main.o: main.cpp
	g++ -c -fopenmp main.cpp

db.o: db.cpp
	g++ -c -fopenmp db.cpp

db_dat.o: db_dat.cpp
	g++ -c -fopenmp db_dat.cpp

db_QAM.o: db_QAM.cpp
	g++ -c -fopenmp db_QAM.cpp

correlation.o: correlation.cpp
	g++ -c correlation.cpp

code.o: code.cpp
	g++ -c code.cpp

replace.o: replace.cpp
	g++ -c replace.cpp

.PHONY: clean
clean:
	rm -f main $(objs)
