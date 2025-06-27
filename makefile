compile: ./viz.cpp
	g++ -std=c++23 ./viz.cpp -o app -lsndfile -lportaudio -lsfml-graphics -lsfml-window -lsfml-system -lfftw3

run:
	./app