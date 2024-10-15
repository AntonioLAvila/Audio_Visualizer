compile: ./viz.cpp
	g++ ./viz.cpp -o app -lsndfile -lportaudio -lsfml-graphics -lsfml-window -lsfml-system -lfftw3

run:
	./app