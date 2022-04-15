default:
	g++ -std=c++11 sat.cpp -o sat
run: default
	./sat
format:
	clang-format -i sat.cpp