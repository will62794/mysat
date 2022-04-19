# default:
# 	g++ -std=c++11 sat.cpp -o sat
test:
	g++ -std=c++11 test_sat.cpp easylogging++.cc -o test_sat
run: default
	./sat
format:
	clang-format -i sat.cpp test_sat.cpp