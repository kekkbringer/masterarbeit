#include "restart.hpp"

#include <fstream>
#include <string>

void setRestart(const int num) {
	const char* file = "restart";
	std::fstream(file, std::ios::out | std::ios::trunc) << num;
}

int readRestart() {
	const char* file = "restart";
	std::string s;
	std::fstream(file, std::ios::in) >> s;

	return stoi(s);
}
