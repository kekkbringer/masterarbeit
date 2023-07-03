#include <iostream>
#include <chrono>
#include <omp.h>
#include <cstdlib>
#include <string>

int main() {
	omp_set_num_threads(std::stoi(std::string(std::getenv("OMP_NUM_THREADS"))));
	std::cout << "Num of threads: " << omp_get_num_threads() << "\n";
}
