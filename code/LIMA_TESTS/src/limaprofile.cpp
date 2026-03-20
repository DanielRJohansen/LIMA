#include "Benchmarks.h"

int main() {
	//Benchmarks::Benchmark("membrane20", "membranesolvated_em", 2);
	//Benchmarks::Psome(EnvMode::ConsoleOnly, 5);

	//Benchmarks::Benchmark("stmv", std::nullopt, 2);
	Benchmarks::STMV();
	return 0;
}