#include "Benchmarks.h"

int main() {
	auto result = Benchmarks::STMV(Environment::Get(), EnvMode::Headless, 300).RunToCompletion();
	result.printStatus();
	return result.success ? 0 : 1;
}
