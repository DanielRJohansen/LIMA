#pragma once

#include <algorithm>
#include <exception>
#include <execution>
#include <mutex>
#include <ranges>

namespace ParallelUtils {
	// Calls f(i) for i in [0, n) in parallel.
	// Exceptions may not escape a parallel algorithm (std::terminate), so we carry the first one out and rethrow it
	template <typename F>
	void ParallelFor(size_t n, F&& f) {
		std::mutex mutex;
		std::exception_ptr error;
		auto indices = std::views::iota(size_t{ 0 }, n);
		std::for_each(std::execution::par, indices.begin(), indices.end(), [&](size_t i) {
			try { f(i); }
			catch (...) {
				std::scoped_lock lock(mutex);
				if (!error) error = std::current_exception();
			}
			});
		if (error)
			std::rethrow_exception(error);
	}

	// Calls f(i) for i in [0, n), in parallel blocks when n is large enough that the parallel dispatch is worth it
	template <typename F>
	void ParallelForBlocked(size_t n, F&& f, size_t minNForParallel = 16384, size_t blockSize = 4096) {
		if (n < minNForParallel) {
			for (size_t i = 0; i < n; i++) f(i);
			return;
		}
		ParallelFor((n + blockSize - 1) / blockSize, [&](size_t block) {
			const size_t end = std::min(n, (block + 1) * blockSize);
			for (size_t i = block * blockSize; i < end; i++) f(i);
			});
	}
}
