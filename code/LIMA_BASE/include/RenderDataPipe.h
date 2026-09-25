#pragma once

#include "LimaTypes.cuh"

#include <atomic>
#include <cstddef>
#include <cstdint>

#include <cuda_runtime.h>

class RenderDataPipe {
public:
	enum class State {
		Uninitialized,
		Empty,
		Writing,
		Ready,
		Reading,
		Stopped
	};

	RenderDataPipe() = default;
	~RenderDataPipe();
	RenderDataPipe(const RenderDataPipe&) = delete;
	RenderDataPipe& operator=(const RenderDataPipe&) = delete;

	void Initialize(std::size_t positionCount);
	Float3* TryBeginWrite();
	void Publish(cudaStream_t producerStream, int64_t step);
	void CancelWrite();
	bool TryCopyToHost(Float3* destination, std::size_t destinationCount, int64_t& step);
	void Stop();

	std::size_t PositionCount() const { return positionCount; }
	State GetState() const { return state.load(std::memory_order_acquire); }

private:
	Float3* positions = nullptr;
	std::size_t positionCount = 0;
	cudaEvent_t writeCompleted = nullptr;
	cudaStream_t readStream = nullptr;
	int deviceId = 0;
	std::atomic<int64_t> publishedStep{ -1 };
	std::atomic<State> state{ State::Uninitialized };
};
