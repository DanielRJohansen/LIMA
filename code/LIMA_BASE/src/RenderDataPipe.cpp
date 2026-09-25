#include "RenderDataPipe.h"

#include <stdexcept>

namespace {
	void CheckCuda(cudaError_t result, const char* operation) {
		if (result != cudaSuccess)
			throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(result));
	}
}

RenderDataPipe::~RenderDataPipe() {
	Stop();
	if (!positions && !writeCompleted && !readStream)
		return;
	cudaSetDevice(deviceId);
	if (readStream) {
		cudaStreamSynchronize(readStream);
		cudaStreamDestroy(readStream);
	}
	if (writeCompleted) {
		cudaEventSynchronize(writeCompleted);
		cudaEventDestroy(writeCompleted);
	}
	if (positions)
		cudaFree(positions);
}

void RenderDataPipe::Initialize(std::size_t count) {
	State expected = State::Uninitialized;
	if (!state.compare_exchange_strong(expected, State::Writing, std::memory_order_acq_rel))
		throw std::logic_error("RenderDataPipe initialized more than once");
	if (count == 0) {
		positionCount = 0;
		state.store(State::Empty, std::memory_order_release);
		return;
	}
	CheckCuda(cudaGetDevice(&deviceId), "Could not get CUDA device for render data pipe");
	CheckCuda(cudaMalloc(&positions, count * sizeof(Float3)), "Could not allocate render data pipe");
	CheckCuda(cudaEventCreateWithFlags(&writeCompleted, cudaEventDisableTiming), "Could not create render data event");
	CheckCuda(cudaStreamCreateWithFlags(&readStream, cudaStreamNonBlocking), "Could not create render data read stream");
	positionCount = count;
	state.store(State::Empty, std::memory_order_release);
}

Float3* RenderDataPipe::TryBeginWrite() {
	State expected = State::Empty;
	if (!state.compare_exchange_strong(expected, State::Writing, std::memory_order_acq_rel))
		return nullptr;
	return positions;
}

void RenderDataPipe::Publish(cudaStream_t producerStream, int64_t step) {
	if (state.load(std::memory_order_acquire) != State::Writing)
		throw std::logic_error("Publishing render data without an active write");
	CheckCuda(cudaEventRecord(writeCompleted, producerStream), "Could not record render data event");
	publishedStep.store(step, std::memory_order_relaxed);
	state.store(State::Ready, std::memory_order_release);
}

void RenderDataPipe::CancelWrite() {
	State expected = State::Writing;
	state.compare_exchange_strong(expected, State::Empty, std::memory_order_release, std::memory_order_relaxed);
}

bool RenderDataPipe::TryCopyToHost(Float3* destination, std::size_t destinationCount, int64_t& step) {
	State expected = State::Ready;
	if (!state.compare_exchange_strong(expected, State::Reading, std::memory_order_acq_rel))
		return false;
	try {
		if (destinationCount < positionCount)
			throw std::invalid_argument("Render destination is smaller than the render data pipe");
		CheckCuda(cudaSetDevice(deviceId), "Could not select render data CUDA device");
		CheckCuda(cudaStreamWaitEvent(readStream, writeCompleted), "Could not wait for render data");
		CheckCuda(cudaMemcpyAsync(destination, positions, positionCount * sizeof(Float3), cudaMemcpyDeviceToHost, readStream),
			"Could not copy render data to host");
		CheckCuda(cudaStreamSynchronize(readStream), "Could not finish copying render data to host");
		step = publishedStep.load(std::memory_order_relaxed);
		state.store(State::Empty, std::memory_order_release);
		return true;
	}
	catch (...) {
		state.store(State::Empty, std::memory_order_release);
		throw;
	}
}

void RenderDataPipe::Stop() {
	auto current = state.load(std::memory_order_acquire);
	while (current != State::Stopped
		&& !state.compare_exchange_weak(current, State::Stopped, std::memory_order_acq_rel)) {}
}

void RenderDataPipe::SetStatus(const SimStatus& value, bool isCompleted) {
	const std::lock_guard lock(statusMutex);
	status = value;
	completed = isCompleted;
}

SimStatus RenderDataPipe::GetStatus(bool& isCompleted) const {
	const std::lock_guard lock(statusMutex);
	isCompleted = completed;
	return status;
}
