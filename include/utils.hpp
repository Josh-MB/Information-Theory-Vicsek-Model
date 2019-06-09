#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <fmt/format.h>
#include <algorithm>
#include <random>
#include <cstring>
#include <chrono>

template <typename T>
class DoubleBuffer {
private:
	T* buffers[2];
	int active = 0;
	size_t sz;
public:
	DoubleBuffer(size_t sz) : sz(sz) {
		buffers[0] = (T*)calloc(sz, sizeof(T));
		buffers[1] = (T*)calloc(sz, sizeof(T));
	}
	DoubleBuffer(const DoubleBuffer<T>& other) : sz(other.sz), active(other.active) {
		for (int i = 0; i < 2; ++i) {
			buffers[i] = (T*)calloc(sz, sizeof(T));
			std::copy(other.buffers[i], other.buffers[i] + sz, buffers[i]);
		}
	}
	DoubleBuffer(DoubleBuffer<T>&& other) : sz{ std::move(other.sz) },
		active{ std::move(other.active) },
		buffers{ std::move(other.buffers) }
	{
		for (int i = 0; i < 2; ++i) {
			other.buffers[i] = nullptr;
		}
	}
	DoubleBuffer<T>& operator=(DoubleBuffer const& rhs) {
		if (this != &rhs) {
			free(buffers[0]);
			free(buffers[1]);
			sz = rhs.sz;
			active = rhs.active;
			buffers[0] = (T*)calloc(sz, sizeof(T));
			buffers[1] = (T*)calloc(sz, sizeof(T));
		}
		return *this;
	}

	DoubleBuffer<T>& operator=(DoubleBuffer && rhs) {
		if (this != &rhs) {
			free(buffers[0]);
			free(buffers[1]);
			sz = rhs.sz; rhs.sz = 0;
			active = rhs.active; rhs.active = 0;
			buffers = rhs.buffers;
			rhs.buffers[0] = rhs.buffers[1] = nullptr;
		}
		return *this;
	}

	~DoubleBuffer() {
		free(buffers[0]);
		free(buffers[1]);
	}
	/** The buffer that is full and ready for reading */
	T const * reader() const {
		return buffers[active];
	}
	/** The buffer that is ready for writing */
	T* writer() {
		return buffers[1 - active];
	}

	/** Swap the buffers. Invalidates handles into the buffers.  */
	void swap() {
		active = 1 - active;
	}

	size_t size() const {
		return sz;
	}
};

// Fills the writer with the data in reader
template <typename T>
void duplicate_reader_to_writer(DoubleBuffer<T> & buf) {
	std::memcpy(buf.writer(), buf.reader(), buf.size() * sizeof(T));
}


struct FlockState {
	FlockState(size_t sz) : h(sz), x(sz), y(sz) {}
	FlockState(FlockState const &other) = default;
	FlockState(FlockState &&other) = default;
	FlockState& operator=(FlockState const& rhs) = default;
	FlockState& operator=(FlockState && rhs) = default;

	DoubleBuffer<double> h;
	DoubleBuffer<double> x;
	DoubleBuffer<double> y;

	void swapBuffers() {
		x.swap(); y.swap(); h.swap();
	}
};

class FlockHistory {
private:
	size_t N;
	size_t U;
	double* h;
	double* x;
	double* y;
public:
	FlockHistory(size_t N, size_t U) : N{ N }, U{ U },
		h{ (double*const)calloc(N*U, sizeof(double)) },
		x{ (double*const)calloc(N*U, sizeof(double)) },
		y{ (double*const)calloc(N*U, sizeof(double)) }
	{

	}
	// No copy, too big
	FlockHistory(FlockHistory const &other) = delete;
	FlockHistory& operator=(FlockHistory const& rhs) = delete;

	FlockHistory(FlockHistory &&other) noexcept : N{ std::move(other.N) },
		U{ std::move(other.U) },
		h{ std::move(other.h) },
		x{ std::move(other.x) },
		y{ std::move(other.y) }
	{
		other.h = other.x = other.y = nullptr;
	}
	FlockHistory& operator=(FlockHistory && rhs) noexcept
	{
		if (this != &rhs) {
			free(h);
			free(x);
			free(y);
			N = rhs.N; rhs.N = 0;
			U = rhs.U; rhs.U = 0;
			h = rhs.h; rhs.h = nullptr;
			x = rhs.x; rhs.x = nullptr;
			y = rhs.y; rhs.y = nullptr;
		}
		return *this;
	}
	~FlockHistory() {
		free(h);
		free(x);
		free(y);
	}

	void recordState(FlockState const& fs, size_t u) {
		size_t fs_N = fs.h.size();
		auto hu = fs.h.reader(), xu = fs.x.reader(), yu = fs.y.reader();
		for (size_t i = 0; i < fs_N; ++i) {
			h[u + i * U] = hu[i];
			x[u + i * U] = xu[i];
			y[u + i * U] = yu[i];
		}
	}

	double const* read_h() const {
		return h;
	}

	double const* read_x() const {
		return x;
	}

	double const* read_y() const {
		return y;
	}
};

/** 
 Helper class for timing how long a scope takes to run
*/
class ScopeTimer {
private:
	std::string name;
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
public:
	ScopeTimer(std::string name) : name(name) {}
	~ScopeTimer() {
		auto end = std::chrono::high_resolution_clock::now();
		int64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		fmt::print("{} took {} ms\n", name, ms);
	}
	ScopeTimer(ScopeTimer const &other) = delete;
	ScopeTimer(ScopeTimer &&other) = delete;
	ScopeTimer& operator=(ScopeTimer const& rhs) = delete;
	ScopeTimer& operator=(ScopeTimer && rhs) = delete;
};

// Returns random number if seed==0, otherwise returns seed
uint64_t calcSeed(uint64_t seed);

/**
 * Following functions and macros thanks to Dr. Lionel Barnett.
 * Modified by Joshua Brown
 */
#define ERRPT       fmt::print(stderr,"ERROR in '{}' [{}:{}]: ",__FUNCTION__,__FILE__,__LINE__)
#define EEXIT(...)  {ERRPT; fmt::print(stderr,__VA_ARGS__); fputc('\n',stderr); exit(EXIT_FAILURE);}
#define PEEXIT(...) {ERRPT; fmt::print(stderr,__VA_ARGS__); fputc('\n',stderr); perror(NULL); exit(EXIT_FAILURE);}

void magic_header(FILE* fp, const char* const hmessage);
bool check_magic_header(FILE* fp);
void progrep(const char* const msg, const size_t u, const size_t U);
void progrep(FILE * fp, const char* const msg, const size_t u, const size_t U);

/**
 * End Barnett functions
 */

#ifdef WIN32
inline void sincos(double h, double* s, double* c) {
	*s = sin(h);
	*c = cos(h);
}
inline FILE* popen(const char* command, const char* mode) {
	return _popen(command, mode);
}

inline int pclose(FILE* stream) {
	return _pclose(stream);
}
#endif

void make_dir(std::string path);

#if (defined _MSVC_LANG && _MSVC_LANG < 201402L) ||  (!defined _MSVC_LANG && __cplusplus < 201402L)
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>
/**
 * Provide the make_unique helper functions for C++11.
 * See: https://stackoverflow.com/questions/17902405/how-to-implement-make-unique-function-in-c11
 */
namespace std {
	template<class T> struct _Unique_if {
		typedef unique_ptr<T> _Single_object;
	};

	template<class T> struct _Unique_if<T[]> {
		typedef unique_ptr<T[]> _Unknown_bound;
	};

	template<class T, size_t N> struct _Unique_if<T[N]> {
		typedef void _Known_bound;
	};

	template<class T, class... Args>
	typename _Unique_if<T>::_Single_object
		make_unique(Args&& ... args) {
		return unique_ptr<T>(new T(std::forward<Args>(args)...));
	}

	template<class T>
	typename _Unique_if<T>::_Unknown_bound
		make_unique(size_t n) {
		typedef typename remove_extent<T>::type U;
		return unique_ptr<T>(new U[n]());
	}

	template<class T, class... Args>
	typename _Unique_if<T>::_Known_bound
		make_unique(Args&& ...) = delete;
}
#endif
#endif // __UTILS_HPP__