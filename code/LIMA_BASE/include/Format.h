#pragma once

// Use Lima::Format instead of std::format. std::format instantiates the whole formatting engine in every translation
// unit that calls it (~3 s of compile time and ~100 KB of code per file at -O3). Lima::Format type-erases the
// arguments, so the engine is compiled once, in Format.cpp. Format strings are still checked at compile time.

#include <format>
#include <string>
#include <string_view>

namespace Lima {
	std::string VFormat(std::string_view fmt, std::format_args args);

	template<typename... Args>
	std::string Format(std::format_string<Args...> fmt, Args&&... args) {
		return VFormat(fmt.get(), std::make_format_args(args...));
	}
}
