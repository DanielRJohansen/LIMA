#include "Format.h"

std::string Lima::VFormat(std::string_view fmt, std::format_args args) {
	return std::vformat(fmt, args);
}
