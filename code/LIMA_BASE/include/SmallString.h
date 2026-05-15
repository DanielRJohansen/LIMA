#include <array>
#include <cstddef>
#include <cstring>
#include <format>
#include <stdexcept>
#include <string_view>

class SmallString {
public:
	static constexpr std::size_t maxSize = 8;

	constexpr SmallString() = default;

	constexpr SmallString(std::string_view text) {
		//Assign(text);
		if (text.size() > maxSize)
			throw std::length_error("SmallString can hold at most 8 chars");
		//strncpy(data_.data(), text.data(), maxSize);
		for (std::size_t i = 0; i < text.size(); ++i)
			data_[i] = text[i];
	}

	constexpr std::size_t Size() const {
		for (std::size_t i = 0; i < maxSize; ++i) {
			if (data_[i] == '\0')
				return i;
		}
		return maxSize;
	}

	constexpr bool Empty() const {
		return data_[0] == '\0';
	}

	constexpr std::size_t Capacity() const {
		return maxSize;
	}

	std::array<char, maxSize>& Data() {
		return data_;
	}
	constexpr std::string_view View() const {
		return { data_.data(), Size()};
	}

	constexpr char operator[](std::size_t index) const {
		return data_[index];
	}

	constexpr char& operator[](std::size_t index) {
		return data_[index];
	}

	constexpr void Clear() {		
		data_[0] = '\0';
	}

	/*constexpr void PushBack(char c) {
		if (size_ == maxSize)
			throw std::length_error("SmallString is full");

		data_[size_++] = c;
		data_[size_] = '\0';
	}*/

	friend constexpr bool operator==(const SmallString& a, const SmallString& b) {
		return a.View() == b.View();
	}
	friend constexpr bool operator==(const SmallString& a, std::string_view b) {
		return a.View() == b;
	}


private:
	std::array<char, maxSize> data_ = {'\0', '\0' , '\0' , '\0','\0', '\0' , '\0' , '\0' };
};

template <>
struct std::formatter<SmallString> : std::formatter<std::string_view> {
	auto format(const SmallString& str, std::format_context& ctx) const {
		return std::formatter<std::string_view>::format(str.View(), ctx);
	}
};