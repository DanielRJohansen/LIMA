#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>
#include <type_traits>
#include <filesystem>



//[numSections:u32]
// section:
//[
//    [ nameLen:u32 ]
//    [name:char[nameLen]]
//    [numBytesOfData:u64]
//    [data:u8[numBytesOfData]]
//]


class UpgradeableFileFormat {
    std::ofstream out_;
    uint32_t numSections_ = 0;
public:
    UpgradeableFileFormat(const std::filesystem::path& path) {
        out_.open(path, std::ios::binary | std::ios::trunc);
        uint32_t header = 0;
        out_.write(reinterpret_cast<const char*>(&header), sizeof(header));
    }

    ~UpgradeableFileFormat() {
        out_.flush();
        out_.seekp(0, std::ios::beg);
        out_.write(reinterpret_cast<const char*>(&numSections_), sizeof(numSections_));
    }

    template<typename T>
    void WriteSection(const std::string& name, const std::vector<T>& vec) {
        uint32_t nameLen = static_cast<uint32_t>(name.size());
        out_.write(reinterpret_cast<const char*>(&nameLen), sizeof(nameLen));
        out_.write(name.data(), nameLen);

        uint64_t dataBytes = static_cast<uint64_t>(vec.size()) * sizeof(T);
        out_.write(reinterpret_cast<const char*>(&dataBytes), sizeof(dataBytes));
        out_.write(reinterpret_cast<const char*>(vec.data()), dataBytes);

        ++numSections_;
    }
};

//
//class UpgradeableFileParser {
//    struct Section {
//        std::string name;
//        std::vector<uint8_t> data;
//    };
//    std::vector<Section> sections_;
//
//public:
//    UpgradeableFileParser(const std::filesystem::path& path) {
//        std::ifstream in(path, std::ios::binary);
//        if (!in) throw std::runtime_error(std::format("Failed to open '{}'", path.string()));
//
//        uint32_t numSections;
//        in.read(reinterpret_cast<char*>(&numSections), sizeof(numSections));
//
//        for (uint32_t i = 0; i < numSections; ++i) {
//            uint32_t nameLen;
//            in.read(reinterpret_cast<char*>(&nameLen), sizeof(nameLen));
//
//            std::string name(nameLen, '\0');
//            in.read(name.data(), nameLen);
//
//            uint64_t dataBytes;
//            in.read(reinterpret_cast<char*>(&dataBytes), sizeof(dataBytes));
//
//            std::vector<uint8_t> data(dataBytes);
//            in.read(reinterpret_cast<char*>(data.data()), dataBytes);
//
//            sections_.push_back({ std::move(name), std::move(data) });
//        }
//    }
//
//    template<typename T>
//    std::vector<T> ReadSection(const std::string& name) const {
//        for (auto const& sec : sections_) {
//            if (sec.name == name) {
//                if (sec.data.size() % sizeof(T) != 0)
//                    throw std::runtime_error(std::format(
//                        "Section '{}' size ({}) is not a multiple of {}",
//                        name, sec.data.size(), sizeof(T)));
//
//                size_t count = sec.data.size() / sizeof(T);
//                std::vector<T> vec(count);
//                std::memcpy(vec.data(), sec.data.data(), sec.data.size());
//                return vec;
//            }
//        }
//        throw std::out_of_range(std::format("Section '{}' not found", name));
//    }
//
//    std::vector<std::string> GetSectionNames() const {
//        std::vector<std::string> names;
//        names.reserve(sections_.size());
//        for (auto const& sec : sections_) names.push_back(sec.name);
//        return names;
//    }
//};