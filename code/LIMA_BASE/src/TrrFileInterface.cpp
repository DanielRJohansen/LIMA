//#include "ChemfilesInterface.h"
//#include "EngineUtils.cuh"
#include "MDFiles.h"

#include <string>
#include <format>

// definitions from the xdrfile library
#define TRR_MAGIC 1993
#define TRR_VERSION "GMX_trn_file"


class TrrWriter {
public:
    TrrWriter(const TrrWriter&) = delete;

    TrrWriter(const fs::path& path) {
        if (path.extension() != ".trr")
            throw std::runtime_error(std::format("Expected extension .trr, got '{}'", path.extension().string()));

        file.open(path, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error(std::format("Failed to open file at '{}'", path.string()));
        }
    }

    ~TrrWriter() {
        if (file.is_open()) {
            file.close();
        }
    }

    template<typename T>
    inline void write_as_big_endian(const T* data, size_t count) {
        const size_t byte_count = sizeof(T) * count;
        write_char(reinterpret_cast<const char*>(data), byte_count);
    }

    void write_u32(const uint32_t* data, size_t count) {
        write_as_big_endian(data, count);
    }

    void write_i32(const int32_t* data, size_t count) {
        write_as_big_endian(data, count);
    }

    void write_f32(const float* data, size_t count) {
        write_as_big_endian(data, count);
    }

    void write_single_u32(uint32_t value) {
        write_u32(&value, 1);
    }

    void write_single_i32(int32_t value) {
        write_i32(&value, 1);
    }

    void write_single_f32(float value) {
        write_f32(&value, 1);
    }

    void write_gmx_string(const std::string& value) {
        const uint32_t len = static_cast<uint32_t>(value.size() + 1);  // length with null terminator
        write_single_u32(len);
        write_opaque(value.c_str(), len - 1);
    }

    void write_opaque(const char* data, uint32_t count) {
        write_single_u32(count);
        write_char(data, count);
        const uint32_t num_filler = (4 - (count % 4)) % 4;
        const std::vector<char> filler(num_filler, 0);
        write_char(filler.data(), num_filler);
    }

    void write_char(const char* data, size_t count) {
        file.write(data, count);
        if (!file) {
            throw std::runtime_error(std::format(
                "Failed to write {} bytes to the file: {}",
                count, std::strerror(errno)
            ));
        }
    }

private:
    std::ofstream file;
};

struct FrameHeader {
    bool use_double;  /* Double precision?                                  */
    size_t ir_size;   /* Backward compatibility                             */
    size_t e_size;    /* Backward compatibility                             */
    size_t box_size;  /* Size in Bytes, non zero if a box is present        */
    size_t vir_size;  /* Backward compatibility                             */
    size_t pres_size; /* Backward compatibility                             */
    size_t top_size;  /* Backward compatibility                             */
    size_t sym_size;  /* Backward compatibility                             */
    size_t x_size;    /* Size in Bytes, non zero if coordinates are present */
    size_t v_size;    /* Size in Bytes, non zero if velocities are present  */
    size_t f_size;    /* Size in Bytes, non zero if forces are present      */

    size_t natoms; /* The total number of atoms                 */
    size_t step;   /* Current step number                       */
    size_t nre;    /* Backward compatibility                    */
    double time;   /* Current time (float or double)            */
    double lambda; /* Current value of lambda (float or double) */
};

void get_cell(std::vector<float>& box, Float3 boxSize) {
    box[0] = 1.f * boxSize.x;
    box[1] = 0.f;
    box[2] = 0.f;
    box[3] = 0.f;
    box[4] = 1.f * boxSize.y;
    box[5] = 0.f;
    box[6] = 0.f;
    box[7] = 0.f;
    box[8] = 1.f * boxSize.z;
}

void write_frame_header(TrrWriter& file, const FrameHeader& header) {
    file.write_single_i32(TRR_MAGIC);

    file.write_gmx_string(TRR_VERSION);

    // use_double is not written and has to be inferred when reading
    file.write_single_i32(static_cast<int32_t>(header.ir_size));
    file.write_single_i32(static_cast<int32_t>(header.e_size));
    file.write_single_i32(static_cast<int32_t>(header.box_size));
    file.write_single_i32(static_cast<int32_t>(header.vir_size));
    file.write_single_i32(static_cast<int32_t>(header.pres_size));
    file.write_single_i32(static_cast<int32_t>(header.top_size));
    file.write_single_i32(static_cast<int32_t>(header.sym_size));
    file.write_single_i32(static_cast<int32_t>(header.x_size));
    file.write_single_i32(static_cast<int32_t>(header.v_size));
    file.write_single_i32(static_cast<int32_t>(header.f_size));

    file.write_single_i32(static_cast<int32_t>(header.natoms));
    file.write_single_i32(static_cast<int32_t>(header.step));
    file.write_single_i32(static_cast<int32_t>(header.nre));
    file.write_single_f32(static_cast<float>(header.time));
    file.write_single_f32(static_cast<float>(header.lambda));
}

void WriteRow(TrrWriter& file, const std::vector<Float3>& positions, int32_t step, Float3 boxSize) {

    size_t boxBytesize = sizeof(float) * 3 * 3;

    const size_t natoms = positions.size();
    const size_t posdata_size = sizeof(float) * natoms * 3;;


    size_t veldata_size = 0;

    FrameHeader header = {
        false,    // use_double
        0,        // ir_size
        0,        // e_size
        boxBytesize, // box_size
        0,        // vir_size
        0,        // pres_size
        0,        // top_size
        0,        // sym_size
        posdata_size,   // x_size
        0,        // v_size
        0,        // f_size

        natoms,                                            // natoms
        static_cast<size_t>(step),                                      // step
        0,                                                 // nre
        0.f,//frame.get("time").value_or(0.0).as_double(),       // time
        0.f//frame.get("trr_lambda").value_or(0.0).as_double(), // lambda
    };
    write_frame_header(file,  header);

    std::vector<float> box(3 * 3);
    if (boxBytesize > 0) {
        get_cell(box, boxSize);
        file.write_f32(box.data(), 9);
    }

    for (auto& position : positions) {
        file.write_f32((float*)&position, 3);
    }
}

void MDFiles::TrrFile::Dump(const fs::path& path) const {
    TrrWriter outFile{ path };

    for (int i = 0; i < positions.size(); i++) {
		WriteRow(outFile, positions[i], i, boxSize);
	}
}

