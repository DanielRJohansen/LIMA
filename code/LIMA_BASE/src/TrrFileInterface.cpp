//#include "ChemfilesInterface.h"
//#include "EngineUtils.cuh"
#include "MDFiles.h"

#include <string>
//#include <format>

// definitions from the xdrfile library
#define TRR_MAGIC 1993
#define TRR_VERSION "GMX_trn_file"

struct ShittyFileFormat {
    ShittyFileFormat(const ShittyFileFormat&) = delete;
    ShittyFileFormat(const std::string& path) {
        file = fopen(path.c_str(), "wb");
    }

    ~ShittyFileFormat() {
        fclose(file);
    }

    FILE* file;

    template<typename T>
    inline void write_as_big_endian(const T* data, size_t count) {
        const size_t byte_count = sizeof(T) * count;
        this->write_char(reinterpret_cast<const char*>(data), byte_count);

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

    void write_gmx_string(const std::string & value) {
        // lenght with null terminator
        const uint32_t len = static_cast<uint32_t>(value.size() + 1);
        write_single_u32(len);
        // next comes XDR string without terminator
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
        auto written = std::fwrite(data, 1, count, file);
        if (written != count) {
        /*    throw std::format(
                "failed to write {} bytes to the file at '{}': {}",
                count, this->path(), std::strerror(errno)
            );*/
        }
    }
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

class TRRFormat {


    ShittyFileFormat file;

    void get_cell(std::vector<float>& box) {
        throw std::runtime_error("Cannot access BOX_LEN, dev please fix");
        //box[0] = 1.f * BOX_LEN_NM;
        //box[1] = 0.f;
        //box[2] = 0.f;
        //box[3] = 0.f;
        //box[4] = 1.f * BOX_LEN_NM;
        //box[5] = 0.f;
        //box[6] = 0.f;
        //box[7] = 0.f;
        //box[8] = 1.f * BOX_LEN_NM;
    }

    void write_frame_header(const FrameHeader& header) {
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

public:
    TRRFormat(const TRRFormat&) = delete;
    TRRFormat(const std::string& path) : file(ShittyFileFormat{ path }) {
    }

    void write(const std::vector<Float3>& positions, int32_t step) {

        size_t box_size = sizeof(float) * 3 * 3;
        
        const size_t natoms = positions.size();
        const size_t posdata_size = sizeof(float) * natoms * 3;;
        

        size_t veldata_size = 0;


        FrameHeader header = {
            false,    // use_double
            0,        // ir_size
            0,        // e_size
            box_size, // box_size
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
        write_frame_header(header);

        std::vector<float> box(3 * 3);
        if (box_size > 0) {
            get_cell(box);
            file.write_f32(box.data(), 9);
        }

        for (auto& position : positions) {
            file.write_f32((float*)&position, 3);
        }
    }
};





void MDFiles::TrrFile::dumpToFile(const Simulation* sim, const std::string& path) {
    TRRFormat trrfile(path);

    const int inc = sim->getStep() > 2000 ? 10 : 1;
    

    std::vector<Float3> positions(sim->boxparams_host.total_particles);
    // TODO: Fix the below, by passing a lambda that gets the entryindex based on the step
    //for (int entryindex = 0; entryindex < LIMALOGSYSTEM::getDataentryIndex(sim->getStep()); entryindex += inc) {

    //    int index = 0; 

    //    // Todo: this can be optimized with some insert magic, but i do not have the brain capacity. Ask gpt?


    //    for (int compound_id = 0; compound_id < sim->boxparams_host.n_compounds; compound_id++) {
    //        for (int i = 0; i < sim->box_host->compounds[compound_id].n_particles; i++) {
    //            positions[index++] = sim->traj_buffer->getCompoundparticleDatapointAtIndex(compound_id, i, entryindex);
    //        }
    //    }

    //    for (int solvent_index = 0; solvent_index < sim->boxparams_host.n_solvents; solvent_index++) {
    //        positions[index++] = sim->traj_buffer->getSolventparticleDatapointAtIndex(solvent_index, entryindex);
    //    }

    //    trrfile.write(positions, entryindex);
    //}
}
