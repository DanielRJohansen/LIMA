#pragma once

#include "Bodies.cuh"
#include "CudaBuffer.h"

#include <algorithm>
#include <bit>
#include <map>
#include <type_traits>

// Groups have at most 64 atoms. Four 6-bit atom ids and an 8-bit parameter id
// fit in one coalesced 32-bit record. Larger tables use dense ids and wider indices.
template <typename Bond>
struct CompactBondTableDevice {
	static_assert(Bond::nAtoms <= 4);
	using Parameters = decltype(Bond::params);
	const std::array<uint8_t, Bond::nAtoms>* atomIds;
	const void* parameterIds;
	const Parameters* parameters;
	int parameterIdBytes;
	const uint32_t* packedBonds;

	__device__ void Load(int index, Bond& bond) const {
		if (packedBonds) {
			const uint32_t record = packedBonds[index];
			bond.params = parameters[record >> 24];
#pragma unroll
			for (int i = 0; i < Bond::nAtoms; i++) {
				const uint8_t id = static_cast<uint8_t>((record >> (i * 6)) & 63u);
				if constexpr (std::is_same_v<Bond, SingleBond>) bond.idInBondgroup[i] = id;
				else bond.atom_indexes[i] = id;
			}
			return;
		}
		const uint32_t parameterId = parameterIdBytes == 2 ? static_cast<const uint16_t*>(parameterIds)[index]
			: static_cast<const uint32_t*>(parameterIds)[index];
		bond.params = parameters[parameterId];
#pragma unroll
		for (int i = 0; i < Bond::nAtoms; i++) {
			if constexpr (std::is_same_v<Bond, SingleBond>) bond.idInBondgroup[i] = atomIds[index][i];
			else bond.atom_indexes[i] = atomIds[index][i];
		}
	}
};

template <typename Bond>
class CompactBondTable {
	using Parameters = decltype(Bond::params);
	CudaBuffer<std::array<uint8_t, Bond::nAtoms>> atomIds;
	CudaBuffer<uint16_t> parameterIds16;
	CudaBuffer<uint32_t> parameterIds32;
	CudaBuffer<Parameters> parameters;
	CudaBuffer<uint32_t> packedBonds;
	int parameterIdBytes = 1;
	size_t nParameters = 0;
	size_t nBonds = 0;

public:
	void SetData(const std::vector<Bond>& bonds) {
		// Compare float bit patterns: preserve signed zero and NaN payloads as well
		// as ordinary values, without relying on padding or approximate equality.
		using Key = std::array<uint32_t, sizeof(Parameters) / sizeof(uint32_t)>;
		static_assert(sizeof(Key) == sizeof(Parameters));
		static_assert(std::has_unique_object_representations_v<Key>);
		std::map<Key, uint32_t> parameterMap;
		std::vector<Parameters> parameterTable;
		std::vector<uint32_t> indices;
		std::vector<std::array<uint8_t, Bond::nAtoms>> ids;
		indices.reserve(bonds.size());
		ids.reserve(bonds.size());
		for (const Bond& bond : bonds) {
			const auto [entry, inserted] = parameterMap.try_emplace(std::bit_cast<Key>(bond.params), static_cast<uint32_t>(parameterTable.size()));
			if (inserted) parameterTable.push_back(bond.params);
			indices.push_back(entry->second);
			std::array<uint8_t, Bond::nAtoms> bondIds;
			if constexpr (std::is_same_v<Bond, SingleBond>) std::copy_n(bond.idInBondgroup, Bond::nAtoms, bondIds.begin());
			else std::copy_n(bond.atom_indexes, Bond::nAtoms, bondIds.begin());
			ids.push_back(bondIds);
		}
		parameters.SetData(parameterTable);
		nParameters = parameterTable.size();
		nBonds = bonds.size();
		if (parameterTable.size() <= 256) {
			parameterIdBytes = 1;
			std::vector<uint32_t> records(bonds.size());
			for (size_t index = 0; index < bonds.size(); index++) {
				uint32_t record = indices[index] << 24;
				for (int atom = 0; atom < Bond::nAtoms; atom++) {
					if (ids[index][atom] >= 64) throw std::invalid_argument("Compact bond atom id exceeds 64-particle group");
					record |= static_cast<uint32_t>(ids[index][atom]) << (atom * 6);
				}
				records[index] = record;
			}
			packedBonds.SetData(records);
		}
		else if (parameterTable.size() <= 65536) {
			parameterIdBytes = 2;
			atomIds.SetData(ids);
			const std::vector<uint16_t> compactIndices(indices.begin(), indices.end());
			parameterIds16.SetData(compactIndices);
		}
		else {
			parameterIdBytes = 4;
			atomIds.SetData(ids);
			parameterIds32.SetData(indices);
		}
	}

	CompactBondTableDevice<Bond> Get() const {
		const void* indices = parameterIdBytes == 2 ? static_cast<const void*>(parameterIds16.Get()) : static_cast<const void*>(parameterIds32.Get());
		return { atomIds.Get(), indices, parameters.Get(), parameterIdBytes, parameterIdBytes == 1 ? packedBonds.Get() : nullptr };
	}

	size_t ParameterCount() const { return nParameters; }
	size_t Bytes() const { return nBonds * (parameterIdBytes == 1 ? 4 : Bond::nAtoms + parameterIdBytes) + nParameters * sizeof(Parameters); }
};

struct CompactBondGroupsDevice {
	const BondGroup* groups;
	CompactBondTableDevice<SingleBond> singlebonds;
	CompactBondTableDevice<PairBond> pairbonds;
	CompactBondTableDevice<AngleUreyBradleyBond> anglebonds;
	CompactBondTableDevice<DihedralBond> dihedralbonds;
	CompactBondTableDevice<ImproperDihedralBond> improperdihedralbonds;
};
