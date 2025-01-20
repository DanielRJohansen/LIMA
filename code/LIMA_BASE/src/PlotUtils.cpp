#include "PlotUtils.h"
#include "Filehandling.h"

#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

void PlotUtils::PlotData(const std::vector<std::vector<float>>& vectors, const std::vector<std::string>& labels, const std::optional<std::vector<float>>& x_axis) {
	// Check if input vectors are valid
	if (vectors.empty() || vectors[0].empty() || labels.size() != vectors.size()) {
		std::cerr << "Error: Vectors or labels are invalid!" << std::endl;
		return;
	}

	size_t size = vectors[0].size();
	for (const auto& vec : vectors) {
		if (vec.size() != size) {
			std::cerr << "Error: All vectors must have the same size!" << std::endl;
			return;
		}
	}

	// Validate the x-axis, if provided
	if (x_axis.has_value() && x_axis->size() != size) {
		std::cerr << "Error: X-axis vector size must match the size of data vectors!" << std::endl;
		return;
	}

	
	const fs::path scriptPath = (FileUtils::GetLimaDir() / "dev" / "PyTools" / "Plotdata.py").string();
	const fs::path dataPath = FileUtils::GetLimaDir() / "dev" / "tmp.bin";

	// Write binary file
	std::ofstream outFile(dataPath, std::ios::binary);
	if (!outFile) {
		std::cerr << "Error: Could not open file for writing!" << std::endl;
		return;
	}

	// Write number of vectors
	int n = static_cast<int>(vectors.size());
	outFile.write(reinterpret_cast<const char*>(&n), sizeof(int));

	// Write size of vectors
	int vectorSize = static_cast<int>(vectors[0].size());
	outFile.write(reinterpret_cast<const char*>(&vectorSize), sizeof(int));

	// Write labels
	for (const auto& label : labels) {
		int labelLength = static_cast<int>(label.size());
		outFile.write(reinterpret_cast<const char*>(&labelLength), sizeof(int));
		outFile.write(label.c_str(), labelLength);
	}

	// Write optional independent variable flag
	int has_x_axis = x_axis.has_value() ? 1 : 0;
	outFile.write(reinterpret_cast<const char*>(&has_x_axis), sizeof(int));

	// Write the x-axis vector, if provided
	if (x_axis.has_value()) {
		outFile.write(reinterpret_cast<const char*>(x_axis->data()), x_axis->size() * sizeof(float));
	}

	// Write data vectors
	for (const auto& vec : vectors) {
		outFile.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(float));
	}

	outFile.close();

	// Call the Python script
	std::string command = "python " + scriptPath.string();
	std::system(command.c_str());
}