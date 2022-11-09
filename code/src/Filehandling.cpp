#include "Filehandling.h"





bool Filehandler::ignoreRow(vector<char> ignores, string line) {
	if (line.length() == 0)
		return true;
	for (char c : ignores) {
		if (line[0] == c)
			return true;
	}
	return false;
}


vector<vector<string>> Filehandler::readFile(string path, int end_at, bool verbose) {
	cout << "Reading file " << path << endl;
	fstream file;
	file.open(path);


	vector<char> ignores = { ';', '/' };

	vector<vector<string>> rows;
	int row_cnt = 0;
	int ignore_cnt = 0;

	string line;
	while (getline(file, line)) {

		//cout << "line:" << line << endl;

		if (ignoreRow(ignores, line)) {
			ignore_cnt++;
			continue;
		}

		vector<string> row;
		stringstream ss(line);
		string element, element_nested;
		while (getline(ss, element, ' ')) {
			stringstream ss2 = stringstream(element);
			while (getline(ss2, element_nested, ';')) {
				if (element_nested.length() > 0)
					row.push_back(element_nested);
			}
		}

		rows.push_back(row);
		row_cnt++;

		if (row_cnt >= end_at)
			break;
	}



	if (verbose) {
		printf("%d rows read. %d rows ignored\n", row_cnt, ignore_cnt);
	}


	return rows;
}