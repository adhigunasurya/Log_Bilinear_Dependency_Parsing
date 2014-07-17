#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <random>
#include <stdlib.h>

int main() {
	std::ifstream input("pairs.txt");
	std::ofstream output("pairs_without_space.txt");
	std::string line;
	if (!input.is_open() || !output.is_open()) {
		std::cout << "either input or output fails to open" << std::endl;
	} else {
		while (getline(input, line)) {
			if (line.size() > 1) {
				output << line << std::endl;
			}
		}
	}
	return 0;

}