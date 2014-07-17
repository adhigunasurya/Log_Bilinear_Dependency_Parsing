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
	std::ofstream output("real_pairs_small.txt");
	if (!input.is_open() || !output.is_open()) {
		std::cout << "fail to open file" << std::endl;
	} else {
		std::string line;
		int no_of_iter = 0;
		int no_of_sentences = 1;
		while (getline(input, line)) {
			if (line.size() > 1) {
				no_of_iter++;
				output << line << std::endl;
			} else { //this is a space, break iteration
				no_of_sentences++;
				output << line << std::endl;
				if (no_of_iter > 40000) {
					break;
				}
			}
		}
		std::cout << "no_of_sentences = " << no_of_sentences << std::endl;
	}
	input.close();
	output.close();
	return 0;
}