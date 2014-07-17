#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

int main() {
	std::ifstream input("english_wsj_dev_edited_small.txt");
	std::ofstream parser_input("parser_input.txt");
	std::ofstream gold_standard("gold_standard.txt");
	std::string line;
	if (!input.is_open() || !parser_input.is_open() || !gold_standard.is_open()) {
		std::cout << "Input or output files fail to open" << std::endl;
	} else {
		while (getline(input, line)) {
			if (line.size() <= 1) {
				parser_input << std::endl;
				gold_standard << std::endl;
			} else {
				std::vector<std::string> temp = split(line, '\t');
				gold_standard << temp[6] << " " << temp[0] << std::endl;
				parser_input << temp[4] << " ";
			}
		}
		input.close();
		parser_input.close();
		gold_standard.close();
	}
	return 0;
}