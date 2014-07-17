#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

std::map<std::string, int> postag_mapping(std::string file_name) {
	std::string line;
	std::map<std::string, int> result;
	std::ifstream fine_pos_tags (file_name);
	int j = -1;
	if (!fine_pos_tags.is_open()) {
		std::cout << "Error: file fails to open" << std::endl;
	} else {
		while (getline(fine_pos_tags, line)) {
			j++;
			result[line] = j;
		}
		j++;
		result["START"] = j++;
		result["END"] = j++;
		result["ROOT"]  = j++;
		result["ROOT+"] = j++;
		result["ROOT-"] = j++; 
	}
	fine_pos_tags.close();
	return result;
}

int main() {
	std::map<std::string, int> mapping = postag_mapping("fine_pos_tags.txt");
	std::cout << mapping["START"] << std::endl;
	return 0;
}