#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <random>
#include <stdlib.h>
#include <cmath>
#include <utility>
#include <ctime>

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

template <typename T>
T StringToNumber (const std::string &Text )//Text not by const reference so that the function can be used with a 
{                               //character array as argument
	std::stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

template <typename T>
std::string NumberToString ( T Number )
{
	std::stringstream ss;
	ss << Number;
	return ss.str();
}


int main() {
	std::ifstream pairs("real_pairs_medium.txt");
	std::map<std::string, int> pairs_frequency_map;
	std::ifstream fine_pos_tags("fine_pos_tags.txt");
	std::string line;
	//initialise
	std::vector<std::string> fine_pos_tags_array;
	while (getline(fine_pos_tags, line)) {
		fine_pos_tags_array.push_back(line);
	}
	fine_pos_tags.close();
	for (int i = 0; i < fine_pos_tags_array.size(); i++) {
		for (int j = 0; j < fine_pos_tags_array.size(); j++) {
				std::string curr_word = fine_pos_tags_array[i] + " " + fine_pos_tags_array[j];
				pairs_frequency_map[curr_word] = 0;
		}
	}
	//get the frequencies
	while (getline(pairs, line)) {
		if (line.size() > 1) {
			std::vector<std::string> temp = split(line, ' ');
			pairs_frequency_map[temp[0]+" "+temp[1]] = pairs_frequency_map.at(temp[0]+" "+temp[1]) + 1;
		}
	}
	std::ofstream counts ("counts.txt");
	for (std::map<std::string, int>::iterator it = pairs_frequency_map.begin(); it != pairs_frequency_map.end(); it++) {
		counts << it->first << " " << it->second << std::endl;
	}
	counts.close();
	pairs.close();
	return 0;
}