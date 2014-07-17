#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <map>

#ifndef DUMMY_HEAD
#define DUMMY_HEAD -1
#endif

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

template <class T>
bool contains(T element, std::vector<T> vect) {
	for (typename std::vector<T>::iterator it = vect.begin(); it != vect.end(); it++) {
		if (element == *it) {
			return true;
		}
	}
	return false;
}

int main() {
    std::ifstream myfile("output_file.txt");
    std::vector<std::string> curr_sentence;
    std::vector<int> heads;
    std::string line;
    std::vector<std::string> temp_heads;
    std::vector<std::string> temp;
    std::ofstream output("pairs.txt");
    std::ofstream output_distance_count("distance_count.txt");
    std::map<int, int> distance_count;
    if (myfile.is_open() && output.is_open()) {
        while (getline (myfile,line)) {
            temp = split(line, '\t');
            curr_sentence = split(temp[0], ' ');
            temp.clear();
            temp = split(split(line, '\t')[1], ' ');
            for (std::vector<std::string>::iterator it = temp.begin(); it != temp.end(); it++) {
                heads.push_back(atoi(it -> c_str()));
            }
            temp.clear();
            int i = -1;
            for (std::vector<std::string>::iterator it = curr_sentence.begin(); it != curr_sentence.end(); it++) {
                i++;
                if (*it != "START" && *it != "END") {
                    //the first component is the modifier's pos tag
                    output << *it << ' ';
                    //the second component is the head's pos tag
                    int head_index = heads[i];
                    if (head_index == 0) {
                        output << "ROOT" << ' ';
                    } else {
                        output << curr_sentence[head_index] << ' ';
                    }
                    //the third component is the pos tag of the token before the modifier
                    output << *(it - 1) << ' ';
                    //the fourth component is the pos tag of the token after the modifier
                    output << *(it + 1) << ' ';
                    //the fifth component is the pos tag of the token before the parent
                    if (head_index == 0) {
                        output << "ROOT-" << ' ';
                    } else {
                        output << curr_sentence[head_index - 1] << ' ';
                    }
                    //the sixth component is the pos tag of the token after the parent
                    if (head_index == 0) {
                        output << "ROOT+" << ' ';
                    } else {
                         output << curr_sentence[head_index + 1] << ' ';
                    }
                    //the seventh component is the distance between the modifier and the head
                    std::map<int, int>::iterator count_iter;
                    count_iter = distance_count.find(abs(i - head_index));
                    if (count_iter == distance_count.end()) {
                        distance_count[abs(i - head_index)] = 1;
                    } else {
                         distance_count[abs(i - head_index)] = distance_count[abs(i - head_index)] + 1;
                    }
                    output << abs(i - head_index) << '\n';
                }
            }
            output << '\n';
            heads.clear();
            temp_heads.clear();
            curr_sentence.clear();
        }
        for (std::map<int, int>::iterator it = distance_count.begin(); it != distance_count.end(); it++) {
            output_distance_count << "Distance " << it->first << " occurs " << it->second << " times" << '\n'; 
        }
        myfile.close();
        output.close();
    } else {
        std::cout << "Either the input or the output file fails to open" << std::endl;
    }

    return 0;
}
