#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <iostream>
#include <stdlib.h>

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
	int no_of_sentences = 0;
	int result = 0;
	std::string line;
	std::vector<std::string> temp;
	std::vector<int> length_of_each_sentence;
	int no_of_sentences_more_than_50 = 0;
	std::ifstream myfile ("english_wsj_dev.conll");
	//get the length of sentences
	if (myfile.is_open()) {
	    while (getline (myfile,line))
	    {
	      if (line.size() > 1) { //if the line is not empty
	      		temp = split(line, '\t'); //split the string with tab (\t) as its delimiter
	      		const char * c = (*(temp.begin())).c_str();
	      		result = atoi(c);
	      		if (result == 1) { //this indicates the beginning of a new sentence
	      			no_of_sentences++;
	      			length_of_each_sentence.push_back(1); 
	      		} else {
	      			length_of_each_sentence[length_of_each_sentence.size() - 1] = result;
	      		}
	      		temp.clear();
	      }
	    }
	    myfile.close();
	 } else {
  		std::cout << "Unable to open file"; 
  	}
	std::ifstream myfile_2 ("english_wsj_dev.conll"); //rewrite only the sentences that include less than 50 tokens
	if (myfile_2.is_open()) {
		int curr_sentence = -1;
		std::ofstream output_file ("english_wsj_dev_edited.conll");
		if (!output_file.is_open()) {
			std::cout << "fail to open output file" << std::endl;
		} else {
			while (getline (myfile_2,line))
			{
			  if (line.size() > 1) { //if the line is not empty
			  		temp = split(line, '\t'); //split the string with tab (\t) as its delimiter
			   		const char * c = (*(temp.begin())).c_str();
			   		result = atoi(c);
			   		if (result == 1) { //this indicates the beginning of a new sentence
			   			curr_sentence++;
			   			if (length_of_each_sentence[curr_sentence] <= 50) {
			   				output_file << '\n' << line << '\n';
			   			}
			   		} else {
			   			if (length_of_each_sentence[curr_sentence] <= 50) {
			   				output_file << line << '\n';
			   			}
			   		}
			   		temp.clear();
			  }
			}
			output_file.close();
			myfile_2.close();
		}
  	} else {
  		std::cout << "Unable to open file"; 
  	} 
	return 0;
}