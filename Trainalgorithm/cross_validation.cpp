#include <Eigen/Dense>
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
#include "graph.cpp" 

#ifndef FEATURE_SIZE
#define FEATURE_SIZE 10
#endif

#ifndef NO_OF_POSTAGS
#define NO_OF_POSTAGS 47
#endif

#ifndef NO_OF_POSTAGS_WITHOUT_DUMMY
#define NO_OF_POSTAGS_WITHOUT_DUMMY 45
#endif

#ifndef MAX_DISTANCE
#define MAX_DISTANCE 49
#endif

#ifndef NO_OF_DISTANCE
#define NO_OF_DISTANCE 28
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

double get_validation_accuracy(std::vector<Eigen::ArrayXXd> Ci, Eigen::ArrayXXd R, std::string validation_input_file, std::string output, std::map<std::string, int> postag_mapping, std::map<std::string, int> big_postag_mapping, Eigen::ArrayXd bv) {
	std::ifstream input_file(validation_input_file);
	std::ofstream output_file(output);
	if (!input_file.is_open() || !output_file.is_open()) {
		std::cout << "either input or output file fails to open" << std::endl;
	} else {
		std::string line;
		//get the input sentences from the validation_input_file
		int j = -1;
		while (getline(input_file, line)) {
			std::vector<std::string> sentence;
			sentence.push_back("START");
			std::vector<std::string> temp = split(line, ' ');
			int length = temp.size();
			for (std::vector<std::string>::iterator it = temp.begin(); it != temp.end(); it++) {
				sentence.push_back(*it);
			}
			sentence.push_back("END");
			//construct a graph with the given number of tokens
			//Graph G(length); 
			//get the probabilities between each token
			for (int i = 0; i < sentence.size(); i++) { //i indicates all the possible modifiers
				for (int j = 0; j < sentence.size(); j++) { //j indicates all the possible heads
					if (sentence[i] != "START" && sentence[i] != "END" && sentence[j] != "START" && sentence[j] != "END" && i != j) {
						//calculate the probability
						int dist = abs(i-j);
						std::ostringstream s;
						s << dist;
						std::string distance = s.str();
						if (big_postag_mapping.find(sentence[i]+distance) != big_postag_mapping.end()) {
							double bias_term = bv(big_postag_mapping[sentence[i] + distance]);
						} else {
							std::cout << "postag not found" << std::endl;
						}
						//connect the graph
					}
				}
			}
			
			//get the MST
		}
		//calculate accuracy (ignoring punctuations)
		input_file.close();
	}
	return 0.0;
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

void get_parameters_without_space(std::string input_file, std::string output_file) {
	std::ifstream input(input_file);
	std::ofstream output(output_file);
	if (!input.is_open() || !output.is_open()) {
		std::cout << "Fail to remove space due to file I/O failures" << std::endl;
	} else {
		std::string line;
		while (getline(input, line)) {
			if (line.size() > 1) {
				std::vector<std::string> temp = split(line, ' ');
				std::vector<std::string> temp_without_space;
				for (std::vector<std::string>::iterator it = temp.begin(); it != temp.end(); it++) {
					if (it -> size() > 1) {
						temp_without_space.push_back(*it);
					}
				}
				for (std::vector<std::string>::iterator it = temp_without_space.begin(); it != temp_without_space.end(); it++) {
					output << *it << std::endl;
				}
			} else {
				output << std::endl;
			}
		}
		input.close();
		output.close();
	}
}

void getMapping(std::map<std::string, int>& big_postag_map, std::map<std::string, int>& small_postag_map, std::string big_postag_file, std::string small_postag_file) {
	std::ifstream big_input(big_postag_file);
	std::ifstream small_input(small_postag_file);
	if (!big_input.is_open() || !small_input.is_open()) {
		std::cout << "Either the big or the small postag mapping cannot be opened" << std::endl;
	} else {
		std::string line;
		while (getline(big_input, line)) {
			std::vector<std::string> temp = split(line, ' ');
			big_postag_map[temp[0]] = StringToNumber<int> (temp[1]);
		}
		while (getline(small_input, line)) {
			std::vector<std::string> temp = split(line, ' ');
			small_postag_map[temp[0]] = StringToNumber<int> (temp[1]);
		}
	}
}

double get_tree_score(std::vector<Edge>& MST) {
	double result = 0.0;
	for (std::vector<Edge>::iterator it = MST.begin(); it != MST.end(); it++) {
		result += it->getWeight();
	}
	return result;
}

int main() {
	//first part of the program: get the trained parameters from an external file
	std::ofstream validation_pairs("validation_pairs_5.txt");
	get_parameters_without_space("cross_validate_5.txt", "validation_pairs_5_removed_space.txt");
	std::ifstream parameters("validation_pairs_5_removed_space.txt");
	if (!parameters.is_open() || !validation_pairs.is_open()) {
		std::cout << "Fail to open either input or output file" << std::endl;
	} else {
		//initialise matrices for the parameters
		Eigen::ArrayXXd R = Eigen::ArrayXXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS);
		Eigen::ArrayXXd Q = Eigen::ArrayXXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS);
		//Eigen::ArrayXd word_bias_pair = Eigen::ArrayXd::Zero(NO_OF_POSTAGS_WITHOUT_DUMMY * NO_OF_DISTANCE);
		std::vector<Eigen::ArrayXXd> Ci;
		for (int i = 0; i < 5; i++) {
			Ci.push_back(Eigen::ArrayXXd::Zero(FEATURE_SIZE, FEATURE_SIZE));
		}
		//First, get all the parameters from the "parameters" file
		int curr_parameter = 0;
		std::string line;
		std::ofstream check_input("check_input.txt");
		//get C1, C2, ..., C5
		while (curr_parameter <= 4) {
			getline(parameters, line);
			if (line.size() > 1) {
				for (int i = 0; i < Ci[curr_parameter].rows(); i++) {
					for (int j = 0; j < Ci[curr_parameter].cols(); j++) {
						if (i == 0 && j == 0) {
							Ci[curr_parameter](i, j) = StringToNumber<double> (line);
						} else {
							getline(parameters, line);
							Ci[curr_parameter](i, j) = StringToNumber<double> (line);
						}
					}
				}
				check_input << Ci[curr_parameter] << std::endl;
				check_input << std::endl;
				curr_parameter++;
			}
		}
		//get R
		getline(parameters, line);
		while (line.size() <= 1) {
			getline(parameters, line);
		}
		for (int i = 0; i < R.rows(); i++) {
			for (int j = 0; j < R.cols(); j++) {
				if (i == 0 && j == 0) {
					R(i, j) = StringToNumber<double> (line);
				} else {
					getline(parameters, line);
					R(i, j) = StringToNumber<double> (line);
				}
			}
		}
		check_input << R << std::endl;
		check_input << std::endl;
		curr_parameter++;
		//get Q
		getline(parameters, line);
		while (line.size() <= 1) {
			getline(parameters, line);
		}
		for (int i = 0; i < Q.rows(); i++) {
			for (int j = 0; j < Q.cols(); j++) {
				if (i == 0 && j == 0) {
					Q(i, j) = StringToNumber<double> (line);
				} else {
					getline(parameters, line);
					Q(i, j) = StringToNumber<double> (line);
				}
			}
		}
		check_input << Q << std::endl;
		check_input << std::endl;
		curr_parameter++;
		//get bias
		/*
		getline(parameters, line);
		while (line.size() <= 1) {
			getline(parameters, line);
		}
		for (int i = 0; i < word_bias_pair.rows(); i++) {
			if (i == 0) {
				word_bias_pair(i) = StringToNumber<double> (line);
			} else {
				getline(parameters, line);
				word_bias_pair(i) = StringToNumber<double> (line);
			}
		} */
		//check_input << word_bias_pair << std::endl;
		//check_input << std::endl;
		//curr_parameter++;
		//close the files that are no longer necessary
		check_input.close();
		parameters.close();
		//get sentences from the cross validation file
		std::ifstream sentences("dev_parser_input_small.txt");
		std::ofstream final_result("parsing_result_cross_validation.txt");
		std::ofstream final_graph("final_graph.txt");
		if (!sentences.is_open()) {
			std::cout << "Fail to open the input from the cross validation file" << std::endl;
		} else {
			//MAIN LOOP
			std::map<std::string, int> big_pos_tag_map;
			std::map<std::string, int> small_pos_tag_map;
			getMapping(big_pos_tag_map, small_pos_tag_map, "big_pos_tag_map.txt", "small_pos_tag_map.txt");
			Eigen::VectorXd prediction = Eigen::VectorXd::Zero(FEATURE_SIZE);
			int no_of_sentence = 0;
			while (getline(sentences, line)) {
				no_of_sentence++;
				std::vector<std::string> temp = split(line, '\t');
				//get only the sentence and ignore the distances
				std::vector<std::string> sentence = split(temp[0], ' ');
				double min = 9999.0;
				std::vector<int> punctuations;
				for (int j = 0; j < sentence.size(); j++) {
					//if punctuations, add to the punctuations vector 
					if (sentence[j] == "``" || sentence[j] == "\'\'" || sentence[j] == "(" || sentence[j] == ")" || sentence[j] == "," || sentence[j] == "." || sentence[j] == ":" || sentence[j] == "$" || sentence[j] == "#") {
						punctuations.push_back(j);
					}
				}
				//ROOT = MIN_TOKEN
				//assign the probabilities of each possible modifier given all other tokens using the pre-defined parameters
				double min_MST = 999999.9;
				Graph * min_graph = nullptr;
				for (int i = 0; i < sentence.size(); i++) {
					if (sentence[i] != "START" && sentence[i] != "END") {
						Graph G(sentence.size() - 2, i-1); //not including "START" and "END"
						for (int i = 0; i < sentence.size(); i++) { //i = head
							if (sentence[i] != "START" && sentence[i] != "END" ) {
								for (int j = 0; j < sentence.size(); j++) { // j = modifier
									if (sentence[j] != "START" && sentence[j] != "END" && i != j) {
										int dist = j - i;
										//std::string target_word = sentence[j] + NumberToString<int> (dist);
										int arr_temp[6] = {small_pos_tag_map[sentence[j]], small_pos_tag_map[sentence[i]], small_pos_tag_map[sentence[j-1]], small_pos_tag_map[sentence[j+1]], small_pos_tag_map[sentence[i-1]], small_pos_tag_map[sentence[i+1]]};
										//compute normalising constant
										double normalising_constant = 0.0;
										//first, compute the prediction
										for (int iter2 = 0; iter2 < 5; iter2++) {
											prediction += (Ci[iter2].matrix() * R.col(arr_temp[iter2+1]).matrix()); //prediction
										}
										//then, iterate through all possible target words using the prediction value
										for (int iter = 0; iter < Q.cols(); iter++) {
											normalising_constant += exp(prediction.transpose() * Q.col(iter).matrix());//dot product
										}
										if (normalising_constant == 0) {
											std::cout << "Error: normalising_constant should not be 0!!!!" << std::endl;
										} 					
										double log_prob = prediction.transpose() * Q.col(arr_temp[0]).matrix();
										log_prob -= log(normalising_constant);
										G.Connect(i-1, j-1, -log_prob);
										//reset prediction
										prediction.setZero();
									}
								}
							}
						}
						std::vector<Edge> curr_MST = G.get_MST();
						if (get_tree_score(curr_MST) < min_MST) {
							min_MST = get_tree_score(curr_MST);
							if (min_graph != nullptr) {
								delete min_graph;
							} 
							min_graph = new Graph(G);					
						}
					}
				}
				std::vector<Edge> result = min_graph->get_MST();
				for (std::vector<Edge>::iterator it = result.begin(); it != result.end(); it++) {
					final_result << it->getDestination()+1 << " " << it->getSource()+1 << std::endl; 
				}
				final_result << std::endl; 
			}
		}

		sentences.close();
		validation_pairs.close();
		final_result.close();
		final_result.close();
	}
	return 0;
}