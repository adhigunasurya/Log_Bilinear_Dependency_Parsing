#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <random>
#include <stdlib.h>

#ifndef FEATURE_SIZE
#define FEATURE_SIZE 5
#endif

#ifndef NO_OF_POSTAGS
#define NO_OF_POSTAGS 45
#endif

#ifndef MAX_DISTANCE
#define MAX_DISTANCE 49
#endif

#ifndef BATCH_SIZE
#define BATCH_SIZE 100
#endif

#ifndef LEARNING_RATE
#define LEARNING_RATE 0.003
#endif

#ifndef LEARNING_DECAY
#define LEARNING_DECAY 0.00001
#endif

#ifndef NO_OF_EPOCH
#define NO_OF_EPOCH 10
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
	}
	fine_pos_tags.close();
	return result;
}

//returns a mapping between the pos tags, concatenated with all possible distance, with integers
//the difference with function postag_mapping is that this time the mapping is concatenated with all the possible distance (from 1-49), 
//while in the function postag_mapping, the pos tags are "vanilla" and not concatenated with the distance
std::map<std::string, int> get_postag_bias(Eigen::ArrayXd& bias, std::string distance_count_file, std::string postags_file) {
	std::string line;
	std::ifstream distance_count(distance_count_file);
	std::ifstream postags(postags_file);
	std::vector<std::string> temp;
	std::map<std::string, int> result_map;
	if (!distance_count.is_open() || !postags.is_open()) {
		std::cout << "Error: file containing distance counts fail to open" << std::endl;
	} else {
		int total = 0;
		//first, create the mapping between list of pos tags (appended with distance) and their counts, initialised to 0
		int j = -1;
		while (getline(postags, line)) {
			for (int i = 1; i <= MAX_DISTANCE; i++) {
				j++;
				std::ostringstream s; 
				s << i;
				std::string result = s.str();
				result_map[line + result] = j;
			}
		}
		//calculate the total occurences of postags
		while (getline(distance_count, line)) {
			temp = split(line, ' ');
			total += std::stoi(temp[4]);
			temp.clear();
		}
		//std::cout << "the total counts for all distance is: " << total << std::endl;
		distance_count.close();
		//re-open the file and assign the probabilities (count of each pos tag divided by total) to the bias array
		distance_count.open(distance_count_file);
		bias = Eigen::ArrayXd::Zero(result_map.size());
		while (getline(distance_count, line)) {
			temp = split(line, ' ');
			int idx = result_map[temp[2]];
			int pos_tag_count = std::stoi(temp[4]);
			bias(idx) = ((double) pos_tag_count) / ((double) total);
			temp.clear();
		}
		distance_count.close();
		postags.close();
	}
	return result_map;
}

int main()
{
	//initialisation
	Eigen::ArrayXXd R = Eigen::ArrayXXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE); //matrix R, dimension: Nw(number of pos tags conjoined with distance) x Nf
	std::default_random_engine generator;
	std::normal_distribution<double> distr(0.0, 0.0005);
	//initialise the matrix R with small, random, non-zero numbers
	for (int i = 0; i < R.rows(); i++) {
		for (int j = 0; j < R.cols(); j++) {
			double d = distr(generator);
			//avoid having 0's as the weight, try to get all the coefficients initialised with small numbers
			while (d == 0.0) {
				d = distr(generator);
			}
			R(i, j) = d;
		}
	}
	//create an array of Ci (where each Ci is a matrix of Nf x Nf that governs the interaction between the ith feature and vn)
	std::vector<Eigen::ArrayXXd> Ci;
	for (int i = 0; i < 5; i++) {
		Eigen::ArrayXXd Ci_temp = Eigen::ArrayXXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
		std::default_random_engine generator;
		std::normal_distribution<double> distr(0.0, 0.0005);
		for (int j = 0; j < Ci_temp.rows(); j++) {
			for (int k = 0; k < Ci_temp.cols(); k++) {
				double d = distr(generator);
				//avoid having 0's as the weight, try to get all the coefficients initialised with small numbers
				while (d == 0.0) {
					d = distr(generator);
				}
				Ci_temp(j, k) = d;
			}
		}
		Ci.push_back(Ci_temp);
	}
	//initialise vi. We use a context of 5 (the postag of the head itself, the postag before the modifier, the postag after the modifier, the postag before the head, and the postag after the head)
	std::vector<Eigen::ArrayXd> vi;
	for (int i = 0; i < 5; i++) {
		Eigen::ArrayXd temp = Eigen::ArrayXd::Zero(NO_OF_POSTAGS);
		vi.push_back(temp);
	}
	//initialise vn. Both vi and vn are column vectors with all 0's and a single 1, which indicates the POS-tag we desire to get
	Eigen::ArrayXd vn = Eigen::ArrayXd::Zero(NO_OF_POSTAGS);
	//FOR NOW, we ignore the term br^transpose * R^transpose*vn, but find out what it is
	//initialise the bias term bv
	Eigen::ArrayXd bv;
	/* print the array of biases into an external file
	std::ofstream final_postag_bias("final_postag_bias.txt"); 
	if (!final_postag_bias.is_open()) {
		std::cout << "fail to create final postag bias file" << std::endl;
	} else {
		std::map<std::string, int> big_pos_tag_map = get_postag_bias(bv, "distance_count.txt", "fine_pos_tags.txt");
		for (std::map<std::string, int>::iterator it = big_pos_tag_map.begin(); it != big_pos_tag_map.end(); it++) {
			final_postag_bias << "Pos tag " << it->first << " has bias value " << bv(it->second) << std::endl;
		}
	} */
	//std::cout << "one example: " << bv(1) << std::endl;
	//main bit
	std::map<std::string, int> mapping = postag_mapping("fine_pos_tags.txt");
	std::ifstream pairs ("real_pairs_small.txt");
	//initialise the gradient matrix
	std::vector<Eigen::MatrixXd> dCi;
	for (int i = 0; i < 5; i++) {
		Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
		dCi.push_back(temp);
	}
	Eigen::MatrixXd dR = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE);
	int no_of_iter = 0;
	std::ofstream parameters ("parameters_each_epoch_small_trial_1.txt");
	if (!parameters.is_open()) {
		std::cout << "Error: output file not opening" << std::endl;
	}
	for (int i = 0; i < NO_OF_EPOCH; i++) {
		if (!pairs.is_open()) {
			std::cout << "Fail to open input file" << std::endl;
		} else {
			std::vector<std::string> temp;
			std::string line;
			while (getline(pairs, line)) {
				if (line.size() > 1) {
					no_of_iter++;
					temp = split(line, ' ');
					//initialise an array with the elements as the indices of: postag of the modifier, postag of the head, postag of before modifier, postag of after modifier, postag of before head, postag of after head
					int arr_temp[6] = {mapping[temp[0]], mapping[temp[1]], mapping[temp[2]], mapping[temp[3]], mapping[temp[4]], mapping[temp[5]]};
					for (int j = 0; j < 5; j++) {
						(vi[j])[arr_temp[j+1]] = 1.0;
					}
					vn[arr_temp[0]] = 1.0; 
					std::cout << no_of_iter << std::endl;
					//all codes here (calculate the gradients)
					Eigen::MatrixXd tempC = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
					Eigen::ArrayXd tempVn = Eigen::ArrayXd::Zero(NO_OF_POSTAGS);
					//calculate all the gradient changes for each of the Ci, from dCi[0] to dCi[4]
					for (int j = 0; j < 5; j++) {
						for (int k = 0; k < NO_OF_POSTAGS; k++) {
							tempVn(k) = 1.0;
							tempC += (R.matrix().transpose() * vi[j].matrix() * tempVn.matrix().transpose() * R.matrix());
							//reset tempVn without re-allocating
							tempVn(k) = 0.0;
						}
						tempC /= NO_OF_POSTAGS;
						dCi[j] += ((R.matrix().transpose() * vi[j].matrix() * vn.matrix().transpose() * R.matrix()) - tempC);
						//reset tempC without re-allocating
						for (int m = 0; m < tempC.rows(); m++) {
							for (int n = 0; n < tempC.cols(); n++) {
								tempC(m, n) = 0.0;
							}
						}
						//end of reset tempC
					}

					//calculate dR 
					Eigen::MatrixXd tempR_D = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE);
					for (int j = 0; j < 5; j++) {
						tempR_D += (vn.matrix() * vi[j].matrix().transpose() * R.matrix() * Ci[j].matrix()) + (vi[j].matrix() * vn.matrix().transpose() * R.matrix() *  Ci[j].matrix().transpose());
					}
					//tempR_D += (vn.matrix() * vi_5.matrix().transpose() * R.matrix() * Ci[4].matrix()) + (vi_5.matrix() * vn.matrix().transpose() * R.matrix() *  Ci[4].matrix().transpose());

					//reinitialise tempVn
					for (int j = 0; j < tempVn.rows(); j++) {
						tempVn(j) = 0.0;
					}

					Eigen::MatrixXd tempR_M = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE);
					for (int j = 0; j < NO_OF_POSTAGS; j++) {
						tempVn(j) = 1.0;
						for (int k = 0; k < 5; k++) {
							tempR_M += (tempVn.matrix() * vi[k].matrix().transpose() * R.matrix() * Ci[k].matrix()) + (vi[k].matrix() * tempVn.matrix().transpose() * R.matrix() *  Ci[k].matrix().transpose());
						}
						//tempR_M += (tempVn.matrix() * vi_5.matrix().transpose() * R.matrix() * Ci[4].matrix()) + (vi_5.matrix() * tempVn.matrix().transpose() * R.matrix() *  Ci[4].matrix().transpose());
						tempVn(j) = 0.0;
					}
					tempR_M /= NO_OF_POSTAGS;
					dR += (tempR_D - tempR_M);
					tempR_D.resize(0, 0);
					tempR_M.resize(0, 0);
					//reset the vectors vi_1 to vi_5 and vn, only 1 element in the vector may be "1" and the rest has to be "0"
					for (int j = 0; j < 5; j++) {
						(vi[j])[arr_temp[j+1]] = 0.0;
					}
					/*
					vi_1[arr_temp[1]] = 0.0;
					vi_2[arr_temp[2]] = 0.0; 
					vi_3[arr_temp[3]] = 0.0; 
					vi_4[arr_temp[4]] = 0.0; 
					vi_5[arr_temp[5]] = 0.0; */
					//reset the vector vn
					vn[arr_temp[0]] = 0.0; 

					//update the gradient if no_of_iter mod mini batch size == 0, reset the gradient matrix to 0
					if (no_of_iter % BATCH_SIZE == 0) {
						double curr_learning_rate = LEARNING_RATE / (1 + no_of_iter * LEARNING_DECAY);
						for (int j = 0; j < 5; j++) {
							Ci[j] = Ci[j].matrix() + curr_learning_rate * dCi[j];
							//reset the matrix dCi[j]
							for (int k = 0; k < dCi[j].rows(); k++) {
								for (int m = 0; m < dCi[j].cols(); m++) {
									(dCi[j])(k, m) = 0.0;
								}
							}
						}
						//Ci[4] = Ci[4].matrix() + curr_learning_rate * dC5;
						R = R.matrix() + curr_learning_rate * dR;

						//reset the matrices dR
						for (int j = 0; j < dR.rows(); j++) {
							for (int k = 0; k < dR.cols(); k++) {
								dR(j, k) = 0.0;
							}
						}

						/*dC1.resize(0, 0); 
						dC2.resize(0, 0);
						dC3.resize(0, 0);
						dC4.resize(0, 0);
						dC5.resize(0, 0); 
						dR.resize(0, 0);*/ 

/*
						dC1 = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
						dC2 = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
						dC3 = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
						dC4 = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
						dC5 = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
						dR = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE); */
					}
				}
				temp.clear();
			}
			pairs.close();
		}
		parameters << "Epoch number: " << i+1 << std::endl;
		for (int j = 0; j < 5; j++) {
			parameters << "Matrix C" << j+1 << std::endl; 
			parameters << Ci[j] << std::endl;
			parameters << std::endl;
			parameters << "End of matrix C" << j+1 << std::endl;
		}
		parameters << "Matrix R" << std::endl;
		parameters << R << std::endl;
		parameters << "End of matrix R" << std::endl;
		pairs.open("real_pairs_small.txt");
	}
	pairs.close();
	parameters.close();
	return 0;
}