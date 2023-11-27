#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <time.h>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <iterator>
using namespace std;
#include "feedbacktrellis.h"
#include "feedforwardtrellis.h"
#include "dualtrellis.h"
#include "listdecoder.h"
#include "lowratelistdecoder.h"
#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"
#include "helper_functions.h"

#include "MinHeap.h"
//static default_random_engine generator;
static std::random_device rd {};
static std::mt19937 generator {rd()};

struct codeInformation{
	int k;
	int n;
	int v;
	int crcDeg;
	int crc;
	int numInfoBits;
	std::vector<int> numerators;
	int denominator;
	std::vector<std::vector<int>> hMatrix;
};

std::vector<double> addNoise(std::vector<int> encodedMsg, double SNR);
void crc_calculation(std::vector<int>& input_data, int crc_bits_num, int crc_dec);

void generateMessages(codeInformation code);

void validateDecoder(codeInformation code);
void testDebugging(codeInformation code);
void speedTest(codeInformation code);
void performanceConsistencyEvaluation(codeInformation code);
void conferencePaperFigure(codeInformation code);
void listSizeHistograms(codeInformation code);
void wavaFERPerformance(codeInformation code);
void ztTrellisFERPerformance(codeInformation code);
void tbTrellisFERPerformance(codeInformation code);

void test(codeInformation code);
void testFERs(codeInformation code);
void pointTest();

void lowrate_test(codeInformation code);
void lowrate_test_variable_listsize(codeInformation code);
void lowrate_test_zero_noise_codewords(codeInformation code);
void lowrate_list_size_for_fixed_dist(codeInformation code);
void lowrate_test1(codeInformation code);
void distance_distribution(codeInformation code);
void q_func_validation(codeInformation code);
void lowrate_find_reciprocity(codeInformation code);

void example_function(codeInformation code);

void tb_nontb_list_sizes(codeInformation code);

void generateNeighborList(codeInformation code);

void linearityDecoding(codeInformation code);
void distanceThreshold_midpoint(codeInformation code);
void minDistPastVoronoi(codeInformation code);
void farthestCwd(codeInformation code);
void generateMiniBallNeighbors(codeInformation code);
void miniBallDecoder(codeInformation code);
void miniBallThresholdTest(codeInformation code);

void ISTC_sim(codeInformation code);

std::vector<int> generateRandomCRCMessage(codeInformation code);
std::vector<double> generateTransmittedMessage(std::vector<int> originalMessage, FeedforwardTrellis encodingTrellis, double snr, std::vector<int> puncturedIndices, bool noiseless = false);

// int main(){
// 	codeInformation code;
// 	// ///	ZTCC CODE INFORMATION
// 	// code.k = 3;						// numerator of the rate
// 	// code.n = 4;						// denominator of the rate
// 	// code.v = 4;						// number of memory elements
// 	// code.crcDeg = 4;				// degree of the crc (one more than the number of crc bits)
// 	// code.crc = 9;					// crc in decimal, convert from hex if needed
// 	// code.numInfoBits = 93;			// number of information bits in the message, 
// 	// 								// actual info length of ZT is code.numInfoBits - 6 terminating bits 

// 	// // optimal code numerators and denominator are known, and are given in octal
// 	// code.numerators = {33, 25, 37};
// 	// code.denominator = 31;

// 	// // hMatrix is {{denom},{numerators[last]},...,{numerators[0]}}, converted to binary then flipped lr
// 	// code.hMatrix = {{1,0,0,1,1}, {1,1,1,1,1}, {1,0,1,0,1}, {1,1,0,1,1}};

// 	code.k = 3;						// numerator of the rate
// 	code.n = 4;						// denominator of the rate
// 	code.v = 6;						// number of memory elements
// 	code.crcDeg = 11;				// degree of the crc (one more than the number of crc bits)
// 	code.crc = 1321;					// crc in decimal, convert from hex if needed
// 	code.numInfoBits = 86;			// number of information bits in the message, 
									
// 	// optimal code numerators and denominator are known, and are given in octal
// 	code.numerators = {107, 135, 133};
// 	code.denominator = 141;

// 	// hMatrix is {{denom},{numerators[last]},...,{numerators[0]}}, converted to binary then flipped lr
// 	code.hMatrix = {{1,0,0,0,0,1,1}, {1,1,0,1,1,0,1}, {1,0,1,1,1,0,1}, {1,1,1,0,0,0,1}};


// 	/*
// 	code.k = 3;				// numerator of the rate
// 	code.n = 4;				// denominator of the rate
// 	code.v = 5;				// number of memory elements
// 	code.crcDeg = 7;			// degree of the crc
// 	code.crc = 83;			// crc is given in decimal, convert from hex if necessary
// 	code.numInfoBits = 42;	// number of information bits
// 	code.numerators = { 47, 73, 57 };
// 	code.denominator = 75;
// 	code.hMatrix = { {1, 0, 1, 1, 1, 1}, {1, 1, 1, 1, 0, 1}, {1,1,0,1,1,1}, {1,1,1,0,0,1} };
// 	*/
// 	// test(code);
// 	tbTrellisFERPerformance(code);
// 	// test(code);
// 	return 0;
// }

int main(){

	std::cout << __cplusplus << std::endl;
	codeInformation code;
	
	code.k = 1;						// numerator of the rate
	code.n = 2;						// denominator of the rate
	code.v = 8;						// number of memory elements
    code.crcDeg = 8;				// degree of the crc (m + 1)
	code.crc = 0xff;				// crc in decimal, convert from hex if needed
	code.numInfoBits = 64;		
	// code.crcDeg = 1;			
	// code.crc = 0x1;				// crc in decimal, convert from hex if needed
	// code.numInfoBits = 71;	
	
	// // from hengjie's journal paper	
	// code.v = 6;		
	// code.crcDeg = 4;			
	// code.crc = 0x1B;				// crc in decimal, convert from hex if needed
	// code.numInfoBits = 64;	// N = 67*2 = 134
    							
	// optimal code numerators and denominator are known, and are given in octal
	code.numerators = {561, 753};
	// code.numerators = {133, 171};

	tbTrellisFERPerformance(code);
	
	// generateNeighborList(code);
	// linearityDecoding(code);
	
	// generateMiniBallNeighbors(code);
	// miniBallDecoder(code);
	// miniBallThresholdTest(code);

	// distanceThreshold_midpoint(code);
	// minDistPastVoronoi(code);
	// farthestCwd(code);
	
	return 0;
}

void ISTC_sim(codeInformation code){
	std::cout << "running the ISTC sim" << std::endl;

	// set random seed for message generation
	srand(time(NULL));

	// check to make sure the code has valid values
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	int listSize = 1e5; 
	std::vector<double> EbN0 = {1, 1.5};
	vector<double> SNR; 										// SNR is required for noise computations
	double offset = 10 * log10((double)2*64 / (2*128)); 		// real rate of this code *2 in the log
	for (int i=0; i< EbN0.size(); i++)
		SNR.push_back(EbN0[i] + offset);
	int maxNumTrials = 1e4;								// the max trials at each point, the necessary actual value will likely be greater

	std::vector<int> puncturedIndices = {4, 10, 21, 24, 31, 37, 42, 48, 59, 62, 69, 75, 80, 86, 97, 100, 107, 113, 118, 124, 135, 138, 145, 151};
	// the below are the relevant initializations for low rate
	// decoding, this will be altered for high rate or ZTCC, for example
	FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
	LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg, code.crc);

	// running the simulations. in this example, we are simulation to collect TFR data
	for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++){
		double snr = SNR[snrIndex];
		double standardMeanListSize = 0;
		int standardNumErrors = 0;
		int standardListSizeExceeded = 0;

		for(int numTrials = 0; numTrials < maxNumTrials; numTrials++){
			if (numTrials % 1000 == 0){
				std::cout << "currently at " << numTrials << std::endl;
			}
			// the following lines generate a message that we can decode
			std::vector<int> originalMessage = generateRandomCRCMessage(code);
			std::vector<double> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices);

			LowRateListDecoder::messageInformation standardDecoding = listDecoder.lowRateDecoding(transmittedMessage);

			if(standardDecoding.listSizeExceeded){
				standardListSizeExceeded++;
			}
			else if (standardDecoding.message != originalMessage){
				standardNumErrors ++;
				standardMeanListSize += (double)standardDecoding.listSize;
			}
			
		}
		// outputting results at each SNR / EbN0 point. in this case, since we're working with TFR, the ratio
		// of incorrect decodings to total decodings is used.
		std::cout << "at Eb/N0 = " << EbN0[snrIndex] << std::endl;
		std::cout << "mean list size: " << (double)standardMeanListSize/(maxNumTrials - standardListSizeExceeded) << std::endl;
		std::cout << "number of times list size exceeded: " << standardListSizeExceeded << std::endl;
		std::cout << "UER: " << (double)standardNumErrors/(maxNumTrials - standardListSizeExceeded) << std::endl;
		std::cout << "TFR: " << (double)(standardNumErrors + standardListSizeExceeded)/maxNumTrials << std::endl;


		// it may be the case that we want to write the outputs to a file, whether it be for further processing
		// or convenience. the following code snipped illustrates how this would be done. typically, the opening
		// and closing of the file, and declaration of the file name would be done once each outside the loop,
		// but they are included here to keep things organized.
		/*
		std::ofstream outputFile;
		filename = "distance_spectra.txt";
		outputFile.open(filename, fstream::app);
		outputFile << "example information to be saved";
		outputFile.close();
		*/
	}
	std::cout << "concluded simulation" << std::endl;

}

// void miniBallThresholdTest(codeInformation code){
// std::cout << "recording mini ball metrics" << std::endl;

// 	// set random seed for message generation
// 	srand(time(NULL));

// 	// check to make sure the code has valid values
// 	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
// 		std::cout << "invalid msg + crc length" << std::endl;
// 		return;
// 	}

// 	int listSize = 1e5;  // for ML decoding result
// 	int maxNumTrials = 2;
// 	// std::vector<double> EbN0 = {1, 1.5, 2, 2.5, 3, 3.5};
// 	std::vector<double> EbN0 = {1};
// 	vector<double> SNR; 										// SNR is required for noise computations
// 	double offset = 10 * log10((double)2*code.numInfoBits / (2*(code.numInfoBits+code.crcDeg-1))); 		// real rate of this code *2 in the log
// 	for (int i=0; i< EbN0.size(); i++)
// 		SNR.push_back(EbN0[i] + offset);

// 	std::vector<int> puncturedIndices = {};
// 	// the below are the relevant initializations for low rate
// 	// decoding, this will be altered for high rate or ZTCC, for example
// 	FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
// 	LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg, code.crc);

// 	std::string dir = "k" + std::to_string(code.k) + "n" + std::to_string(code.n) + "v" + std::to_string(code.v) +
// 					 "m" + std::to_string(code.crcDeg - 1) + "K" + std::to_string(code.numInfoBits) + "/";
// 	std::cout << dir << std::endl;
// 	listDecoder.readMiniBallNeighborList(dir,1);
// 	std::cout << "neighbors read in, beginning trials" << std::endl;

// 	// std::ofstream outputFile_decoded_correct; // distance between the received word and the decoded codeword
// 	// string filename_correct = "miniBall_L2048_correctdecoding_EbNo1.txt";
// 	// outputFile_decoded_correct.open(filename_correct);

// 	// std::ofstream outputFile_decoded_incorrect; // distance between the received word and the decoded codeword
// 	// string filename_incorrect = "miniBall_L2048_incorrectdecoding_EbNo1.txt";
// 	// outputFile_decoded_incorrect.open(filename_incorrect);

// 	// running the simulations. in this example, we are simulation to collect TFR data
// 	for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++){
// 		double snr = SNR[snrIndex];
// 		double MeanListSize = 0;
// 		int NumErrors = 0;
// 		int ListSizeExceeded = 0;
// 		int maxListSize = 0;

// 		double standardMeanListSize = 0;
// 		int standardNumErrors = 0;
// 		int standardListSizeExceeded = 0;
// 		int standardMaxListSize = 0;

// 		// apply SLVD on the first codeword of each ESD
// 		for(int numTrials = 0; numTrials < maxNumTrials; numTrials++){
// 			// the following lines generate a message that we can decode
// 			std::vector<int> originalMessage = generateRandomCRCMessage(code);
// 			// encode the message
// 			std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);
// 			// std::cout << "original codeword: " << std::endl;
// 			// print_int_vector(encodedMessage);
// 			// add noise
// 			std::vector<double> transmittedMessage = addNoise(encodedMessage, snr);

// 			// std::vector<double> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices);
// 			// project onto sphere with squared distance n
// 			double cur_energy = 0;
// 			for(int i = 0; i < transmittedMessage.size(); i++){
// 				cur_energy += std::pow(transmittedMessage[i], 2);
// 			}
// 			double scale = double(transmittedMessage.size())/cur_energy;

// 			for(int i = 0; i < transmittedMessage.size(); i++){
// 				transmittedMessage[i] = transmittedMessage[i] * sqrt(scale);
// 			}
// 			// print_double_vector(transmittedMessage);

// 			LowRateListDecoder::messageInformation standardDecoding = listDecoder.lowRateDecoding(transmittedMessage);
// 			// if(!standardDecoding.listSizeExceeded){
// 			// 	std::cout << "SLVD decoded! metric: " << standardDecoding.metric << std::endl;
// 			// }
// 			// std::cout << "SLVD codeword: " << std::endl;
// 			// print_int_vector(standardDecoding.codeword);
			
// 			LowRateListDecoder::messageInformation miniBallDecoding = listDecoder.miniBallDecoder(transmittedMessage, INT_MAX);
// 			// std::cout << "miniball codeword metric: " << miniBallDecoding.metric << std::endl;
// 			// print_int_vector(miniBallDecoding.codeword);

// 			listDecoder.miniBallMetrics(encodedMessage, transmittedMessage, standardDecoding.codeword);

// 			// // mini-ball decoder
// 			// if(miniBallDecoding.listSizeExceeded){
// 			// 	ListSizeExceeded++;
// 			// }
// 			// else{
// 			// 	MeanListSize += (double)miniBallDecoding.listSize;
// 			// 	if (miniBallDecoding.message != originalMessage){
// 			// 		NumErrors ++;
// 			// 	}
// 			// 	if (miniBallDecoding.listSize > maxListSize){
// 			// 		maxListSize = miniBallDecoding.listSize;
// 			// 	}
// 			// }

// 			// // SLVD
// 			// if(standardDecoding.listSizeExceeded){
// 			// 	standardListSizeExceeded++;
// 			// }
// 			// else{
// 			// 	standardMeanListSize += (double)standardDecoding.listSize;
// 			// 	if (standardDecoding.message != originalMessage){
// 			// 		standardNumErrors ++;
// 			// 	}
// 			// 	if (standardDecoding.listSize > standardMaxListSize){
// 			// 		standardMaxListSize = standardDecoding.listSize;
// 			// 	}
// 			// }

// 			// // if the decoder made the ML decision regardless of whether it was correct
// 			// if(miniBallDecoding.message != multiTrellisDecoding.message){ // not ML decision
// 			// 	// distance between the received word and the decoded codeword
// 			// 	outputFile_decoded_incorrect << miniBallDecoding.metric << std::endl;
// 			// }
// 			// else{ //correct decoding
// 			// 	outputFile_decoded_correct<< miniBallDecoding.metric << std::endl;
// 			// }
			
// 		}
// 		// outputting results at each SNR / EbN0 point. in this case, since we're working with TFR, the ratio
// 		// of incorrect decodings to total decodings is used.
// 		// std::cout << "at Eb/N0 = " << EbN0[snrIndex] << std::endl;
		
// 		// std::cout << "mini ball linearity decoder" << std::endl;
// 		// std::cout << "mean list size: " << (double)MeanListSize/(maxNumTrials - ListSizeExceeded) << std::endl;
// 		// std::cout << "max list size: " << (double)maxListSize << std::endl;
// 		// std::cout << "number of times list size exceeded: " << ListSizeExceeded << std::endl;
// 		// std::cout << "UER: " << (double)NumErrors/(maxNumTrials - ListSizeExceeded) << std::endl;
// 		// std::cout << "TFR: " << (double)(NumErrors + ListSizeExceeded)/maxNumTrials << std::endl;

// 		// std::cout << "standard SLVD decoder" << std::endl;
// 		// std::cout << "mean list size: " << (double)standardMeanListSize/(maxNumTrials - standardListSizeExceeded) << std::endl;
// 		// std::cout << "max list size: " << (double)standardMaxListSize << std::endl;
// 		// std::cout << "number of times list size exceeded: " << standardListSizeExceeded << std::endl;
// 		// std::cout << "UER: " << (double)standardNumErrors/(maxNumTrials - standardListSizeExceeded) << std::endl;
// 		// std::cout << "TFR: " << (double)(standardNumErrors + standardListSizeExceeded)/maxNumTrials << std::endl;
// 	}
// 	std::cout << "concluded simulation" << std::endl;
// 	// outputFile_decoded_correct.close();
// 	// outputFile_decoded_incorrect.close();	

// }

// void generateMiniBallNeighbors(codeInformation code){
// std::cout << "finding miniball neighbors" << std::endl;

// 	// set random seed for message generation
// 	srand(time(NULL));

// 	// check to make sure the code has valid values
// 	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
// 		std::cout << "invalid msg + crc length" << std::endl;
// 		return;
// 	}

// 	int listSize = 1; 
// 	int num_neighbors = 3; // neighbors with the closest 2 distances

// 	// the below are the relevant initializations for low rate
// 	// decoding, this will be altered for high rate or ZTCC, for example
// 	FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
// 	LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg, code.crc);

// 	std::string dir = "k" + std::to_string(code.k) + "n" + std::to_string(code.n) + "v" + std::to_string(code.v) +
// 					 "m" + std::to_string(code.crcDeg - 1) + "K" + std::to_string(code.numInfoBits) + "/";
// 	std::cout << dir << std::endl;
// 	listDecoder.readNeighborList(dir);
// 	std::cout << "neighbors read in, beginning trials" << std::endl;

// 	// apply SLVD on the first codeword of each ESD, starting at 1
// 	for(int esd = 1; esd < listDecoder.codewordNeighborList.size(); esd++){
// 		std::vector<int> curNeighbor = listDecoder.codewordNeighborList[esd][0];
// 		std::vector<double> transmittedMessage;
// 		for (int i=0; i<curNeighbor.size(); i++){
// 			transmittedMessage.push_back((double)curNeighbor[i]);
// 		}
// 		listDecoder.findMiniBallNeighbors(transmittedMessage,listDecoder.messageNeighborList[esd][0], dir, esd, num_neighbors);
// 	}
// 	std::cout << "concluded simulation" << std::endl;
// }


void miniBallDecoder(codeInformation code){
std::cout << "running the mini multi ball decoder (MMD)" << std::endl;

	// set random seed for message generation
	srand(time(NULL));

	// check to make sure the code has valid values
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}
	// std::vector<double> thresholds = {8.313}; // radius around single cwd for N=134 code
	// std::vector<double> thresholds = {5,8.72,10}; // radius around single cwd for N=134 code
	std::vector<double> thresholds = {5}; // radius around single cwd for N=134 code
	for (int t=0; t< thresholds.size(); t++){
		double threshold = thresholds[t];
		std::cout << "current threshold is: " << threshold << std::endl;

		int listSize = 2048;  // for ML decoding result
		int num_neighbors = 1;
		int maxNumTrials = 1e4;
		std::vector<double> EbN0 = {2, 2.5};
		// std::vector<double> EbN0 = {3};
		vector<double> SNR; 										// SNR is required for noise computations
		double offset = 10 * log10((double)2*code.numInfoBits / (2*(code.numInfoBits+code.crcDeg-1))); 		// real rate of this code *2 in the log
		for (int i=0; i< EbN0.size(); i++)
			SNR.push_back(EbN0[i] + offset);

		std::vector<int> puncturedIndices = {};
		// the below are the relevant initializations for low rate
		// decoding, this will be altered for high rate or ZTCC, for example
		FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
		LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg, code.crc);

		
		std::string dir = "k" + std::to_string(code.k) + "n" + std::to_string(code.n) + "v" + std::to_string(code.v) +
						"m" + std::to_string(code.crcDeg - 1) + "K" + std::to_string(code.numInfoBits) + "/";
		std::cout << dir << std::endl;
		listDecoder.readMiniBallNeighborList(dir,num_neighbors);
		std::cout << "neighbors read in, beginning trials" << std::endl;
		

		// std::ofstream outputFile_decoded_correct; // distance between the received word and the decoded codeword
		// string filename_correct = "miniBall_L2048_correctdecoding_EbNo3.txt";
		// outputFile_decoded_correct.open(filename_correct);

		// std::ofstream outputFile_decoded_incorrect; // distance between the received word and the decoded codeword
		// string filename_incorrect = "miniBall_L2048_incorrectdecoding_EbNo3.txt";
		// outputFile_decoded_incorrect.open(filename_incorrect);

		// running the simulations. in this example, we are simulation to collect TFR data
		for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++){
			double snr = SNR[snrIndex];
			double MeanListSize = 0;
			int NumErrors = 0;
			int ListSizeExceeded = 0;
			int maxListSize = 0;

			// double standardMeanListSize = 0;
			// int standardNumErrors = 0;
			// int standardListSizeExceeded = 0;
			// int standardMaxListSize = 0;

			// apply SLVD on the first codeword of each ESD
			for(int numTrials = 0; numTrials < maxNumTrials; numTrials++){
				// the following lines generate a message that we can decode
				std::vector<int> originalMessage = generateRandomCRCMessage(code);
				std::vector<double> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices);
				// project onto sphere with squared distance n
				double cur_energy = 0;
				for(int i = 0; i < transmittedMessage.size(); i++){
					cur_energy += std::pow(transmittedMessage[i], 2);
				}
				double scale = double(transmittedMessage.size())/cur_energy;

				for(int i = 0; i < transmittedMessage.size(); i++){
					transmittedMessage[i] = transmittedMessage[i] * sqrt(scale);
				}
				// std::cout<< "original msg" << std::endl;
				// print_int_vector(originalMessage);
				// std::cout<< "transmitted msg" << std::endl;
				// print_double_vector(transmittedMessage);

				LowRateListDecoder::messageInformation miniBallDecoding = listDecoder.miniBallDecoder(transmittedMessage, threshold);
				// LowRateListDecoder::messageInformation multiTrellisDecoding = listDecoder.multiTrellisLowRateDecoding(transmittedMessage);
				// LowRateListDecoder::messageInformation standardDecoding = listDecoder.lowRateDecoding(transmittedMessage);

				// mini-ball decoder
				if(miniBallDecoding.listSizeExceeded){
					ListSizeExceeded++;
				}
				else{
					MeanListSize += (double)miniBallDecoding.listSize;
					if (miniBallDecoding.message != originalMessage){
						NumErrors ++;
					}
					if (miniBallDecoding.listSize > maxListSize){
						maxListSize = miniBallDecoding.listSize;
					}
				}

				// // SLVD
				// if(standardDecoding.listSizeExceeded){
				// 	standardListSizeExceeded++;
				// }
				// else{
				// 	standardMeanListSize += (double)standardDecoding.listSize;
				// 	if (standardDecoding.message != originalMessage){
				// 		standardNumErrors ++;
				// 	}
				// 	if (standardDecoding.listSize > standardMaxListSize){
				// 		standardMaxListSize = standardDecoding.listSize;
				// 	}
				// }

				// // if the decoder made the ML decision regardless of whether it was correct
				// if(miniBallDecoding.message != standardDecoding.message){ // not ML decision
				// 	// distance between the received word and the decoded codeword
				// 	outputFile_decoded_incorrect << miniBallDecoding.metric << std::endl;
				// }
				// else{ //correct decoding
				// 	outputFile_decoded_correct<< miniBallDecoding.metric << std::endl;
				// }
				
			}
			// outputting results at each SNR / EbN0 point. in this case, since we're working with TFR, the ratio
			// of incorrect decodings to total decodings is used.
			std::cout << "at Eb/N0 = " << EbN0[snrIndex] << std::endl;
			
			std::cout << "mini ball linearity decoder" << std::endl;
			std::cout << "with number of nearest distances: " << num_neighbors << std::endl;
			std::cout << "mean list size: " << (double)MeanListSize/(maxNumTrials - ListSizeExceeded) << std::endl;
			std::cout << "max list size: " << (double)maxListSize << std::endl;
			std::cout << "number of times list size exceeded: " << ListSizeExceeded << std::endl;
			std::cout << "UER: " << (double)NumErrors/(maxNumTrials - ListSizeExceeded) << std::endl;
			std::cout << "TFR: " << (double)(NumErrors + ListSizeExceeded)/maxNumTrials << std::endl;

			// std::cout << "standard SLVD decoder" << std::endl;
			// std::cout << "mean list size: " << (double)standardMeanListSize/(maxNumTrials - standardListSizeExceeded) << std::endl;
			// std::cout << "max list size: " << (double)standardMaxListSize << std::endl;
			// std::cout << "number of times list size exceeded: " << standardListSizeExceeded << std::endl;
			// std::cout << "total number of errors: " << standardNumErrors + standardListSizeExceeded << std::endl;
			// std::cout << "UER: " << (double)standardNumErrors/(maxNumTrials - standardListSizeExceeded) << std::endl;
			// std::cout << "TFR: " << (double)(standardNumErrors + standardListSizeExceeded)/maxNumTrials << std::endl;
		}
	}
	std::cout << "concluded simulation" << std::endl;
	// outputFile_decoded_correct.close();
	// outputFile_decoded_incorrect.close();	

}

void farthestCwd(codeInformation code){
	std::cout << "running the standard decoder for farthest codeword" << std::endl;

	// set random seed for message generation
	srand(time(NULL));

	// check to make sure the code has valid values
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}
	
	int listSize = 2048;
	
	std::vector<int> puncturedIndices = {};
	// the below are the relevant initializations for low rate
	// decoding, this will be altered for high rate or ZTCC, for example
	FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
	LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg, code.crc);

	// all-zeros message
	std::vector<int> originalMessage;
	for(int i = 0; i < code.numInfoBits; i++)
		originalMessage.push_back(0);
	// compute the CRC
	crc_calculation(originalMessage, code.crcDeg, code.crc);
	// encode the message
	std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);
	// test to see if we find the farthest codeword
	std::vector<double> transmittedMsg;
	for(int i = 0; i < encodedMessage.size(); i++){
		transmittedMsg.push_back(double(encodedMessage[i]));
	}
	
	// LowRateListDecoder::messageInformation standardDecoding = listDecoder.minimumLikelihoodLowRateDecoding(transmittedMsg);
	
	
	// read in farthest codewords
	std::ifstream inputFile;
    inputFile.open("farthest_codewords.txt");
	std::string codewordLine;
	std::vector<std::vector<int>> codewordList;
	// reads the file in line by line
	while(std::getline(inputFile, codewordLine)){
		std::vector<int> codeword;
		std::stringstream stream(codewordLine);
		std::string val;
		while(stream >> val)
			codeword.push_back(std::stoi(val));
		codewordList.push_back(codeword);
	}

	inputFile.close();
	for (int i=0; i<codewordList.size(); i++){


		bool correct_flag = false;
		int counter = 0;
		// the farthest TB codeword is the all-ones codeword in this case
		std::vector<int> currCodeword = codewordList[i];
		
		// loop until we reach the first cwd whose midpoint decodes back to all-zeros cwd
		while(!correct_flag){	
			// std::cout << "current codeword is: " << std::endl;
			// print_int_vector(currCodeword);
			// find midpoint codeword
			std::vector<double> justPastMidpointCodeword;
			for(int i = 0; i < encodedMessage.size(); i++){
				justPastMidpointCodeword.push_back(0.499*double(currCodeword[i]) + 0.501*double(encodedMessage[i]));
			}

			// double cur_energy = 0;
			// for(int i = 0; i < justPastMidpointCodeword.size(); i++){
			// 	cur_energy += std::pow(justPastMidpointCodeword[i], 2);
			// }
			// double scale = double(justPastMidpointCodeword.size())/cur_energy;

			// for(int i = 0; i < justPastMidpointCodeword.size(); i++){
			// 	justPastMidpointCodeword[i] = justPastMidpointCodeword[i] * sqrt(scale);
			// }
			
			LowRateListDecoder::messageInformation standardDecoding = listDecoder.lowRateDecoding(justPastMidpointCodeword);
			counter ++;
			// std::cout << "decoded codeword is: " << std::endl;
			// print_int_vector(standardDecoding.codeword);
			// midpoint distance from the all-zeros cwd
			double midMetric = 0;
			double metric = 0;
			for(int i = 0; i < justPastMidpointCodeword.size(); i++){
				midMetric += std::pow(justPastMidpointCodeword[i] - encodedMessage[i], 2);
				metric += std::pow(currCodeword[i] - encodedMessage[i], 2);
			}
			std::cout << "midpoint metric: " << midMetric <<  ", codeword metric: " << metric << std::endl;
			std::cout << "this is the " << counter << "th iteration" << std::endl;

			// if there is decoding error
			if(standardDecoding.message == originalMessage){
				correct_flag = true;
			}
			else{
				currCodeword = standardDecoding.codeword;
			}
		}
	}

	std::cout << "concluded simulation" << std::endl;

}

void minDistPastVoronoi(codeInformation code){
	std::cout << "running the linearity decoder" << std::endl;

	// set random seed for message generation
	srand(time(NULL));

	// check to make sure the code has valid values
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	int listSize = 2048;
	
	double EbN0 = 1;
	double offset = 10 * log10((double)2*64 / (2*71)); 		// real rate of this code *2 in the log
	double snr = EbN0 + offset;								// SNR is required for noise computation

	std::vector<int> puncturedIndices = {};
	// the below are the relevant initializations for low rate
	// decoding, this will be altered for high rate or ZTCC, for example
	FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
	// LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg, code.crc);

	// all-zeros message
	std::vector<int> originalMessage;
	for(int i = 0; i < code.numInfoBits; i++)
		originalMessage.push_back(0);

	// compute the CRC
	crc_calculation(originalMessage, code.crcDeg, code.crc);
	// encode the message
	std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);
	// add noise
	std::vector<double> finalMessage = addNoise(encodedMessage, snr);

	// apply arbitary noise, evaluate energy, scale to dist = 90
	double setDist = 30;
	double cur_dist = 0;
	for(int i = 0; i < finalMessage.size(); i++){
		cur_dist += std::pow(finalMessage[i] - encodedMessage[i], 2);
	}
	double scale =  setDist/cur_dist;

	std::vector<double> transmittedMessage;
	for(int i = 0; i < finalMessage.size(); i++){
		transmittedMessage.push_back(finalMessage[i] * sqrt(scale) - sqrt(scale) + 1);
	}
	// note: transmitted msg at this point is of distance 90 from the all-zeros cwd, 
	// but it's not necessarily on the sphere with radius n
	
	// // to check if scale value is correct
	// double scaled_energy = 0;
	// for(int i = 0; i < transmittedMessage.size(); i++){
	// 	scaled_energy += std::pow(transmittedMessage[i] - encodedMessage[i], 2);
	// }
	// std::cout << scaled_energy << std::endl;

	// remove the projection onto all-zeros cwd
	// vector between all-zeros and current point

	std::vector<double> no0ProjMsg;
	for(int i = 0; i < transmittedMessage.size(); i++){
		no0ProjMsg.push_back(transmittedMessage[i] - encodedMessage[i]);
	}

	// project onto sphere with squared distance n
	std::vector<double> projectedMsg;
	double cur_energy = 0;
	for(int i = 0; i < no0ProjMsg.size(); i++){
		cur_energy += std::pow(no0ProjMsg[i], 2);
	}
	double scale2 = double(no0ProjMsg.size())/cur_energy;

	for(int i = 0; i < no0ProjMsg.size(); i++){
		projectedMsg.push_back(no0ProjMsg[i] * sqrt(scale2));
	}
	// to check if scale value is correct
	double scaled_energy = 0;
	double scaled_dist = 0;
	for(int i = 0; i < projectedMsg.size(); i++){
		scaled_energy += std::pow(projectedMsg[i], 2);
		scaled_dist += std::pow(projectedMsg[i] - encodedMessage[i], 2);
	}
	std::cout << "projected onto sphere" << std::endl;
	std::cout << scaled_energy << std::endl;
	std::cout << scaled_dist << std::endl;

	// we want the distance between the final projected point and all-zeros cwd to be 90.
	// test if we have ML decoding.. 


	std::cout << "concluded simulation" << std::endl;

}

void distanceThreshold_midpoint(codeInformation code){
	std::cout << "running the linearity decoder" << std::endl;

	// set random seed for message generation
	srand(time(NULL));

	// check to make sure the code has valid values
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}
	
	int listSize = 2048;
	
	std::vector<int> puncturedIndices = {};
	// the below are the relevant initializations for low rate
	// decoding, this will be altered for high rate or ZTCC, for example
	FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
	LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg, code.crc);

	std::string dir = "k" + std::to_string(code.k) + "n" + std::to_string(code.n) + "v" + std::to_string(code.v) +
					 "m" + std::to_string(code.crcDeg - 1) + "K" + std::to_string(code.numInfoBits) + "/";
	std::cout << dir << std::endl;
	listDecoder.readNeighborList(dir);
	std::cout << "neighbors read in, beginning trials" << std::endl;

	// all-zeros message
	std::vector<int> originalMessage;
	for(int i = 0; i < code.numInfoBits; i++)
		originalMessage.push_back(0);
	// compute the CRC
	crc_calculation(originalMessage, code.crcDeg, code.crc);
	// encode the message
	std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);
	
	// loop over neighbors until we reach the first cwd whose midpoint does not immediately decode back to all-zeros cwd
	for(int neighbor = 1; neighbor < listDecoder.codewordNeighborList[0].size(); neighbor++){	
		// current neighbor codeword - it's on the sphere
		std::vector<int> neighborCodeword = listDecoder.codewordNeighborList[0][neighbor];

		// find midpoint codeword
		std::vector<double> justPastMidpointCodeword;
		for(int i = 0; i < encodedMessage.size(); i++){
			justPastMidpointCodeword.push_back(0.499*double(neighborCodeword[i]) + 0.501*double(encodedMessage[i]));
		}

		// double cur_energy = 0;
		// for(int i = 0; i < justPastMidpointCodeword.size(); i++){
		// 	cur_energy += std::pow(justPastMidpointCodeword[i], 2);
		// }
		// double scale = double(justPastMidpointCodeword.size())/cur_energy;

		// for(int i = 0; i < justPastMidpointCodeword.size(); i++){
		// 	justPastMidpointCodeword[i] = justPastMidpointCodeword[i] * sqrt(scale);
		// }
		
		LowRateListDecoder::messageInformation standardDecoding = listDecoder.lowRateDecoding(justPastMidpointCodeword);

		// midpoint distance from the all-zeros cwd
		double midMetric = 0;
		double neighborMetric = 0;
		for(int i = 0; i < justPastMidpointCodeword.size(); i++){
			midMetric += std::pow(justPastMidpointCodeword[i] - encodedMessage[i], 2);
			neighborMetric += std::pow(justPastMidpointCodeword[i] - neighborCodeword[i], 2);
		}
		std::cout << midMetric << std::endl;
		std::cout << neighborMetric << std::endl;
		std::cout << "and this is the " << neighbor << "th neighbor" << std::endl;
		
		// if there is decoding error
		if(standardDecoding.message != originalMessage){
			std::cout << "no longer decoding to the all-zeros codeword at distance:" << std::endl;
			// midpoint distance from the all-zeros cwd
			// double midMetric = 0;
			// for(int i = 0; i < justPastMidpointCodeword.size(); i++){
			// 	midMetric += std::pow(justPastMidpointCodeword[i] - encodedMessage[i], 2);
			// }
			// std::cout << midMetric << std::endl;
			// std::cout << "and this is the " << neighbor << "th neighbor" << std::endl;
			break;
		}

	}
	std::cout << "concluded simulation" << std::endl;

}

void linearityDecoding(codeInformation code){
	std::cout << "running the linearity decoder" << std::endl;

	// set random seed for message generation
	srand(time(NULL));

	// check to make sure the code has valid values
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	int listSize = 2048; 
	std::vector<double> EbN0 = {1, 1.5, 2, 2.5, 3, 3.5};
	// std::vector<double> EbN0 = {1};
	vector<double> SNR; 										// SNR is required for noise computations
	double offset = 10 * log10((double)2*code.numInfoBits / (2*(code.numInfoBits+code.crcDeg-1))); 		// real rate of this code *2 in the log
	for (int i=0; i< EbN0.size(); i++)
		SNR.push_back(EbN0[i] + offset);
	int maxNumTrials = 1e4;								// the max trials at each point, the necessary actual value will likely be greater

	std::vector<int> puncturedIndices = {};
	// the below are the relevant initializations for low rate
	// decoding, this will be altered for high rate or ZTCC, for example
	FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
	LowRateListDecoder listDecoder(encodingTrellis, 1e5, code.crcDeg, code.crc);
	LowRateListDecoder listDecoder_linearity(encodingTrellis, listSize, code.crcDeg, code.crc);
	LowRateListDecoder listDecoder_viterbi(encodingTrellis, 1, code.crcDeg, code.crc);

	std::string dir = "k" + std::to_string(code.k) + "n" + std::to_string(code.n) + "v" + std::to_string(code.v) +
					 "m" + std::to_string(code.crcDeg - 1) + "K" + std::to_string(code.numInfoBits) + "/";
	std::cout << dir << std::endl;
	listDecoder.readNeighborList(dir);
	std::cout << "neighbors read in, beginning trials" << std::endl;


	// std::ofstream outputFile_offset_correct; // distance between the received word and the decoding center
	// string filename_offset = "ML_L2048_offset_correctdecoding_EbNo1.txt";
	// outputFile_offset_correct.open(filename_offset, fstream::app);

	// std::ofstream outputFile_decoded_correct; // distance between the received word and the decoded codeword
	// string filename_correct = "ML_L2048_decoded_correctdecoding_EbNo1.txt";
	// outputFile_decoded_correct.open(filename_correct, fstream::app);

	// std::ofstream outputFile_offset_incorrect;
	// string filename_offset2 = "ML_L2048_offset_incorrectdecoding_EbNo1.txt";
	// outputFile_offset_incorrect.open(filename_offset2, fstream::app);

	// std::ofstream outputFile_decoded_incorrect;
	// string filename_correct2 = "ML_L2048_decoded_incorrectdecoding_EbNo1.txt";
	// outputFile_decoded_incorrect.open(filename_correct2, fstream::app);


	// running the simulations. in this example, we are simulation to collect TFR data
	for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++){
		double snr = SNR[snrIndex];
		double linearityMeanListSize = 0;
		double standardMeanListSize = 0;
		double standardMeanListSize2048 = 0;
		double multiMeanListSize = 0;

		int linearityNumErrors = 0;
		int standardNumErrors = 0;
		int standardNumErrors2048 = 0;
		int multiNumErrors = 0;

		int linearityListSizeExceeded = 0;
		int standardListSizeExceeded = 0;
		int standardListSizeExceeded2048 = 0;
		int multiListSizeExceeded = 0;

		int nontrivialLinearityDecoding = 0;

		for(int numTrials = 0; numTrials < maxNumTrials; numTrials++){
			// the following lines generate a message that we can decode
			std::vector<int> originalMessage = generateRandomCRCMessage(code);
			std::vector<double> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices);
			// project onto sphere with squared distance n
			double cur_energy = 0;
			for(int i = 0; i < transmittedMessage.size(); i++){
				cur_energy += std::pow(transmittedMessage[i], 2);
			}
			double scale = double(transmittedMessage.size())/cur_energy;

			for(int i = 0; i < transmittedMessage.size(); i++){
				transmittedMessage[i] = transmittedMessage[i] * sqrt(scale);
			}
			// std::cout << transmittedMessage.size() << std::endl;
			/// to check if scale value is correct
			// double scaled_energy = 0;
			// for(int i = 0; i < transmittedMessage.size(); i++){
			// 	scaled_energy += std::pow(transmittedMessage[i], 2);
			// }


			// std::cout << "We are at: " << numTrials << "th trial" << std::endl;
			// decoding the transmitted message. a different decoding function may be 
			// employed depending on what information is desired
			//std::cout << "begin linearity decoding" << std::endl;
			LowRateListDecoder::messageInformation linearityDecoding = listDecoder_linearity.linearityDecoder(transmittedMessage);
			LowRateListDecoder::messageInformation standardDecoding = listDecoder.lowRateDecoding(transmittedMessage);
			// 
			LowRateListDecoder::messageInformation standardDecoding_viterbi = listDecoder_viterbi.lowRateDecoding(transmittedMessage);
			int numBall = 256;
			LowRateListDecoder::messageInformation multiLinearityDecoding = listDecoder.multiLinearityDecoder(transmittedMessage, numBall);

			// linearity decoding results
			// if(linearityDecoding.listSizeExceeded){
			// 	linearityListSizeExceeded++;
			// }
			// else if (linearityDecoding.message != originalMessage){
			// 	linearityNumErrors ++;
			// 	linearityMeanListSize += (double)linearityDecoding.listSize;
			// }

			// // list size = 1 standard decoding results
			// if(standardDecoding.listSizeExceeded){
			// 	standardListSizeExceeded++;
			// }
			// else if (standardDecoding.message != originalMessage){
			// 	standardNumErrors ++;
			// 	standardMeanListSize += (double)standardDecoding.listSize;
			// }

			// // when both standard and linearity make the same error
			// if ((linearityDecoding.message != originalMessage) && (standardDecoding.message != originalMessage)){
			// 	nontrivialLinearityDecoding ++;
			// 	// std::vector<int> foundCodeword = linearityDecoding.path; // yoinking path since we're not using it- it is actually the codeword
			// 	// double foundMetric = 0;
			// 	// for(int i = 0; i < transmittedMessage.size(); i++){
			// 	// 	foundMetric += std::pow(transmittedMessage[i] - foundCodeword[i], 2);
			// 	// }
			// 	// std::cout << "nontrivial linearity decoding metric is: " << foundMetric << std::endl;
			// }

			// // list size = 2048 standard decoding results
			// if(standardDecoding2048.listSizeExceeded){
			// 	standardListSizeExceeded2048++;
			// }
			// else if (standardDecoding2048.message != originalMessage){
			// 	standardNumErrors2048 ++;
			// 	standardMeanListSize2048 += (double)standardDecoding2048.listSize;
			// }
			
			// multi-ball decoder
			if(multiLinearityDecoding.listSizeExceeded){
				multiListSizeExceeded++;
			}
			else if (multiLinearityDecoding.message != originalMessage){
				multiNumErrors ++;
				multiMeanListSize += (double)multiLinearityDecoding.listSize;
			}


			/*
			// //// if there is decoding error
			// if(linearityDecoding.message != originalMessage){
			// if the decoder made the ML decision regardless of whether it was correct
			// if(linearityDecoding.message != standardDecoding.message){ // not ML decision
				// if(standardDecoding.message != originalMessage)
				// 	std::cout << "the ML decoder also made an error" << std::endl;
				// else
				// 	std::cout << "the ML decoder did NOT make an error here" << std::endl;
				// numErrors++;
				// std::vector<int> foundCodeword = linearityDecoding.path; // yoinking path since we're not using it- it is actually the codeword
				// std::vector<int> originalCodeword = encodingTrellis.encoder(originalMessage);
				// double foundMetric = 0;
				// // double originalMetric = 0;
				// for(int i = 0; i < transmittedMessage.size(); i++){
				// 	foundMetric += std::pow(transmittedMessage[i] - foundCodeword[i], 2);
				// 	// originalMetric += std::pow(transmittedMessage[i] - originalCodeword[i], 2);
				// }
				// std::cout << "foundMetric: " << foundMetric << std::endl;
				// std::cout << "original metric: " << originalMetric << std::endl;
				// std::cout << "distance to the center of the decoding sphere: " << linearityDecoding.metric << std::endl;

				// // distance between the received word and the decoded codeword
				// outputFile_decoded_incorrect << foundMetric << std::endl;
				// // offset distance between the received word and the decoding center
				// outputFile_offset_incorrect << linearityDecoding.offset << std::endl;

				// used to confirm that the codewords and messages being returned are internally consistent- which they are
				//if(codeword != testCodeword){
				//	std::cout << "issue- the codeword corresponding to the message according to the decoder and encoder differ" << std::endl;
				//	std::cout << "encoder: " << std::endl;
				//	print_int_vector(testCodeword);
				//	std::cout << "decoder: " << std::endl;
				//	print_int_vector(codeword);
				//}
				
			// }
			// else{ //correct decoding
			// 	outputFile_decoded_correct<< linearityDecoding.metric << std::endl;
			// 	outputFile_offset_correct << linearityDecoding.offset << std::endl;
			// }
			// meanListSize += (double)linearityDecoding.listSize;
			*/
			
		}
		// outputting results at each SNR / EbN0 point. in this case, since we're working with TFR, the ratio
		// of incorrect decodings to total decodings is used.
		std::cout << "at Eb/N0 = " << EbN0[snrIndex] << std::endl;
		
		// std::cout << "Linearity Decoder" << std::endl;
		// std::cout << "mean list size: " << (double)linearityMeanListSize/(maxNumTrials - linearityListSizeExceeded) << std::endl;
		// std::cout << "number of times list size exceeded: " << linearityListSizeExceeded << std::endl;
		// std::cout << "UER: " << (double)linearityNumErrors/(maxNumTrials - linearityListSizeExceeded) << std::endl;
		// std::cout << "TFR: " << (double)(linearityNumErrors + linearityListSizeExceeded)/maxNumTrials << std::endl;
		
		// std::cout << "Standard Decoder with L=1" << std::endl;
		// std::cout << "mean list size: " << (double)standardMeanListSize/(maxNumTrials - standardListSizeExceeded) << std::endl;
		// std::cout << "number of times list size exceeded: " << standardListSizeExceeded << std::endl;
		// std::cout << "UER: " << (double)standardNumErrors/(maxNumTrials - standardListSizeExceeded) << std::endl;
		// std::cout << "TFR: " << (double)(standardNumErrors + standardListSizeExceeded)/maxNumTrials << std::endl;
		
		// std::cout << "nontrivial Linearity decoding frequency: "<< (double)nontrivialLinearityDecoding/(maxNumTrials - linearityListSizeExceeded)  << std::endl;

		// std::cout << "Standard Decoder with L=2048" << std::endl;
		// std::cout << "mean list size: " << (double)standardMeanListSize2048/(maxNumTrials - standardListSizeExceeded2048) << std::endl;
		// std::cout << "number of times list size exceeded: " << standardListSizeExceeded2048 << std::endl;
		// std::cout << "UER: " << (double)standardNumErrors2048/(maxNumTrials - standardListSizeExceeded2048) << std::endl;
		// std::cout << "TFR: " << (double)(standardNumErrors2048 + standardListSizeExceeded2048)/maxNumTrials << std::endl;
		
		std::cout << "multi-ball linearity decoder" << std::endl;
		// std::cout << "mean list size: " << (double)multiMeanListSize/(maxNumTrials - multiListSizeExceeded) << std::endl;
		std::cout << "number of times list size exceeded: " << multiListSizeExceeded << std::endl;
		std::cout << "UER: " << (double)multiNumErrors/(maxNumTrials - multiListSizeExceeded) << std::endl;
		std::cout << "TFR: " << (double)(multiNumErrors + multiListSizeExceeded)/maxNumTrials << std::endl;
		

		// it may be the case that we want to write the outputs to a file, whether it be for further processing
		// or convenience. the following code snipped illustrates how this would be done. typically, the opening
		// and closing of the file, and declaration of the file name would be done once each outside the loop,
		// but they are included here to keep things organized.
		/*
		std::ofstream outputFile;
		filename = "distance_spectra.txt";
		outputFile.open(filename, fstream::app);
		outputFile << "example information to be saved";
		outputFile.close();
		*/
	}
	std::cout << "concluded simulation" << std::endl;
	// outputFile_decoded_correct.close();
	// outputFile_offset_correct.close();
	// outputFile_decoded_incorrect.close();	
	// outputFile_offset_incorrect.close();

}

// describe function here
void generateNeighborList(codeInformation code){
	std::cout << "generating the list of neighbors" << std::endl;

	// set random seed for message generation
	srand(time(NULL));

	// check to make sure the code has valid values
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	// this will decide how long we want each list to be
	int listSize = 1;


	// the below are the relevant initializations for low rate
	// decoding, this will be altered for high rate or ZTCC, for example
	FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
	LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg, code.crc);

	// initializing the zeros message to find neighbors
	std::vector<double> allZerosMessage;
	int messageLength = (code.numInfoBits + code.crcDeg - 1) * code.n / code.k;
	std::cout << "message length: " << messageLength << std::endl;
	for(int i = 0; i < messageLength; i++)
		allZerosMessage.push_back(1);
	std::string dir = "k" + std::to_string(code.k) + "n" + std::to_string(code.n) + "v" + std::to_string(code.v) +
					 "m" + std::to_string(code.crcDeg - 1) + "K" + std::to_string(code.numInfoBits) + "/";
	std::cout << dir << std::endl;
	listDecoder.generateNeighborList(allZerosMessage, dir, listSize);
	//listDecoder.readNeighborList("test/");

}

// finds both the list size on the standard trellis and the list size 
// when considering only tailbiting codewords
void tb_nontb_list_sizes(codeInformation code){
	std::cout << "beginning tb_nontb_list_sizes" << std::endl;

	// set random seed for message generation
	srand(time(NULL));

	// check to make sure the code has valid values
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	// initalizing the values to be used for simulation
	int listSize = 1e7;
	std::vector<double> EbN0 = {0.5, 0, -0.5, -1, -1.5, -2, -2.5, -3, -3.5, -4};
	vector<double> SNR; 										// SNR is required for noise computations
	double offset = 10 * log10((double)2*64 / (2*71)); 		// real rate of this code is 32/512
	for (int i=0; i< EbN0.size(); i++)
		SNR.push_back(EbN0[i] + offset);
	int maxNumTrials = 1e5;								// the max trials at each point, the necessary actual value will likely be greater

	std::vector<int> puncturedIndices = {};
	// the below are the relevant initializations for low rate
	// decoding, this will be altered for high rate or ZTCC, for example
	FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
	LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg, code.crc);

	// running the simulations. in this example, we are simulation to collect TFR data
	for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++){
		double snr = SNR[snrIndex];
		double meanListSize = 0;
		double meanTBListSize = 0;
		for(int numTrials = 0; numTrials < maxNumTrials; numTrials++){
			// the following lines generate a message that we can decode
			std::vector<int> originalMessage = generateRandomCRCMessage(code);
			std::vector<double> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices);
			// std::cout << "We are at: " << numTrials << "th trial" << std::endl;
			// decoding the transmitted message. a different decoding function may be 
			// employed depending on what information is desired
			//std::cout << "begin linearity decoding" << std::endl;
			LowRateListDecoder::messageInformation serialDecoding = listDecoder.lowRateDecoding(transmittedMessage);

			if(serialDecoding.listSizeExceeded)
				std::cout << "MAJOR ISSUE- LIST SIZE EXCEEDED" << std::endl;

			meanListSize += (double)serialDecoding.listSize/maxNumTrials;
			meanTBListSize += (double)serialDecoding.TBListSize/maxNumTrials;

			
			
		}
		// outputting results at each SNR / EbN0 point. in this case, since we're working with TFR, the ratio
		// of incorrect decodings to total decodings is used.
		std::cout << "at Eb/N0 = " << EbN0[snrIndex] << std::endl;
		std::cout << "mean list size: " << meanListSize << std::endl;
		std::cout << "TB mean list size: " << meanTBListSize << std::endl;

		// it may be the case that we want to write the outputs to a file, whether it be for further processing
		// or convenience. the following code snipped illustrates how this would be done. typically, the opening
		// and closing of the file, and declaration of the file name would be done once each outside the loop,
		// but they are included here to keep things organized.
		/*
		std::ofstream outputFile;
		filename = "distance_spectra.txt";
		outputFile.open(filename, fstream::app);
		outputFile << "example information to be saved";
		outputFile.close();
		*/
	}
}

// describe function here
void example_function(codeInformation code){
	std::cout << "beginning example_function" << std::endl;

	// set random seed for message generation
	srand(time(NULL));

	// check to make sure the code has valid values
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	// initalizing the values to be used for simulation
	int listSize = 512;
	std::vector<double> EbN0 = {1, 1.5 , 2, 2.5, 3, 3.5};
	vector<double> SNR; 										// SNR is required for noise computations
	double offset = 10 * log10((double)2*32 / (double)512); 		// real rate of this code is 32/512
	for (int i=0; i< EbN0.size(); i++)
		SNR.push_back(EbN0[i] + offset);
	std::vector<int> puncturedIndices = {47, 60, 129, 504};
	int maxNumTrials = 1e4;								// the max trials at each point, the necessary actual value will likely be greater


	// the below are the relevant initializations for low rate
	// decoding, this will be altered for high rate or ZTCC, for example
	FeedforwardTrellis encodingTrellis(code.k, code.n, code.v, code.numerators);
	LowRateListDecoder highRateListFinder(encodingTrellis, listSize, code.crcDeg, code.crc);
	LowRateListDecoder TBGuaranteeListFinder(encodingTrellis, listSize, code.crcDeg, code.crc);
	LowRateListDecoder fullNeighborFinder(encodingTrellis, 16, code.crcDeg, code.crc);
	LowRateListDecoder fullNeighborFinder2(encodingTrellis, 16, code.crcDeg, code.crc);

	// puncture the 4 bits from length-516 noisy message
	std::vector<double> allZeroMessage;
	for (int i = 0; i < 516; i++){
		if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
			allZeroMessage.push_back(1);
		}
		else{
			allZeroMessage.push_back(0);
		}
	}

	std::vector<double> allOneMessage = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 0, 1, -1, -1, -1, -1, 1, 1, 1, -1, 1, 1, 1, 0, 1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 0, -1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	// for (int i = 0; i < 516; i++){
	// 	if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
	// 		allOneMessage.push_back(-1);
	// 	}
	// 	else{
	// 		allOneMessage.push_back(0);
	// 	}
	// }


	
	// Find the closest ${listSize} high-rate codewords
	LowRateListDecoder::messageInformation allZeroNeighbors = highRateListFinder.linearityLowRateListSearching(allZeroMessage);
	// Run the list decoder once to identify all-zero neighbor list
	std::vector<std::vector<int>> neighbor_list(listSize, vector<int>(516, 0)); // A list of codewords neighboring all-zero codeword
	std::vector<std::vector<int>> neighbor_msgs(listSize, vector<int>(43, 0)); // A list of codewords neighboring all-zero codeword
	std::vector<std::vector<int>> path_ie(listSize, vector<int>(2,0));
	neighbor_list = allZeroNeighbors.neighbor_cwds;
	neighbor_msgs = allZeroNeighbors.neighbor_msgs;
	path_ie = allZeroNeighbors.path_ie;

	// Find the closest ${listSize} TB passing codewords
 	LowRateListDecoder::messageInformation TBguaranteeNeighbors = TBGuaranteeListFinder.TBGuaranteedListSearching(allZeroMessage);
	std::vector<std::vector<int>> TBguarantee_neighbor_list(listSize, vector<int>(516, 0)); // A list of codewords neighboring all-zero codeword
	std::vector<std::vector<int>> TBguarantee_neighbor_msgs(listSize, vector<int>(43, 0)); // A list of codewords neighboring all-zero codeword
	std::vector<std::vector<int>> TBguarantee_path_ie(listSize, vector<int>(2,0));
	TBguarantee_neighbor_list = TBguaranteeNeighbors.neighbor_cwds;
	TBguarantee_neighbor_msgs = TBguaranteeNeighbors.neighbor_msgs;
	TBguarantee_path_ie = TBguaranteeNeighbors.path_ie;

	// Find all neighbor lists 
	std::vector<LowRateListDecoder::messageInformation> fullNeighbors = fullNeighborFinder.fullListSearching(allZeroMessage);
	std::vector<LowRateListDecoder::messageInformation> fullNeighbors_ones = fullNeighborFinder2.fullListSearching(allOneMessage);




	/////////// Initialize the list decoder with neighbor_list
	LowRateListDecoder linearityListDecoder(encodingTrellis, 1, code.crcDeg, code.crc, neighbor_list, neighbor_msgs, path_ie);
	LowRateListDecoder TBlinearityListDecoder(encodingTrellis, listSize, code.crcDeg, code.crc, TBguarantee_neighbor_list, TBguarantee_neighbor_msgs, TBguarantee_path_ie);
	LowRateListDecoder listDecoder(encodingTrellis, listSize, code.crcDeg, code.crc);

	// running the simulations. in this example, we are simulation to collect TFR data
	for(int snrIndex = 0; snrIndex < SNR.size(); snrIndex++){
		double snr = SNR[snrIndex];
		int linearErrors = 0;
		int TBlinearErrors = 0;
		int serialErrors = 0;
		int offsetSerialErrors = 0;
		int linearListSizeExceeded = 0;
		int fullNeighborErrors = 0;
		for(int numTrials = 0; numTrials < maxNumTrials; numTrials++){
			// the following lines generate a message that we can decode
			std::vector<int> originalMessage = generateRandomCRCMessage(code);
			std::vector<double> transmittedMessage = generateTransmittedMessage(originalMessage, encodingTrellis, snr, puncturedIndices);
			// std::cout << "We are at: " << numTrials << "th trial" << std::endl;
			// decoding the transmitted message. a different decoding function may be 
			// employed depending on what information is desired
			//std::cout << "begin linearity decoding" << std::endl;
			LowRateListDecoder::messageInformation decoding = linearityListDecoder.linearityLowRateDecoding(transmittedMessage);			// call listDecoder (with listsize 1) to denoise
			LowRateListDecoder::messageInformation TBguaranteeDecoding = TBlinearityListDecoder.TBGuaranteedLinearityDecoding(transmittedMessage);
			LowRateListDecoder::messageInformation serialDecoding = listDecoder.lowRateDecoding(transmittedMessage);
			LowRateListDecoder::messageInformation offsetSerialDecoding = listDecoder.offsetLowRateDecoding(transmittedMessage);
			LowRateListDecoder::messageInformation fullNeighborDecoding = fullNeighborFinder.fullListDecoding(transmittedMessage);
			

			// since we're interested in TFR, we validate that the decoded message matches the original message
			if(decoding.message != originalMessage)
				linearErrors++;
			if (TBguaranteeDecoding.message != originalMessage)
				TBlinearErrors++;
			if(decoding.listSizeExceeded)
				linearListSizeExceeded++;
			if(serialDecoding.message != originalMessage)
				serialErrors++;
			if(offsetSerialDecoding.message != originalMessage)
				offsetSerialErrors++;
			if(fullNeighborDecoding.message != originalMessage)
				fullNeighborErrors++;

			if((numTrials+ 1)%1000 == 0){
				std::cout << "in " << numTrials + 1 << "trials, we have: ---------------------------------" << std::endl;
				std::cout << "linear errors: " << linearErrors << std::endl;
				std::cout << "serial errors: " << serialErrors << std::endl;
				std::cout << "offset serial errors: " << offsetSerialErrors << std::endl;
				std::cout << "TB guarantee linear errors: " << TBlinearErrors << std::endl;
				std::cout << "Full List errors: " << fullNeighborErrors << std::endl;
			}
			
		}
		// outputting results at each SNR / EbN0 point. in this case, since we're working with TFR, the ratio
		// of incorrect decodings to total decodings is used.
		std::cout << "at Eb/N0 = " << EbN0[snrIndex] << ", we achieved " << linearErrors << " incorrect decodings (from the linearity decoder) in " << maxNumTrials << " trials" << std::endl;
		std::cout << "this yields a TFR of " << (double)linearErrors/maxNumTrials << " for the linear decoder" << std::endl;
		std::cout << "this yields a TFR of " << (double)TBlinearErrors/maxNumTrials << " for the TB guaranteed linear decoder" << std::endl;
		std::cout << "this yields a TFR of " << (double)serialErrors/maxNumTrials << " for the serial decoder" << std::endl;
		std::cout << "this yields a TFR of " << (double)offsetSerialErrors/maxNumTrials << " for the offset serial decoder" << std::endl;
		std::cout << "this yields a TFR of " << (double)fullNeighborErrors/maxNumTrials << " for the full neighbor decoder" << std::endl;
		std::cout << "we had a total of " << linearListSizeExceeded << " trials where the linear list size was exceeded" << std::endl;

		// it may be the case that we want to write the outputs to a file, whether it be for further processing
		// or convenience. the following code snipped illustrates how this would be done. typically, the opening
		// and closing of the file, and declaration of the file name would be done once each outside the loop,
		// but they are included here to keep things organized.
		/*
		std::ofstream outputFile;
		filename = "distance_spectra.txt";
		outputFile.open(filename, fstream::app);
		outputFile << "example information to be saved";
		outputFile.close();
		*/
	}
}

// this generates a random binary string of length code.numInfoBits, and appends the appropriate CRC bits
std::vector<int> generateRandomCRCMessage(codeInformation code){
	std::vector<int> message;
	for(int i = 0; i < code.numInfoBits; i++)
		message.push_back(rand()%2);
	// compute the CRC
	crc_calculation(message, code.crcDeg, code.crc);
	return message;
}

// this takes the message bits, including the CRC, and encodes them using the trellis,
// then adds noise and punctures the result
std::vector<double> generateTransmittedMessage(std::vector<int> originalMessage, FeedforwardTrellis encodingTrellis, double snr, std::vector<int> puncturedIndices, bool noiseless){
	// encode the message
	std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);
	// add noise
	std::vector<double> finalMessage;
	if(noiseless){
		for(int i = 0; i < encodedMessage.size(); i++)
			finalMessage.push_back((double)encodedMessage[i]);
	}
	else
		finalMessage = addNoise(encodedMessage, snr);
	// puncture the bits. it is more convenient to puncture on this side than on the 
	// decoder, so we insert zeros which provide no information to the decoder
	for(int index = 0; index < puncturedIndices.size(); index++)
		finalMessage[puncturedIndices[index]] = 0;
	return finalMessage;
}

// adds noise at the given snr to each bit, based on a normal distribution
std::vector<double> addNoise(std::vector<int> encodedMsg, double SNR){
	std::vector<double> noisyMsg;

	/*
	// the below lines break the noise generator for some compilers
	std::random_device rd;
	std::default_random_engine generator;
	generator.seed( rd() ); //Now this is seeded differently each time.
	*/
	double variance = pow(10.0, -SNR/10.0);
	double sigma = sqrt(variance);
	normal_distribution<double> distribution(0.0, sigma);
	
	// cout << variance << endl;

	for (int i = 0; i < encodedMsg.size(); i++) {
		noisyMsg.push_back(encodedMsg[i] + distribution(generator));
	}
	return noisyMsg;
}

void lowrate_find_reciprocity(codeInformation code){
	// for a given codeword, find the list size it takes to go back to all-zero codeword
	std::cout << "working on rate-1/12 SLVD" << std::endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "TBCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	std::cout << "working on code " << filename.substr(0, filename.size() - 4) << std::endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	FeedforwardTrellis lowrate_trellis(code.k, code.n, code.v, code.numerators);
	std::cout << "trellis construction complete" << std::endl;

	int listSize = std::pow(2,18);
	std::vector<int> listSizes;

	LowRateListDecoder LowRateListDecoder(lowrate_trellis, listSize, code.crcDeg, code.crc); 

	std::vector<int> originalMessage;
	for (int i = 0; i < code.numInfoBits; i++){
		originalMessage.push_back(0);
	}
	
	crc_calculation(originalMessage, code.crcDeg, code.crc);

	std::vector<int> encodedMessage = lowrate_trellis.encoder(originalMessage);		

	// puncture the 4 bits from length-516 noisy message
	std::vector<double> punc_msg;
	for (int i = 0; i < encodedMessage.size(); i++){
		if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
			punc_msg.push_back(encodedMessage[i]);
		}
		else{
			punc_msg.push_back(0);
		}
	}

	// find 2^18 nearest codewords
	// a new decoding function that prints the codeword at listSize = i

	// NOTE: WE ARE NOW COLLECTING THE MIDPOINT BETWEEN THE CODEWORDS AND THE ALL ZEROS CODEWORDS
	LowRateListDecoder::messageInformation lowrate_codeword = LowRateListDecoder.lowRateDecoding_listSize(punc_msg);
	cout << "Collected " << listSize << " codewords." << endl;

	// read in the codewords we just saved
	std::vector<std::vector<double>> codewords;
	std::ifstream file_in("codewords.txt");
	std::string line;
	while (std::getline(file_in, line))
    {
        std::istringstream ss(line);
        codewords.emplace_back(std::istream_iterator<double>(ss), std::istream_iterator<double>());
    }
	cout << "Done reading in codewords " << endl;

	std::ofstream outputFile;
	filename = "reciprocity_list_sizes.txt";
	outputFile.open(filename);

	for(int i = 0; i < 18; i++){
		std::vector<double> codeword = codewords[(int)std::pow(2,i)];
		double diff = 0;
		for(int j = 0; j < codeword.size(); j++)
			diff += std::pow(1-codeword[j],2);
		std::cout << "at " << (int)std::pow(2,i) << ", we had diff " << diff << std::endl;
	}

	// find reciporcity list size of each codeword
	for (int i=0; i<listSize; i++){
		std::vector<double> cur_codeword;
		for (int j=0; j<encodedMessage.size(); j++){
			cur_codeword.push_back(codewords[i][j]);
		}
		// print_double_vector(cur_codeword);
		LowRateListDecoder::messageInformation lowrate_decoded = LowRateListDecoder.lowRateDecoding_recip(cur_codeword);
		
		outputFile << lowrate_decoded.listSize - 1 << std::endl;
		//std::cout << "Reciprocol list size is : " << lowrate_decoded.listSize << std::endl;
	}


	outputFile.close();
	
}



void q_func_validation(codeInformation code){
	// observe the distribution of decoding weights to
	// verify how we are computing the q_function
	std::cout << "working on rate-1/12 SLVD" << std::endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "TBCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	std::cout << "working on code " << filename.substr(0, filename.size() - 4) << std::endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	FeedforwardTrellis lowrate_trellis(code.k, code.n, code.v, code.numerators);
	std::cout << "trellis construction complete" << std::endl;

	int listSize = std::pow(2,18);
	
	listSize = 1e3;
	std::cout << "list size: " << listSize << std::endl;

	LowRateListDecoder LowRateListDecoder(lowrate_trellis, listSize, code.crcDeg, code.crc); 

	vector<int> breakpoints = {4, 72, 108};
	
	vector<double> EbN0 = {1, 1.5, 2, 2.5, 3, 3.5, 50};			// will be running from 1 to 3.5 in increments of 0.5
	vector<double> SNR; 
	double offset = 10 * log10((double)32 / (double)512); // real rate of this code is 32/512
	//double offset = 10 * log10((double)2*32 / (double)512); // ?? should the rate be doubled? seems to agree better with jacob's results
	std::cout << "offset: " << offset << std::endl;
	for (int i=0; i< EbN0.size(); i++){
		SNR.push_back(EbN0[i] + offset);
	}
	for(int i = 0; i < SNR.size(); i++)
		std::cout << "at EbN0 " << EbN0[i] << ", we have SNR: " << SNR[i] << std::endl;
	int maxNumTrials = 1e4;
	/*
	for(int i = 0; i < EbN0.size(); i++){
		vector<int> dist = {0, 0, 0};
		double meanNoise = 0;
		double meanSquaredNoise = 0;
		for(int numTrials = 0; numTrials < maxNumTrials; numTrials++){
			std::vector<int> originalMessage;
			for (int i = 0; i < code.numInfoBits; i++){
				originalMessage.push_back(rand()%2);
			}
			
			crc_calculation(originalMessage, code.crcDeg, code.crc);

			std::vector<int> encodedMessage = lowrate_trellis.encoder(originalMessage);		
			std::vector<double> noisyMessage = addNoise(encodedMessage, SNR[i]);
			
			// puncture the 4 bits from length-516 noisy message
			std::vector<double> punc_msg;
			for (int i = 0; i < noisyMessage.size(); i++){
				if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
					punc_msg.push_back(noisyMessage[i]);
				}
				else{
					punc_msg.push_back(0);
				}
			}
			double localMeanNoise = 0;
			double localMeanSquaredNoise = 0;
			for(int i = 0; i < punc_msg.size(); i++){
				localMeanNoise += std::abs(punc_msg[i] - encodedMessage[i]);
				localMeanSquaredNoise += std::pow(punc_msg[i] - encodedMessage[i],2);
			}
			meanNoise += localMeanNoise / maxNumTrials;
			meanSquaredNoise += localMeanSquaredNoise / maxNumTrials;
			for(int j = 0; j < breakpoints.size(); j++){
				if(localMeanNoise > breakpoints[j])
					dist[j]++;
			}
			
		}
		std::cout << "at ebno: " << EbN0[i] << std::endl;
		for(int j = 0; j < dist.size(); j++){
			std::cout << "we had " << dist[j] << " messages exceed the weight " << breakpoints[j] << std::endl;
		}
		std::cout << "when using the euclidean distance, the distance is " << meanNoise << std::endl;
		std::cout << "when using the squared distance, the distance is " << meanSquaredNoise << std::endl;

	}
	*/
	std::vector<int> zerosCodeword;
	for (int i = 0; i < code.numInfoBits; i++){
		zerosCodeword.push_back(0);
	}
	
	crc_calculation(zerosCodeword, code.crcDeg, code.crc);

	std::vector<int> zerosEncodedMessage = lowrate_trellis.encoder(zerosCodeword);

	std::vector<double> zerosNoisyMessage = addNoise(zerosEncodedMessage, 1);
	std::vector<double> allZerosMessage;
	for(int i = 0; i < zerosEncodedMessage.size(); i++)
		allZerosMessage.push_back((double)zerosEncodedMessage[i]);

	std::cout << "check 1" << std::endl;
	// puncture the 4 bits from length-516 noisy message
	std::vector<double> noisy_punc_msg;
	std::vector<double> punc_msg;
	for (int i = 0; i < zerosNoisyMessage.size(); i++){
		if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
			noisy_punc_msg.push_back(zerosNoisyMessage[i]);
		}
		else{
			noisy_punc_msg.push_back(0);
		}
	}
	for (int i = 0; i < zerosNoisyMessage.size(); i++){
		if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
			punc_msg.push_back(allZerosMessage[i]);
		}
		else{
			punc_msg.push_back(0);
		}
	}
	std::cout << "check 1.1" << std::endl;
	std::vector<double> noisyDistances = LowRateListDecoder.lowRateDistanceDistribution(noisy_punc_msg);
	std::cout << "check 1.10" << std::endl;
	std::vector<double> distances = LowRateListDecoder.lowRateDistanceDistribution(punc_msg);
	std::cout << "distance distribution of the all zeros message with noise------------------------------------------------" << std::endl;
	print_double_vector(noisyDistances);
	std::cout << "distance distribution of the all zeros codeword with no noise------------------------------------------------" << std::endl;
	print_double_vector(distances);
	std::cout << "printing diff---------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "check 1.11" << std::endl;
	for(int i = 0; i < noisyDistances.size(); i++)
		noisyDistances[i] -= distances[i];
	print_double_vector(noisyDistances);
	std::cout << "check 1.2" << std::endl;
}

void distance_distribution(codeInformation code){
	// get the distance distribution from the all-zeros codeword
	// of all codewords up to the list size, including
	// non-tb, non-crc codewords
	std::cout << "working on rate-1/12 SLVD" << std::endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "TBCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	std::cout << "working on code " << filename.substr(0, filename.size() - 4) << std::endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	FeedforwardTrellis lowrate_trellis(code.k, code.n, code.v, code.numerators);
	std::cout << "trellis construction complete" << std::endl;

	int listSize = std::pow(2,18)+1;
	std::cout << "list size: " << listSize << std::endl;

	LowRateListDecoder LowRateListDecoder(lowrate_trellis, listSize, code.crcDeg, code.crc); 

	std::vector<int> originalMessage;
	for (int i = 0; i < code.numInfoBits; i++){
		originalMessage.push_back(0);
	}
	
	crc_calculation(originalMessage, code.crcDeg, code.crc);

	std::vector<int> encodedMessage = lowrate_trellis.encoder(originalMessage);		

	// puncture the 4 bits from length-516 noisy message
	std::vector<double> punc_msg;
	for (int i = 0; i < encodedMessage.size(); i++){
		if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
			punc_msg.push_back(encodedMessage[i]);
		}
		else{
			punc_msg.push_back(0);
		}
	}

	std::vector<double> distance_distribution = LowRateListDecoder.lowRateDistanceDistribution(punc_msg);
	std::cout << "simulation complete- printing distances" << std::endl;
	std::ofstream outputFile;
	filename = "distance_spectra.txt";
	outputFile.open(filename, fstream::app);
	for(int i = 0; i < distance_distribution.size() - 1; i++){
		outputFile << distance_distribution[i] << ", ";
	}
	outputFile << distance_distribution[distance_distribution.size() - 1];
	outputFile.close();
}

void pointTest(){
	std::vector<int> nums = get_point(8, 12);
	for(int i = 0; i < 12; i++){
		std::cout << nums[i] << std::endl;
	}
}

void lowrate_test1(codeInformation code){
	// run the list decoder on two different noiseless TB and CRC-checking codewords up to list size 2048  
	// record the 2048 sequences that are decoded.  
	// In particular, prepare the list of 2048 distances.  
	// Plot how distance grows as a function of list size.
	std::cout << "working on rate-1/12 SLVD" << std::endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "TBCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	std::cout << "working on code " << filename.substr(0, filename.size() - 4) << std::endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	FeedforwardTrellis lowrate_trellis(code.k, code.n, code.v, code.numerators);
	std::cout << "trellis construction complete" << std::endl;

	int listSize = 1e5;

	LowRateListDecoder LowRateListDecoder(lowrate_trellis, listSize, code.crcDeg, code.crc); 

	std::vector<int> originalMessage;
	for (int i = 0; i < code.numInfoBits; i++){
		originalMessage.push_back(0);
	}
	
	crc_calculation(originalMessage, code.crcDeg, code.crc);

	std::vector<int> encodedMessage = lowrate_trellis.encoder(originalMessage);		

	// puncture the 4 bits from length-516 noisy message
	std::vector<double> punc_msg;
	for (int i = 0; i < encodedMessage.size(); i++){
		if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
			punc_msg.push_back(encodedMessage[i]);
		}
		else{
			punc_msg.push_back(0);
		}
	}
	std::vector<std::vector<double>> codewords;

	codewords.push_back({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, -1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 0, -1, -1, 1, 1, -1, -1, 1, 1, 1, -1, 1});
	codewords.push_back({-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 0, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1});
	codewords.push_back({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 0, -1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1});
	codewords.push_back({1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1});
	codewords.push_back({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, 1, -1, -1, 1, 0, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1});
	codewords.push_back({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, -1, -1, 0, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1});
	codewords.push_back({-1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, 0, -1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 0, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
	codewords.push_back({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 0, 1, -1, 1, -1, -1, 1, -1, -1, -1, 1, 1});
	codewords.push_back({1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, 0, -1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 0, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
	codewords.push_back({-1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 0, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
	codewords.push_back({-1, 1, -1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, 0, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
	codewords.push_back({1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1});
	codewords.push_back({1, 1, -1, 1, -1, 1, -1, 1, 1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 0, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
	codewords.push_back({-1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 0, -1, -1, -1, -1, 1, 1, 1, -1, 1, 1, 1});
	codewords.push_back({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, 1, 1, -1, 0, 1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1});
	codewords.push_back({1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 0, -1, -1, -1, -1, 1, 1, 1, -1, 1, 1, 1});
	codewords.push_back({-1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, 0, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});


	for(int i = 0; i < codewords.size(); i++){
		int dist = 0;
		for(int j = 0; j < punc_msg.size(); j++){
			dist += pow(punc_msg[j] - codewords[i][j], 2);
		}
		std::cout << "distance: " << dist << std::endl;
	}

	std::vector<int> listSizes;
	for(int i = 0; i < codewords.size(); i++){
		LowRateListDecoder::messageInformation lowrate_decoded = LowRateListDecoder.lowRateDecoding(codewords[i]);
		std::cout << "list size: " << lowrate_decoded.listSize << std::endl;
		std::cout << "metric: " << lowrate_decoded.metric << std::endl;
		if(lowrate_decoded.listSizeExceeded){
			std::cout << "list size exceeded- bug likely" << std::endl;
		}
		if(lowrate_decoded.message != originalMessage){
			std::cout << "our messages do NOT match- bug likely" << std::endl;
		}
		while(listSizes.size() < lowrate_decoded.listSize){
			listSizes.push_back(0);
		}
		listSizes[lowrate_decoded.listSize - 1]++;
	}
	std::cout << "simulation complete- printing list sizes" << std::endl;
	print_int_vector(listSizes);
}

void lowrate_list_size_for_fixed_dist(codeInformation code){
	// run the list decoder on codewords with bits flipped to yield a distance metric of 168, and 
	// record the resulting list size. preliminarily, we expect this to be between 128 and 144
	std::cout << "working on rate-1/12 SLVD" << std::endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "TBCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	std::cout << "working on code " << filename.substr(0, filename.size() - 4) << std::endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	FeedforwardTrellis lowrate_trellis(code.k, code.n, code.v, code.numerators);
	std::cout << "trellis construction complete" << std::endl;

	int listSize = 1e5;
	int fixedDistance = 168;
	int numBitsToFlip = (fixedDistance / 4) - 1;
	std::cout << "bit flips: " << numBitsToFlip << std::endl;
	std::vector<int> listSizes;
	int maxNumTrials = 1e5;

	LowRateListDecoder LowRateListDecoder(lowrate_trellis, listSize, code.crcDeg, code.crc); 

	for(int i = 0; i < maxNumTrials; i++){
		std::vector<int> originalMessage;
		for (int i = 0; i < code.numInfoBits; i++){
			originalMessage.push_back(rand() % 2);
		}
		
		crc_calculation(originalMessage, code.crcDeg, code.crc);

		std::vector<int> encodedMessage = lowrate_trellis.encoder(originalMessage);

		std::vector<bool> bitFlipped;
		for(int i = 0; i < encodedMessage.size(); i++)
			bitFlipped.push_back(false);

		// puncture the 4 bits from length-516 noisy message
		std::vector<double> punc_msg;
		for (int i = 0; i < encodedMessage.size(); i++){
			if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
				punc_msg.push_back(encodedMessage[i]);
			}
			else{
				punc_msg.push_back(0);
				bitFlipped[i] = true; 		// we don't want to bit flip punctured bits
			}
		}

		int count = 0;
		while(count < numBitsToFlip){
			int index = rand()%encodedMessage.size();
			if(!bitFlipped[index]){
				bitFlipped[index] = true;
				punc_msg[index] *= -1;
				count++;
			}
		}
		/*
		int dist = 0;
		for(int i = 0; i < punc_msg.size(); i++){
			dist += pow(punc_msg[i] - encodedMessage[i], 2);
		}
		std::cout << "distance: " << dist << std::endl;
		*/
		LowRateListDecoder::messageInformation lowrate_decoded = LowRateListDecoder.lowRateDecoding(punc_msg);
		//std::cout << "list size: " << lowrate_decoded.listSize << std::endl;
		//std::cout << "metric: " << lowrate_decoded.metric << std::endl;
		if(lowrate_decoded.listSizeExceeded){
			std::cout << "list size exceeded- bug likely" << std::endl;
		}
		if(lowrate_decoded.message != originalMessage){
			std::cout << "our messages do NOT match- bug likely" << std::endl;
		}
		while(listSizes.size() < lowrate_decoded.listSize){
			listSizes.push_back(0);
		}
		listSizes[lowrate_decoded.listSize - 1]++;
	}
	std::cout << "simulation complete- printing list sizes" << std::endl;
	print_int_vector(listSizes);
}

void lowrate_test_zero_noise_codewords(codeInformation code){
	// run the list decoder on two different noiseless TB and CRC-checking codewords up to list size 2048  
	// record the 2048 sequences that are decoded.  
	// In particular, prepare the list of 2048 distances.  
	// Plot how distance grows as a function of list size.
	std::cout << "working on rate-1/12 SLVD" << std::endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "TBCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	std::cout << "working on code " << filename.substr(0, filename.size() - 4) << std::endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	FeedforwardTrellis lowrate_trellis(code.k, code.n, code.v, code.numerators);
	std::cout << "trellis construction complete" << std::endl;

	int listSize = 128;
	std::vector<int> listSizes;
	int maxNumTrials = 1e6;

	LowRateListDecoder LowRateListDecoder(lowrate_trellis, listSize, code.crcDeg, code.crc); 

	for(int i = 0; i < maxNumTrials; i++){
		std::vector<int> originalMessage;
		for (int i = 0; i < code.numInfoBits; i++){
			originalMessage.push_back(rand() % 2);
		}
		
		crc_calculation(originalMessage, code.crcDeg, code.crc);

		std::vector<int> encodedMessage = lowrate_trellis.encoder(originalMessage);		

		// puncture the 4 bits from length-516 noisy message
		std::vector<double> punc_msg;
		for (int i = 0; i < encodedMessage.size(); i++){
			if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
				punc_msg.push_back(encodedMessage[i]);
			}
			else{
				punc_msg.push_back(0);
			}
		}
		
		LowRateListDecoder::messageInformation lowrate_decoded = LowRateListDecoder.lowRateReciprocityCheck(punc_msg);
		//std::cout << "list size: " << lowrate_decoded.listSize << std::endl;
		if(lowrate_decoded.listSizeExceeded){
			std::cout << "list size exceeded- bug likely" << std::endl;
		}
		if(lowrate_decoded.message != originalMessage){
			std::cout << "our messages do NOT match- bug likely" << std::endl;
		}
		while(listSizes.size() < lowrate_decoded.listSize){
			listSizes.push_back(0);
		}
		listSizes[lowrate_decoded.listSize - 1]++;
	}
	std::cout << "simulation complete- printing list sizes" << std::endl;
	print_int_vector(listSizes);
	
}


void lowrate_test_variable_listsize(codeInformation code){
	std::cout << "working on rate-1/12 SLVD" << std::endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "TBCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	std::cout << "working on code " << filename.substr(0, filename.size() - 4) << std::endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	FeedforwardTrellis lowrate_trellis(code.k, code.n, code.v, code.numerators);
	std::cout << "trellis construction complete" << std::endl;

	vector<double> EbN0 = {4,5};			// will be running from 1 to 3.5 in increments of 0.5
	vector<double> SNR;  
	double offset = 10 * log10((double)2*32 / (double)512); // ?? should the rate be doubled? seems to agree better with jacob's results
	for (int i=0; i< EbN0.size(); i++){
		SNR.push_back(EbN0[i] + offset);
		std::cout << "SNR is: " << SNR[i] << std::endl;
	}

	for(int j = 0; j < SNR.size(); j++){
		int numErrors = 0;
		int maxNumTrials = 1e5;
		// int listSize = std::pow(2,18);					// limited to 2048- we already have the 1e5 data, no need to go above that
		int listSize = 1;					// limited to 2048- we already have the 1e5 data, no need to go above that
		int numTrials = 0;
		double listsize_sum = 0;
		std::vector<int> errorHistogram;
		std::vector<int> successHistogram;
		int numNACK = 0;
		std::cout << "beginning simulations at Eb/N0 = " << EbN0[j] << endl;

		LowRateListDecoder LowRateListDecoder(lowrate_trellis, listSize, code.crcDeg, code.crc); 

		// while(numTrials < maxNumTrials ){
		while(numTrials < maxNumTrials|| numErrors < 200){				//usable for running many parallel trials, but less useful if we want to hit a particular number of errors (a la 200)
			std::vector<int> originalMessage;
			for (int i = 0; i < code.numInfoBits; i++){
				originalMessage.push_back(rand() % 2);
			}

			crc_calculation(originalMessage, code.crcDeg, code.crc);

			std::vector<int> encodedMessage = lowrate_trellis.encoder(originalMessage);
			
			std::vector<double> noisyMessage = addNoise(encodedMessage, SNR[j]);

			// puncture the 4 bits from length-516 noisy message
			std::vector<double> punc_msg;
			for (int i = 0; i < noisyMessage.size(); i++){
				if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
					punc_msg.push_back(noisyMessage[i]);
				}
				else{
					punc_msg.push_back(0);
				}
			}
			LowRateListDecoder::messageInformation lowrate_decoded = LowRateListDecoder.lowRateDecoding(punc_msg);
			if(lowrate_decoded.listSizeExceeded){
				numNACK++;
				continue;
			}
			while(successHistogram.size() < lowrate_decoded.listSize){
				successHistogram.push_back(0);
				errorHistogram.push_back(0);
			}
			if(lowrate_decoded.message != originalMessage){
				numErrors++;
				errorHistogram[lowrate_decoded.listSize - 1]++;
			}
			else{
				successHistogram[lowrate_decoded.listSize - 1]++;
			}
			numTrials ++;
			
			
			if(numTrials%10000 == 0){
				std::cout << "we are at " << numTrials << " trials, with " << numErrors << " errors" << std::endl;
				std::cout << "and we have " << numNACK << " NACKs" << std::endl;
			}
			
		}
		std::cout << "Total num trials is " << numTrials << " with " << numErrors << " errors, FER is " << (double)numErrors/numTrials << std::endl;
		std::cout << "numNACK (where listSize > 2048): " << numNACK << std::endl;
		std::cout << "histogram of successful decodings: " << std::endl;
		print_int_vector(successHistogram);
		std::cout << "histogram of unsuccessful decodings: " << std::endl;
		print_int_vector(errorHistogram);
	}
	
}


void lowrate_test(codeInformation code){
	std::cout << "working on rate-1/12 SLVD" << std::endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "TBCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	std::cout << "working on code " << filename.substr(0, filename.size() - 4) << std::endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	FeedforwardTrellis lowrate_trellis(code.k, code.n, code.v, code.numerators);
	std::cout << "trellis construction complete" << std::endl;

	// vector<double> SNR = {-7};		// tune values for finding 10^-4 FER point
	vector<double> EbN0 = {11};			// will be running from 1 to 3.5 in increments of 0.5
	vector<double> SNR;  
	double offset = 10 * log10((double)2*32 / (double)512); // ?? should the rate be doubled? seems to agree better with jacob's results
	for (int i=0; i< EbN0.size(); i++){
		SNR.push_back(EbN0[i] + offset);
		std::cout << SNR[i] << std::endl;
	}
	
	for(int j = 0; j < SNR.size(); j++){
		int numErrors = 0;
		int maxNumTrials = 1e5;
		// int maxNumTrials = 1;
		int listSize = 1;
		int numTrials = 0;
		double listsize_sum = 0;
		// std::cout << "beginning simulations at Eb/N0 = " << EbN0[j] << endl;
		std::cout << "beginning simulations at SNR = " << SNR[j] << endl;

		LowRateListDecoder LowRateListDecoder(lowrate_trellis, listSize, code.crcDeg, code.crc); 

		while(numTrials < maxNumTrials || numErrors < 50){
		// while(numTrials < maxNumTrials){
			std::vector<int> originalMessage;
			for (int i = 0; i < code.numInfoBits; i++){
				originalMessage.push_back(rand() % 2);
			}
			// for(int i = 0; i < originalMessage.size(); i++)
			// 	std::cout << originalMessage[i];
			// std::cout << std::endl;

			crc_calculation(originalMessage, code.crcDeg, code.crc);

			std::vector<int> encodedMessage = lowrate_trellis.encoder(originalMessage);
			
			std::vector<double> noisyMessage = addNoise(encodedMessage, SNR[j]);

			// puncture the 4 bits from length-516 noisy message
			std::vector<double> punc_msg;
			for (int i = 0; i < noisyMessage.size(); i++){
				if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
					punc_msg.push_back(noisyMessage[i]);
				}
				else{
					punc_msg.push_back(0);
				}
			}
			/*
			for(int i = 0; i < punc_msg.size(); i++)
				std::cout << punc_msg[i];
			std::cout << std::endl;
			*/
			LowRateListDecoder::messageInformation lowrate_decoded = LowRateListDecoder.lowRateDecoding(punc_msg);
			// if(lowrate_decoded.listSize != 1)
			// 	std::cout << "non-one list size--------------------------------" << lowrate_decoded.listSize << std::endl;
			if(lowrate_decoded.message != originalMessage){
				// std::cout << "error" << std::endl;
				numErrors++;
			}
			numTrials ++;
			if(numTrials%10000 == 0){
				std::cout << "we are at " << numTrials << " trials, with " << numErrors << " errors" << std::endl;
			}		
		}
		std::cout << "Total num trials is " << numTrials << ", FER is " << (double)numErrors/numTrials << std::endl;
	}

	// // simulate the comms system
	// while(numTrials < 1e6){
	// 	numTrials++;
	// 	double snr = SNR[0];
		
	// 	std::vector<int> originalMessage;
	// 	for (int i = 0; i < code.numInfoBits; i++)
	// 		originalMessage.push_back(rand() % 2);
		
	// 	crc_calculation(originalMessage, code.crcDeg, code.crc);

	// 	std::vector<int> encodedMessage = lowrate_trellis.encoder(originalMessage);

	// 	// puncture the 4 bits from length-516 message
	// 	std::vector<int> punc_msg;
	// 	for (int i=0; i<encodedMessage.size(); i++){
	// 		if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
	// 			punc_msg.push_back(encodedMessage[i]);
	// 		}
	// 	}

	// 	std::vector<double> noisyMessage = addNoise(encodedMessage, snr);

	// 	LowRateListDecoder::messageInformation lowrate_decoded = LowRateListDecoder.lowRateDecoding(noisyMessage);
		
	// 	if(lowrate_decoded.message != originalMessage)
	// 		numErrors++;
	// 	if(numTrials % (int)1e4 == 0){
	// 		std::cout << "numTrials: " << numTrials << std::endl;
	// 		std::cout << "numErrors: " << numErrors << std::endl;
	// 		std::cout << "FER: " << (double)numErrors/numTrials << std::endl;
	// 	}
	// }
	
}

void testFERs(codeInformation code){
	std::cout << "testing the wava trellis for FER against the one trellis" << std::endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "TBCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	std::cout << "working on code " << filename.substr(0, filename.size() - 4) << std::endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators, code.denominator);
	std::cout << "feedback trellis construction complete" << std::endl;

	vector<double> SNR = {5};		// tune values for finding 10^-4 FER point

	DualTrellis decodingTrellis(code.hMatrix);
	std::cout << "dual trellis construction complete" << endl;

	int numWAVAErrors = 0;
	int numModifiedWAVAErrors = 0;
	int numOneTrellisErrors = 0;
	int maxNumTrials = 1e5;
	int listSize = 1e5;
	int numTrials = 0;
	double listsize_sum = 0;

	ListDecoder listDecoder(decodingTrellis, listSize, code.crcDeg, code.crc); 

	std::cout << "beginning simulations" << endl;
	// simulate the comms system
	while(numTrials < 1e6){
		numTrials++;
		double snr = SNR[0];
		
		std::vector<int> originalMessage;
		for (int i = 0; i < code.numInfoBits; i++)
			originalMessage.push_back(rand() % 2);
		
		crc_calculation(originalMessage, code.crcDeg, code.crc);

		std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);

		std::vector<double> noisyMessage = addNoise(encodedMessage, snr);
		// noisyMessage = {1.2854, -0.0264, -2.1992, -0.5423, -0.8308, -1.6942, -1.2302, -0.8181,  2.8997,  0.4703, -1.7166,  0.6112,  1.3851, -1.0335, -0.6206, -1.1088, -1.0659, -0.2091, -0.2520, -0.2476, -0.6435, -1.6410, -0.6192, -0.1345, -0.7405, -0.4507, -0.6141, -1.1611, -0.8440, -1.4180, -0.5284, -1.6090,  0.4326, -1.4298, -2.5631, -0.2364, -0.8274, -1.4008, -0.2725, -1.9086,  0.9457, -1.1282, -0.8305, -0.8339,  0.5408, -1.0160, -1.0875, -0.6668,  1.5804, -0.4111, -1.4585, -0.9589, -1.6446, -1.5911, -1.0036, -0.1864, -1.4086, -0.8028, -1.1198, -0.4068, -1.5782, -0.9827, -0.7067, -0.4157,  1.8198, -0.9544, -1.7919, -1.3941,  0.4364,  0.2478, -1.3268, -0.6029, -1.1022, -0.5283, -1.4060, -1.7444,  0.2449, -0.7408, -1.0942, -1.1041,  1.7535, -0.8452, -0.8950, -0.1571, -1.4271, -0.6302, -0.5567, -1.1294,  1.1145, -1.6189, -1.6094, -0.9443, -0.6166,  0.3726, -1.3540, -0.9005, -1.0438, -2.0262, -1.2330, -1.9528,  1.4461, -1.4714, -0.9469, -1.2891, -0.8389, -1.3187, -0.7399, -0.6075, -0.0912, -1.1031, -2.1352, -1.4457, -0.2809, -1.5692, -0.4898, -0.9341,  1.7627, -2.0410, -1.1050, -1.6412,  0.5438, -0.5619, -0.2679, -1.5618,  0.7512, -1.1446,  1.5831, -1.1475};

		ListDecoder::messageInformation wavaDecoding = listDecoder.wavaDecoding(noisyMessage);
		//std::cout << "before modified wava" << std::endl;
		ListDecoder::messageInformation modifiedWAVADecoding = listDecoder.modifiedWAVADecoding(noisyMessage);
		//std::cout << "after modified wava" << std::endl;
		ListDecoder::messageInformation oneTrellisDecoding = listDecoder.oneTrellisDecoding(noisyMessage);
		
		if(wavaDecoding.message != originalMessage)
			numWAVAErrors++;
		if(modifiedWAVADecoding.message != originalMessage)
			numModifiedWAVAErrors++;
		if(oneTrellisDecoding.message != originalMessage)
			numOneTrellisErrors++;
		if(numTrials % (int)1e4 == 0){
			std::cout << "numTrials: " << numTrials << std::endl;
			std::cout << "numWAVAErrors: " << numWAVAErrors << ", numModifiedWAVAErrors: " << numModifiedWAVAErrors << ", numOneTrellisErrors: " << numOneTrellisErrors << std::endl;
			std::cout << "wava FER: " << (double)numWAVAErrors/numTrials << ", modified WAVA FER: " << (double)numModifiedWAVAErrors/numTrials <<", oneTrellis FER: " << (double)numOneTrellisErrors/numTrials << std::endl;
		}
	}
	
}

void test(codeInformation code){
	std::cout << "used to generate wava points RCU bound gap vs complexity" << std::endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "ZTCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	std::cout << "working on code " << filename.substr(0, filename.size() - 4) << std::endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << std::endl;
		return;
	}

	FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators, code.denominator);
	std::cout << "feedback trellis construction complete" << std::endl;

	vector<double> SNR = {6, 6.5};		// tune values for finding 10^-4 FER point

	DualTrellis decodingTrellis(code.hMatrix);
	std::cout << "dual trellis construction complete" << endl;

	int numErrors = 0;
	int maxNumTrials = 1e5;
	int listSize = 1e5;
	int numTrials = 0;
	double listsize_sum = 0;


	ListDecoder listDecoder(decodingTrellis, listSize, code.crcDeg, code.crc); 

	std::cout << "beginning simulations" << endl;
	// simulate the comms system
	for(int i = 0; i < SNR.size(); i++){
		numErrors = 0;
		numTrials = 0;
		double snr = SNR[i];
		std::vector<int> listSizes;
		while(numTrials < maxNumTrials || numErrors < 200){
			
			std::vector<int> originalMessage;
			for (int i = 0; i < code.numInfoBits; i++)
				originalMessage.push_back(rand() % 2);
			
			crc_calculation(originalMessage, code.crcDeg, code.crc);

			std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);

			std::vector<double> noisyMessage = addNoise(encodedMessage, snr);

			ListDecoder::messageInformation wavaDecoding = listDecoder.modifiedWAVADecoding(noisyMessage);
			
			if(wavaDecoding.message != originalMessage)
				numErrors++;

			while(wavaDecoding.listSize >= listSizes.size())
				listSizes.push_back(0);
			listSizes[wavaDecoding.listSize]++;

			if(numTrials%(int)1e4 == 0)
				std::cout << "at trial " << numTrials << ", with num errors: " << numErrors << std::endl;
			numTrials++;
		}
		int nonZeroTrials = 0;
		for(int i = 1; i < listSizes.size(); i++)
			nonZeroTrials += listSizes[i];
		double meanNonZeroListSize = 0;
		for(int i = 1; i < listSizes.size(); i++)
			meanNonZeroListSize += (double)(listSizes[i]*i)/nonZeroTrials;
		std::cout << "------------------------------------------------------" << std::endl;
		std::cout << "at SNR: " << snr << std::endl;
		std::cout << "total number of trials: " << numTrials << std::endl;
		std::cout << "trials with list size zero: " << listSizes[0] << std::endl;
		std::cout << "mean list size for non-zero trials: " << meanNonZeroListSize << std::endl;
		std::cout << "------------------------------------------------------" << std::endl;

	}
	
}

// ZTCC simulations
void ztTrellisFERPerformance(codeInformation code){
	cout << "running 1-trellis ZTCC FER performance simulations" << endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "ZTCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	cout << "working on code " << filename.substr(0, filename.size() - 4) << endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		cout << "invalid msg + crc length" << endl;
		return;
	}

	std::ofstream outputFile;
	outputFile.open(filename, fstream::app);
	
	FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators, code.denominator);
	
	cout << "feedback trellis construction complete" << endl;

	vector<double> SNR = {3, 5};		// tune values for finding 10^-4 FER point

	DualTrellis decodingTrellis(code.hMatrix);
	// get dual trellis termination bits
	
	std::cout << "dual trellis construction complete" << endl;

	int numZTErrors = 0;
	int maxNumTrials = 1e5;
	int listSize = 1e5;
	int numTrials = 0;
	double listsize_sum = 0;

	ListDecoder listDecoder(decodingTrellis, listSize, code.crcDeg, code.crc); 

	std::cout << "beginning simulations" << endl;
	// simulate the comms system
	for (int i = 0; i < SNR.size(); i++) {
		std::cout << "beginning simulations at SNR " << SNR[i] << "____________________" << std::endl; 
		numTrials = 0;
		numZTErrors = 0;
		listsize_sum = 0;
		while(numTrials < maxNumTrials || numZTErrors < 200){
		// while(numTrials < 1){
			double snr = SNR[i];
			
			std::vector<int> originalMessage;
			for (int i = 0; i < code.numInfoBits; i++)
				// originalMessage.push_back(rand() % 2);
				originalMessage.push_back(1);
			
			crc_calculation(originalMessage, code.crcDeg, code.crc);
			// std::cout<< "CRC message: " << std::endl;
			// for (int i=0; i<originalMessage.size(); i++){
			// 	std::cout<< originalMessage[i];
			// }
			// std::cout<< "end of message" << std::endl;

			std::vector<int> terminatedMessage = encodingTrellis.terminateMsg(originalMessage);
			// std::cout<< "ZT message: " << std::endl;
			// for (int i=0; i<terminatedMessage.size(); i++){
			// 	std::cout<< terminatedMessage[i];
			// }
			// std::cout<< "end of message" << std::endl;


			std::vector<int> encodedMessage = encodingTrellis.ztencoder(terminatedMessage);
			// std::cout<< "encoded message: " << std::endl;
			// for (int i=0; i<encodedMessage.size(); i++){
			// 	std::cout<< encodedMessage[i];
			// }
			// std::cout<< "end of message" << std::endl;

			std::vector<double> noisyMessage = addNoise(encodedMessage, snr);
			// noisyMessage = {1.2854, -0.0264, -2.1992, -0.5423, -0.8308, -1.6942, -1.2302, -0.8181,  2.8997,  0.4703, -1.7166,  0.6112,  1.3851, -1.0335, -0.6206, -1.1088, -1.0659, -0.2091, -0.2520, -0.2476, -0.6435, -1.6410, -0.6192, -0.1345, -0.7405, -0.4507, -0.6141, -1.1611, -0.8440, -1.4180, -0.5284, -1.6090,  0.4326, -1.4298, -2.5631, -0.2364, -0.8274, -1.4008, -0.2725, -1.9086,  0.9457, -1.1282, -0.8305, -0.8339,  0.5408, -1.0160, -1.0875, -0.6668,  1.5804, -0.4111, -1.4585, -0.9589, -1.6446, -1.5911, -1.0036, -0.1864, -1.4086, -0.8028, -1.1198, -0.4068, -1.5782, -0.9827, -0.7067, -0.4157,  1.8198, -0.9544, -1.7919, -1.3941,  0.4364,  0.2478, -1.3268, -0.6029, -1.1022, -0.5283, -1.4060, -1.7444,  0.2449, -0.7408, -1.0942, -1.1041,  1.7535, -0.8452, -0.8950, -0.1571, -1.4271, -0.6302, -0.5567, -1.1294,  1.1145, -1.6189, -1.6094, -0.9443, -0.6166,  0.3726, -1.3540, -0.9005, -1.0438, -2.0262, -1.2330, -1.9528,  1.4461, -1.4714, -0.9469, -1.2891, -0.8389, -1.3187, -0.7399, -0.6075, -0.0912, -1.1031, -2.1352, -1.4457, -0.2809, -1.5692, -0.4898, -0.9341,  1.7627, -2.0410, -1.1050, -1.6412,  0.5438, -0.5619, -0.2679, -1.5618,  0.7512, -1.1446,  1.5831, -1.1475};

			ListDecoder::messageInformation ztDecoding = listDecoder.ztListDecoding(noisyMessage);
			
			// if(ztDecoding.message != originalMessage){
			if(ztDecoding.message != terminatedMessage){
				numZTErrors++;
			}
			listsize_sum += ztDecoding.listSize;
			numTrials++;
			if(numTrials%1000 == 0){
				std::cout << "currently at " << numTrials << " trials with " << numZTErrors << " errors with zt decoding" << std::endl;
			}
		}
		outputFile << "current SNR: " << SNR[i] << ". FER: " << (double)numZTErrors/(double)numTrials << ". list size: " << listsize_sum/(double)numTrials << std::endl; 
	}
		
	outputFile.close();
}

// TBCC simulations
void tbTrellisFERPerformance(codeInformation code){
	cout << "running 1-trellis TBCC FER performance simulations" << endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "TBCC rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	cout << "working on code " << filename.substr(0, filename.size() - 4) << endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		cout << "invalid msg + crc length" << endl;
		return;
	}

	std::ofstream outputFile;
	outputFile.open(filename, fstream::app);
	
	FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators, code.denominator);
	
	cout << "feedback trellis construction complete" << endl;

	vector<double> SNR = {3,3.5};		// tune values for finding 10^-4 FER point

	DualTrellis decodingTrellis(code.hMatrix);
	// get dual trellis termination bits
	
	std::cout << "dual trellis construction complete" << endl;

	int numTBErrors = 0;
	int maxNumTrials = 1e5;
	int listSize = 2048;
	int numTrials = 0;
	double listsize_sum = 0;

	ListDecoder rbListDecoder(decodingTrellis, listSize, code.crcDeg, code.crc); 

	std::cout << "beginning simulations" << endl;
	// simulate the comms system
	for (int i = 0; i < SNR.size(); i++) {
		std::cout << "beginning simulations at SNR " << SNR[i] << "____________________" << std::endl; 
		numTrials = 0;
		numTBErrors = 0;
		listsize_sum = 0;
		// while(numTrials < maxNumTrials || numTBErrors < 2){
		while(numTrials < maxNumTrials){
			double snr = SNR[i];
			
			std::vector<int> originalMessage;
			for (int i = 0; i < code.numInfoBits; i++)
				originalMessage.push_back(rand() % 2);
			
			crc_calculation(originalMessage, code.crcDeg, code.crc);


			std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);

			std::vector<double> noisyMessage = addNoise(encodedMessage, snr);

			// noisyMessage = {1.2854, -0.0264, -2.1992, -0.5423, -0.8308, -1.6942, -1.2302, -0.8181,  2.8997,  0.4703, -1.7166,  0.6112,  1.3851, -1.0335, -0.6206, -1.1088, -1.0659, -0.2091, -0.2520, -0.2476, -0.6435, -1.6410, -0.6192, -0.1345, -0.7405, -0.4507, -0.6141, -1.1611, -0.8440, -1.4180, -0.5284, -1.6090,  0.4326, -1.4298, -2.5631, -0.2364, -0.8274, -1.4008, -0.2725, -1.9086,  0.9457, -1.1282, -0.8305, -0.8339,  0.5408, -1.0160, -1.0875, -0.6668,  1.5804, -0.4111, -1.4585, -0.9589, -1.6446, -1.5911, -1.0036, -0.1864, -1.4086, -0.8028, -1.1198, -0.4068, -1.5782, -0.9827, -0.7067, -0.4157,  1.8198, -0.9544, -1.7919, -1.3941,  0.4364,  0.2478, -1.3268, -0.6029, -1.1022, -0.5283, -1.4060, -1.7444,  0.2449, -0.7408, -1.0942, -1.1041,  1.7535, -0.8452, -0.8950, -0.1571, -1.4271, -0.6302, -0.5567, -1.1294,  1.1145, -1.6189, -1.6094, -0.9443, -0.6166,  0.3726, -1.3540, -0.9005, -1.0438, -2.0262, -1.2330, -1.9528,  1.4461, -1.4714, -0.9469, -1.2891, -0.8389, -1.3187, -0.7399, -0.6075, -0.0912, -1.1031, -2.1352, -1.4457, -0.2809, -1.5692, -0.4898, -0.9341,  1.7627, -2.0410, -1.1050, -1.6412,  0.5438, -0.5619, -0.2679, -1.5618,  0.7512, -1.1446,  1.5831, -1.1475};
			
			ListDecoder::messageInformation decodedMessage = rbListDecoder.oneTrellisDecoding(noisyMessage);
			if(!decodedMessage.listSizeExceeded){
				if(decodedMessage.message != originalMessage){
					numTBErrors++;
				}
				listsize_sum += decodedMessage.listSize;
				numTrials++;
			}
			// else{
			// 	std::cout << "exceed list size" << std::endl;
			// }
			if(numTrials%10000 == 0){
				std::cout << "currently at " << numTrials << " trials with " << numTBErrors << " errors with TB decoding" << std::endl;
			}
		}
		outputFile << "current SNR: " << SNR[i] << ". FER: " << (double)numTBErrors/(double)numTrials << ". list size: " << listsize_sum/(double)numTrials << std::endl; 
	}
		
	outputFile.close();
}




void speedTest(){
	cout << "initializing" << endl;
	srand(time(NULL));
	string filename = "";
	
	int k = 3;				// numerator of the rate
	int n = 4;				// denominator of the rate
	int v = 5;				// number of memory elements
	int crcDeg = 7;			// degree of the crc
	int crc = 83;			// crc is given in decimal, convert from hex if necessary
	int infoLength = 42;	// number of information bits
	
	
	filename += "rate" + to_string(k) + "-" + to_string(n) + ",v" + to_string(v) + ",crcdeg" + to_string(crcDeg) + ",crc" + to_string(crc) + ",infolen" + to_string(infoLength) + ".txt";
	cout << "working on code " << filename.substr(0, filename.size() - 4) << endl;
	if ((infoLength + crcDeg - 1) % k != 0) {
		cout << "invalid msg + crc length" << endl;
		return;
	}

	// optimal numerator, denominator, and crc are known- this is one example
	/*std::vector<int> numerators = {33, 25, 37};
	int denominator = 31;*/
	
	std::vector<int> numerators = { 47, 73, 57 };
	int denominator = 75;
	
	
	FeedbackTrellis test(k, n, v, numerators, denominator);
	
	std::cout << "trellis construction complete" << std::endl;

	// Htest is {{denom},{numerators[last]},...,{numerators[0]}}, converted to binary then flipped lr

	//vector<vector<int>> Htest = { {1,0,0,1,1}, {1,1,1,1,1}, {1,0,1,0,1}, {1,1,0,1,1} };
	std::vector<std::vector<int>> Htest = { {1, 0, 1, 1, 1, 1}, {1, 1, 1, 1, 0, 1}, {1,1,0,1,1,1}, {1,1,1,0,0,1} };
	//vector<double> SNR = {2,2.5,3,3.5,4,4.5,5,5.5,6,6.5};
	//vector<double> SNR = {3,3.5,4,4.5,5};
	vector<double> SNR = {2,4,6,8};
	vector<int> numTrialsVector;
	vector<int> numErrorsVector;
	vector<double> timeTakenVector;


	DualTrellis Mtest(Htest);
	
	int numErrorsOriginal = 0;
	int numErrorsRB = 0;
	int maxNumTrials = 1e4;
	int listSize = 1e5;
	double encodingTime = 0;
	double decodingTime = 0;
	int numTrials = 0;

	ListDecoder rbListDecoder(Mtest, listSize, crcDeg, crc); 

	std::cout << "beginning simulations" << endl;
	// simulate the comms system
	for (int i = 0; i < SNR.size(); i++) {
		numTrials = 0;
		numErrorsOriginal = 0;
		numErrorsRB = 0;
		
		double snr = SNR[i];
		std::cout << "begin generating messages for snr " << snr << std::endl;
		std::vector<std::vector<double>> noisyEncodedMessages;
		std::vector<std::vector<int>> originalMessages;
		while(numTrials < maxNumTrials){
			std::vector<int> input;
			for(int i = 0; i < infoLength; i++)
				input.push_back(rand() % 2);
			crc_calculation(input, crcDeg, crc);
			std::vector<int> encodedMessage = test.encoder(input);
			std::vector<double> noisyEncodedMessage = addNoise(encodedMessage, snr);
			originalMessages.push_back(encodedMessage);
			noisyEncodedMessages.push_back(noisyEncodedMessage);
			numTrials++;
		}
		std::cout << "message generation finished" << std::endl;
		
		auto start = chrono::high_resolution_clock::now();
		numTrials = 0;
		while(numTrials < maxNumTrials){
			ListDecoder::messageInformation decodedMessage = rbListDecoder.nTrellisDecoding(noisyEncodedMessages[numTrials]);
			if ((numTrials+1) % 100 == 0)
				std::cout << "curr snr: " << snr << "\tnum trials : " << numTrials + 1 << std::endl;
			numTrials++;
		}
		auto end = chrono::high_resolution_clock::now();
		double timeTaken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
		timeTaken *= 1e-9;
		std::cout << "at snr " << snr << " the n trellis decoding took " << timeTaken << " secs" << "for an average of " << timeTaken/maxNumTrials << std::endl;
		
		start = chrono::high_resolution_clock::now();
		numTrials = 0;
		while(numTrials < maxNumTrials){
			ListDecoder::messageInformation decodedMessage = rbListDecoder.oneTrellisDecoding(noisyEncodedMessages[numTrials]);
			if ((numTrials+1) % 100 == 0)
				std::cout << "curr snr: " << snr << "\tnum trials : " << numTrials + 1 << std::endl;
			numTrials++;
		}
		end = chrono::high_resolution_clock::now();
		timeTaken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
		timeTaken *= 1e-9;
		std::cout << "at snr " << snr << " the one trellis decoding took " << timeTaken << " secs" << "for an average of " << timeTaken/maxNumTrials << std::endl;

	}
		
}

void wavaFERPerformance(codeInformation code){
	cout << "running wava FER performance simulations" << endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "WAVA rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	cout << "working on code " << filename.substr(0, filename.size() - 4) << endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		cout << "invalid msg + crc length" << endl;
		return;
	}

	std::ofstream outputFile;
	outputFile.open(filename, fstream::app);
	
	FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators, code.denominator);
	
	cout << "feedback trellis construction complete" << endl;

	vector<double> SNR = {5.5, 6, 6.5};		// tune values for finding 10^-4 FER point

	DualTrellis decodingTrellis(code.hMatrix);
	
	std::cout << "dual trellis construction complete" << endl;

	int numWAVAErrors = 0;
	int maxNumTrials = 1e5;
	int listSize = 1e5;
	int numTrials = 0;

	ListDecoder rbListDecoder(decodingTrellis, listSize, code.crcDeg, code.crc); 

	std::cout << "beginning simulations" << endl;
	// simulate the comms system
	for (int i = 0; i < SNR.size(); i++) {
		std::cout << "beginning simulations at SNR " << SNR[i] << "____________________" << std::endl; 
		numTrials = 0;
		numWAVAErrors = 0;
		while(numTrials < maxNumTrials || numWAVAErrors < 200){
			double snr = SNR[i];
			
			std::vector<int> originalMessage;
			for (int i = 0; i < code.numInfoBits; i++)
				originalMessage.push_back(rand() % 2);
			
			crc_calculation(originalMessage, code.crcDeg, code.crc);

			std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);

			std::vector<double> noisyMessage = addNoise(encodedMessage, snr);
			ListDecoder::messageInformation wavaDecoding = rbListDecoder.wavaDecoding(noisyMessage);
			
			if(wavaDecoding.message != originalMessage){
				numWAVAErrors++;
			}
			numTrials++;
			if(numTrials%10000 == 0){
				std::cout << "currently at " << numTrials << " trials with " << numWAVAErrors << " with wava decoding" << std::endl;
			}
		}
		outputFile << "current SNR: " << SNR[i] << ". FER: " << (double)numWAVAErrors/(double)numTrials << std::endl; 
	}
		
	outputFile.close();
}



void validateDecoder(codeInformation code){
	cout << "running decoder validation" << endl;
	srand(time(NULL));
	string filename = "";
	
	filename += "rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) + ".txt";
	cout << "working on code " << filename.substr(0, filename.size() - 4) << endl;
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		cout << "invalid msg + crc length" << endl;
		return;
	}
	
	FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators, code.denominator);
	
	cout << "feedback trellis construction complete" << endl;

	vector<double> SNR = {4,6,8};
	vector<int> numTrialsVector;
	vector<int> numErrorsVector;
	vector<double> timeTakenVector;

	DualTrellis decodingTrellis(code.hMatrix);
	
	std::cout << "dual trellis construction complete" << endl;

	int numMLErrors = 0;
	int numWAVAErrors = 0;
	int maxNumTrials = 1e5;
	int listSize = 1e5;
	double encodingTime = 0;
	double decodingTime = 0;
	int numTrials = 0;

	ListDecoder rbListDecoder(decodingTrellis, listSize, code.crcDeg, code.crc); 

	std::cout << "beginning simulations" << endl;
	// simulate the comms system
	for (int i = 0; i < SNR.size(); i++) {
		std::cout << "beginning simulations at " << SNR[i] << "____________________" << std::endl; 
		numTrials = 0;
		numMLErrors = 0;
		numWAVAErrors = 0;
		while(numTrials < maxNumTrials){
			double snr = SNR[i];
			
			std::vector<int> originalMessage;
			for (int i = 0; i < code.numInfoBits; i++)
				originalMessage.push_back(rand() % 2);
			
			crc_calculation(originalMessage, code.crcDeg, code.crc);

			std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);

			std::vector<double> noisyMessage = addNoise(encodedMessage, snr);
			//rblistdecoder::messageInformation nTrellisDecoding = rbListDecoder.nTrellisDecoding(noisyMessage);
			ListDecoder::messageInformation oneTrellisDecoding = rbListDecoder.oneTrellisDecoding(noisyMessage);
			ListDecoder::messageInformation wavaDecoding = rbListDecoder.wavaDecoding(noisyMessage);
			/*
			if(nTrellisDecoding.message != oneTrellisDecoding.message){
				std::cout << "the two ML approaches do not match" << std::endl;
				for(int i = 0; i < nTrellisDecoding.path.size(); i++){
					std::cout << nTrellisDecoding.path[i] << ", ";
				}
				std::cout << std::endl;
				for(int i = 0; i < oneTrellisDecoding.path.size(); i++){
					std::cout << oneTrellisDecoding.path[i] << ", ";
				}
				std::cout << std::endl;
			}
			if(oneTrellisDecoding.message != wavaDecoding.message){
				std::cout << "wava decoding disagrees (this is not necessarily a concern)" << std::endl;
			}*/
			if(oneTrellisDecoding.message != originalMessage){
				numMLErrors++;
			}
			if(wavaDecoding.message != originalMessage){
				numWAVAErrors++;
			}
			/*
			std::cout << "ntrellis list size: " << nTrellisDecoding.listSize << ", 1trellis list size: " << oneTrellisDecoding.listSize << std::endl;
			std::cout << "list size diff: " << oneTrellisDecoding.listSize - nTrellisDecoding.listSize << std::endl;*/
			if(numTrials%1000 == 0){
				std::cout << "currently at " << numTrials << " trials with " << numMLErrors << " errors with ML decoding, and " << numWAVAErrors << " with wava decoding" << std::endl;
			}
			numTrials++;
		}
	}
		
}


void listSizeHistograms(codeInformation code){
	std::cout << "Calling listSizeHistograms" << std::endl;
	string inputFileName = "pre_generated_messagesrate3-4,v4,crcdeg4,crc9,infolen93.txt";
	string outputFileName = "list_size_histograms_";
	
	std::ifstream inputFile;
	inputFile.open(inputFileName);

	string temp;
	std::getline(inputFile, temp);
	outputFileName += temp + ".txt";

	std::ofstream outputFile;
	outputFile.open(outputFileName);
	
	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		std::cout << "invalid msg + crc length" << endl;
		return;
	}
	FeedbackTrellis encodingTrellis(code.k, code.n, code.v, code.numerators, code.denominator);
	DualTrellis decodingTrellis(code.hMatrix);
	
	int numTrials = 1e6;
	int listSize = 1e5;
	int snr = 2;

	std::vector<int> oneTrellisListSizes;
	std::vector<int> wavaListSizes;
	std::vector<int> nTrellisListSizes;

	//outputFile << "decodingsPerTrial: " << decodingsPerTrial << std::endl;
	ListDecoder rbListDecoder(decodingTrellis, listSize, code.crcDeg, code.crc);

	std::cout << "beginning simulations" << endl;
	
	// simulate the comms system
	for (int i = 0; i < numTrials; i++) {
		std::vector<int> originalMessage;
		for (int i = 0; i < code.numInfoBits; i++)
			originalMessage.push_back(rand() % 2);
		
		crc_calculation(originalMessage, code.crcDeg, code.crc);

		std::vector<int> encodedMessage = encodingTrellis.encoder(originalMessage);

		std::vector<double> noisyMessage = addNoise(encodedMessage, snr);
		ListDecoder::messageInformation oneTrellisDecoding = rbListDecoder.oneTrellisDecoding(noisyMessage);
		ListDecoder::messageInformation wavaDecoding = rbListDecoder.wavaDecoding(noisyMessage);
		ListDecoder::messageInformation nTrellisDecoding = rbListDecoder.nTrellisDecoding(noisyMessage);

		if(oneTrellisDecoding.listSizeExceeded || wavaDecoding.listSizeExceeded || nTrellisDecoding.listSizeExceeded){
			std::cout << "list size exceeded issue" << std::endl;
			continue;
		}

		while(oneTrellisListSizes.size() < oneTrellisDecoding.listSize)
			oneTrellisListSizes.push_back(0);
		while(nTrellisListSizes.size() < nTrellisDecoding.listSize)
			nTrellisListSizes.push_back(0);
		while(wavaListSizes.size() < wavaDecoding.listSize)
			wavaListSizes.push_back(0);


		// we never have a list size of zero, so we offset indexing sucht that listSizes[0] 
		// is the number of trials with list size = 1
		oneTrellisListSizes[oneTrellisDecoding.listSize - 1]++;
		nTrellisListSizes[nTrellisDecoding.listSize - 1]++;
		wavaListSizes[wavaDecoding.listSize - 1]++;


		if(i%1000 == 0)
			std::cout << "at trial " << i << std::endl;
	}


	// outputing the results to the file
	outputFile << "one trellis list sizes" << std::endl;
	for(int i = 0; i < oneTrellisListSizes.size(); i++)
		outputFile << oneTrellisListSizes[i] << " ";
	outputFile << std::endl << "n trellis list sizes:" << std::endl;
	for(int i = 0; i < nTrellisListSizes.size(); i++)
		outputFile << nTrellisListSizes[i] << " ";
	outputFile << std::endl << "wava list sizes:" << std::endl;
	for(int i = 0; i < wavaListSizes.size(); i++)
		outputFile << wavaListSizes[i] << " ";


	inputFile.close();
	outputFile.close();
}

void generateMessages(codeInformation code){
	srand(time(NULL));
	ofstream myfile;
	string filename = "pre_generated_messages";
	string codeString = "rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits);
	filename += codeString + ".txt";

	std::vector<int> breakPoints;
	std::vector<bool> visitedPoints;
	for(int i = 0; i < 4; i++){
		for(int j = 1; j < 10; j++){
			breakPoints.push_back(j * std::pow(10, i));
			visitedPoints.push_back(false);
		}
	}
	for(int i = 0; i < breakPoints.size(); i++)
		std::cout << breakPoints[i] << ",";
	std::cout << std::endl;
	int count = 0;

	if ((code.numInfoBits + code.crcDeg - 1) % code.k != 0) {
		cout << "invalid msg + crc length" << endl;
		return;
	}

	myfile.open(filename, ofstream::app);
	myfile << "rate" + to_string(code.k) + "-" + to_string(code.n) + ",v" + to_string(code.v) + ",crcdeg" + to_string(code.crcDeg) + ",crc" + to_string(code.crc) + ",infolen" + to_string(code.numInfoBits) << std::endl;
	
	FeedbackTrellis encodingTrellis = FeedbackTrellis(code.k, code.n, code.v, code.numerators, code.denominator);
	
	vector<int> numTrialsVector;
	vector<int> numErrorsVector;
	vector<double> timeTakenVector;


	DualTrellis decodingTrellis(code.hMatrix);

	ListDecoder rbListDecoder(decodingTrellis, 1e5, code.crcDeg, code.crc); 

	int numErrors = 0;
	int maxNumTrials = 1e7;
	double encodingTime = 0;
	double decodingTime = 0;
	int numTrials = 0;
	cout << "beginning simulations" << endl;
	// simulate the comms system
	double snr = 3;
	numTrials = 0;
	while(numTrials < maxNumTrials){
		std::vector<int> input;
		for (int i = 0; i < code.numInfoBits; i++)
			input.push_back(rand() % 2);
		
		crc_calculation(input, code.crcDeg, code.crc);

		std::vector<int> encodedMessage = encodingTrellis.encoder(input);

		std::vector<double> noisyMessage = addNoise(encodedMessage, snr);
		
		ListDecoder::messageInformation oneTrellisDecoding = rbListDecoder.oneTrellisDecoding(noisyMessage);
		for(int i = 0; i < breakPoints.size() - 1; i++){
			if(!visitedPoints[i] && oneTrellisDecoding.listSize >= breakPoints[i] && oneTrellisDecoding.listSize < breakPoints[i + 1]){
				std::cout << "new breakpoint found: " << breakPoints[i] << std::endl;
				myfile << oneTrellisDecoding.listSize << ": " << std::endl;
				for(int i = 0; i < noisyMessage.size(); i++){
					myfile << noisyMessage[i] << ",";
				}
				myfile << std::endl;
				visitedPoints[i] = true;
				count++;
				break;
			}
		}
		if(count == breakPoints.size())
			break;
		numTrials++;
	}	
	myfile.close();
}