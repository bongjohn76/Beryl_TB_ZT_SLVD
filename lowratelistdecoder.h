#ifndef LOWRATELISTDECODER_H
#define LOWRATELISTDECODER_H

#include "dualtrellis.h"
#include "MinHeap.h"
#include "maxheap.h"
#include "helper_functions.h"
#include "feedforwardtrellis.h"
#include "set"
#include <climits>

class LowRateListDecoder{
public:
	LowRateListDecoder(FeedforwardTrellis FT, int listSize, int crcDegree, int crc);
	LowRateListDecoder(FeedforwardTrellis FT, int listSize, int crcDegree, int crc, std::vector<std::vector<int>> neighboring_cwds, std::vector<std::vector<int>> neighboring_msgs, std::vector<std::vector<int>> path_ie_state);

	struct messageInformation{
		messageInformation(LowRateListDecoder * decoder){
			neighbor_cwds = std::vector<std::vector<int>>(decoder->listSize, std::vector<int>(516,0));
			neighbor_msgs = std::vector<std::vector<int>>(decoder->listSize, std::vector<int>(43,0));
			path_ie = std::vector<std::vector<int>>(decoder->listSize, std::vector<int>());
		};
		messageInformation(){};
		bool operator<(const messageInformation &y) const { return metric < y.metric; }
		std::vector<std::vector<int>> neighbor_cwds; // ${listSize} x 516 matrix
		std::vector<std::vector<int>> neighbor_msgs; // ${listSize} x 43 matrix
		std::vector<std::vector<int>> path_ie; // initial and ending states of path
		std::vector<int> message;
		std::vector<int> codeword;
		std::vector<int> path;
		int listSize;
		int TBListSize = -1;
		bool listSizeExceeded = false;
		double metric = -1;
		double offset = -1;
		//can potentially add more information as needed for debugging
	};

	messageInformation lowRateDecoding(std::vector<double> receivedMessage);
	messageInformation offsetLowRateDecoding(std::vector<double> receivedMessage);
	messageInformation linearityLowRateListSearching(std::vector<double> receivedMessage);
	messageInformation linearityLowRateDecoding(std::vector<double> receivedMessage);
	messageInformation TBGuaranteedListSearching(std::vector<double> receivedMessage);
	messageInformation TBGuaranteedLinearityDecoding(std::vector<double> receivedMessage);
	messageInformation lowRateDecoding_listSize(std::vector<double> receivedMessage);
	messageInformation lowRateDecoding_recip(std::vector<double> receivedMessage);
	messageInformation lowRateReciprocityCheck(std::vector<double> originalPoint, std::vector<double> distantPoint = {});
	std::vector<double> lowRateDistanceDistribution(std::vector<double> receivedMessage);
	std::vector<messageInformation> fullListSearching(std::vector<double> receivedMessage); // find all possible lists each with same path_ie condiitons
	std::vector<std::set<LowRateListDecoder::messageInformation>> multiTrellisFullListFinder(std::vector<double> receivedMessage);
	messageInformation fullListDecoding(std::vector<double> receivedMessage);
	void printRangeOfCodewords(std::vector<double> receivedMessage, int start, int end);
	bool MIListisCompleted(std::vector<messageInformation>fullNeighborMI, int v);
	messageInformation multiTrellisLowRateDecoding(std::vector<double> receivedMessage);

	void generateNeighborList(std::vector<double> allZerosMessage, std::string dir, int listSize);
	messageInformation linearityDecoder(std::vector<double> receivedMessage);
	void readNeighborList(std::string path);
	void readMiniBallNeighborList(std::string path, int num_neighbors);

	messageInformation minimumLikelihoodLowRateDecoding(std::vector<double> receivedMessage);

	messageInformation multiLinearityDecoder(std::vector<double> receivedMessage, int NumBall);
	void findMiniBallNeighbors(std::vector<double> ESDcodeword, std::vector<int> ESDmessage, std::string dir, int ESD, int num_neighbors);
	messageInformation miniBallDecoder(std::vector<double> receivedMessage, double threshold);
	void miniBallMetrics(std::vector<int> originalCwd,std::vector<double> receivedMessage,std::vector<int> MLcwd);

	// both neighbor lists are indexed [ESD][place in list][bits / symbols]
	std::vector<std::vector<std::vector<int>>> codewordNeighborList;
	std::vector<std::vector<std::vector<int>>> messageNeighborList;

	std::vector<std::vector<std::vector<int>>> miniCodewordNeighborList;
	std::vector<std::vector<std::vector<int>>> miniMessageNeighborList;


private:
	int numForwardPaths;
	int listSize;
	int crcDegree;
	int crc;
	int n;

	std::vector<std::vector<int>> lowrate_nextStates;
	std::vector<std::vector<int>> lowrate_outputs;
	std::vector<std::vector<int>> neighboring_cwds; // ${listSize} x 516 matrix
	std::vector<std::vector<int>> neighboring_msgs;  // ${listSize} x 43 matrix
	std::vector<std::vector<int>> path_ie_state;
	std::vector<messageInformation> fullNeighborMI;
	int lowrate_numStates;
	int lowrate_symbolLength;
	int lowrate_pathLength;

	// // both neighbor lists are indexed [ESD][place in list][bits / symbols]
	// std::vector<std::vector<std::vector<int>>> codewordNeighborList;
	// std::vector<std::vector<std::vector<int>>> messageNeighborList;

	struct cell {
		int optimalFatherState = -1;
		int suboptimalFatherState = -1;
		double pathMetric = INT_MAX;
		double suboptimalPathMetric = INT_MAX;
		bool init = false;
	};
    std::vector<int> pathToMessage(std::vector<int>); 
    std::vector<int> pathToCodeword(std::vector<int>); 
	std::vector<std::vector<cell>> constructLowRateTrellis(std::vector<double> receivedMessage);
	std::vector<std::vector<std::vector<cell>>> constructLowRateMultiTrellis(std::vector<double> receivedMessage);
	std::vector<std::vector<cell>> constructMinimumLikelihoodLowRateTrellis(std::vector<double> receivedMessage);
};


#endif