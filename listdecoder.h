#ifndef LISTDECODER_H
#define LISTDECODER_H

#include "dualtrellis.h"
#include "MinHeap.h"
#include "helper_functions.h"
// #include "feedbacktrellis.h"

#include <climits>

class ListDecoder{
public:
	ListDecoder(DualTrellis DT, int listSize, int crcDegree, int crc);

	struct messageInformation{
		messageInformation(): listSizeExceeded(false) {};
		std::vector<int> message;
		std::vector<int> path;
		int listSize;
		bool listSizeExceeded;
		//can potentially add more information as needed for debugging
	};

	messageInformation nTrellisDecoding(std::vector<double> receivedMessage);

	messageInformation oneTrellisDecoding(std::vector<double> receivedMessage);

	messageInformation wavaDecoding(std::vector<double> receivedMessage);
	messageInformation modifiedWAVADecoding(std::vector<double> receivedMessage);

	messageInformation ztListDecoding(std::vector<double> receivedMessage);

private:
	std::vector<std::vector<std::vector<int> > > nextStates;
	int numTrellisSegLength;
	int numStates;
	int numForwardPaths;
	int pathLength;
	int listSize;
	int crcDegree;
	int crc;
	int ztTrailingBits;
	int n;

	struct cell {
		cell(): optimalFatherState(-1), suboptimalFatherState(-1), 
				pathMetric(INT_MAX), suboptimalPathMetric(INT_MAX), init(false) {};
		int optimalFatherState;
		int suboptimalFatherState;
		double pathMetric;
		double suboptimalPathMetric;
		bool init;
	};

	std::vector<int> pathToMessage(std::vector<int>); 
	std::vector<int> ztPathToMessage(std::vector<int> path);

	std::vector<std::vector<std::vector<cell>>> constructNTrellis(std::vector<double> receivedMessage);

	std::vector<std::vector<cell>> constructOneTrellis(std::vector<double> receivedMessage);

	std::vector<std::vector<cell>> constructWAVATrellis(std::vector<double> receivedMessage, std::vector<std::vector<cell>> oneTrellis);

	std::vector<std::vector<cell>> constructZTTrellis(std::vector<double> receivedMessage);

	
};


#endif