#include "lowratelistdecoder.h"
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <climits>
#include <fstream>
#include <iterator>

LowRateListDecoder::LowRateListDecoder(FeedforwardTrellis feedforwardTrellis, int listSize, int crcDegree, int crc) {
    this->lowrate_nextStates = feedforwardTrellis.getNextStates();
	this->lowrate_outputs = feedforwardTrellis.getOutputs();
	this->lowrate_numStates = feedforwardTrellis.getNumStates();
	this->lowrate_symbolLength = feedforwardTrellis.getN();
	this->numForwardPaths = lowrate_nextStates[0].size();
	// this->lowrate_pathLength = lowrate_nextStates;
    this->listSize = listSize;
    this->crcDegree = crcDegree;
    this->crc = crc;
	
	int v = feedforwardTrellis.getV();
}

LowRateListDecoder::LowRateListDecoder(FeedforwardTrellis feedforwardTrellis, int listSize, int crcDegree, int crc, std::vector<std::vector<int>> neighboring_cwds, std::vector<std::vector<int>> neighboring_msgs, std::vector<std::vector<int>> path_ie_state) {
	this->lowrate_nextStates = feedforwardTrellis.getNextStates();
	this->lowrate_outputs = feedforwardTrellis.getOutputs();
	this->lowrate_numStates = feedforwardTrellis.getNumStates();
	this->lowrate_symbolLength = feedforwardTrellis.getN();
	this->numForwardPaths = lowrate_nextStates[0].size();
	// this->lowrate_pathLength = lowrate_nextStates;
    this->listSize = listSize;
    this->crcDegree = crcDegree;
    this->crc = crc;
	this->neighboring_cwds = neighboring_cwds;
	this->neighboring_msgs = neighboring_msgs;
	this->path_ie_state = path_ie_state;
	
	int v = feedforwardTrellis.getV();
}

std::vector<std::vector<LowRateListDecoder::cell>> LowRateListDecoder::constructLowRateTrellis(std::vector<double> receivedMessage){
	std::vector<std::vector<cell>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	trellisInfo = std::vector<std::vector<cell>>(lowrate_numStates, std::vector<cell>(lowrate_pathLength));

	// initializes all the valid starting states
	for(int i = 0; i < lowrate_numStates; i++){
		trellisInfo[i][0].pathMetric = 0;
		trellisInfo[i][0].init = true;
	}
	

	
	// building the trellis
	for(int stage = 0; stage < lowrate_pathLength - 1; stage++){
		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// otherwise, we compute the relevent information
			for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
				// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
				// beyond indexing the forward path

				int nextState = lowrate_nextStates[currentState][forwardPathIndex];
				
				// if the nextState is invalid, we move on
				if(nextState < 0)
					continue;
				
				double branchMetric = 0;
				std::vector<int> output_point = get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
				
				for(int i = 0; i < lowrate_symbolLength; i++){
					branchMetric += std::pow(receivedMessage[lowrate_symbolLength * stage + i] - (double)output_point[i], 2);
					// branchMetric += std::abs(receivedMessage[lowrate_symbolLength * stage + i] - (double)output_point[i]);
				}
				double totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
				
				// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
				if(!trellisInfo[nextState][stage + 1].init){
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
					trellisInfo[nextState][stage + 1].init = true;
				}
				else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				}
				else{
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
				}
			}

		}
	}
	return trellisInfo;
}

LowRateListDecoder::messageInformation LowRateListDecoder::lowRateDecoding(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	messageInformation output;
	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int TBPathsSearched = 0;
	while(numPathsSearched < this->listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);
		
		// one trellis decoding requires both a tb and crc check
		if(path[0] == path[lowrate_pathLength - 1] && crc_check(message, crcDegree, crc)){
			output.message = message;
			output.codeword = codeword;
			output.path = path;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			output.TBListSize = TBPathsSearched + 1;
		 	return output;
		}

		numPathsSearched++;
		if(path[0] == path[lowrate_pathLength - 1])
			TBPathsSearched++;
	}
	output.listSizeExceeded = true;
	return output;
}

LowRateListDecoder::messageInformation LowRateListDecoder::offsetLowRateDecoding(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// first, we perform a single SSV pass, then use that codeword as
	// the input of the standard low rate decoder
	std::vector<int> path(lowrate_pathLength);

	int newTracebackStage = lowrate_pathLength - 1;
	double forwardPartialPathMetric = 0;
	int currentState = -1;

	double minMetric = INT_MAX;
	for(int i = 0; i < lowrate_numStates; i++){
		if(trellisInfo[i][lowrate_pathLength - 1].pathMetric < minMetric){
			currentState = i;
			minMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		}
	}

	path[newTracebackStage] = currentState;

	// actually tracing back
	for(int stage = newTracebackStage; stage > 0; stage--){
		double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
		double currPathMetric = trellisInfo[currentState][stage].pathMetric;

		currentState = trellisInfo[currentState][stage].optimalFatherState;
		double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
		forwardPartialPathMetric += currPathMetric - prevPathMetric;
		path[stage - 1] = currentState;
	}

	std::vector<int> offsetCodeword = pathToCodeword(path);
	std::vector<double> doubleOffset(offsetCodeword.begin(), offsetCodeword.end());
	
	return lowRateDecoding(doubleOffset);
}

LowRateListDecoder::messageInformation LowRateListDecoder::linearityLowRateListSearching(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	messageInformation output(this);

	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0; 
	while(numPathsSearched < (this->listSize)){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);
		
		for (int j=0; j<output.neighbor_cwds[0].size(); j++){
			output.neighbor_cwds[numPathsSearched][j] = codeword[j];
		}
		for (int k=0; k<output.neighbor_msgs[0].size(); k++){
			output.neighbor_msgs[numPathsSearched][k] = message[k];
		}
		// print_int_vector(output.neighbors[numPathsSearched]);
		
		output.path_ie[numPathsSearched].push_back(path[0]);
		output.path_ie[numPathsSearched].push_back(path[lowrate_pathLength-1]);
		output.message = message;
		output.path = path;
		output.listSize = numPathsSearched;
		output.metric = forwardPartialPathMetric;

		numPathsSearched++;
	}
	output.listSizeExceeded = true;
	return output;

}

LowRateListDecoder::messageInformation LowRateListDecoder::linearityLowRateDecoding(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	messageInformation output(this);

	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0; 
	while(numPathsSearched < (this->listSize)){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);
		

		if (path[0] == path[lowrate_pathLength - 1] && crc_check(message, crcDegree, crc)){
			//std::cout << "trivially correct - linear" << std::endl;
			output.message = message;
			output.path = path;
			output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			return output;
		}

		// TB and CRC check
		std::vector<std::vector<int>> passingCodewords;
		std::vector<std::vector<int>> passingMessages;
		std::vector<std::vector<int>> xor_neighbor_cwds(this->neighboring_cwds.size(), std::vector<int>(this->neighboring_cwds[0].size(),0));
		std::vector<std::vector<int>> xor_neighbor_msgs(this->neighboring_msgs.size(), std::vector<int>(this->neighboring_msgs[0].size(),0));
		std::vector<std::vector<int>> xor_path_ie(this->path_ie_state.size(), std::vector<int>(this->path_ie_state[0].size(),0));

		for (int i=0; i < this->path_ie_state.size(); i++){
			//std::cout << "neighbor message path_ie: " << this->path_ie_state[i][0] << " , " << this->path_ie_state[i][1] << std::endl;
			xor_path_ie[i][0] = this->path_ie_state[i][0] ^ path[0];
			xor_path_ie[i][1] = this->path_ie_state[i][1] ^ path[lowrate_pathLength - 1];
			//std::cout << "updated neighbor message path_ie: " << this->path_ie_state[i][0] << " , " << this->path_ie_state[i][1] << std::endl;
			if (xor_path_ie[i][0] == xor_path_ie[i][1]){
				//std::cout << "We are at: " << i << "th neighbor" << std::endl;
				// XOR the codewords
				for (int j=0; j < this->neighboring_cwds[0].size(); j++){
					xor_neighbor_cwds[i][j] = this->neighboring_cwds[i][j] ^ codeword[j];
				}
				for (int k=0; k < this->neighboring_msgs[0].size(); k++){
					xor_neighbor_msgs[i][k] = this->neighboring_msgs[i][k] ^ message[k];
				}
				// codeword to message
				// messageInformation decoding;
				// std::vector<double> temp(this->neighboring_cwds[i].size());
				// for (int k=0; k<temp.size(); k++){
				// 	temp[k] = static_cast<double>(this->neighboring_cwds[i][k]);
				// }
				// decoding = this->lowRateDecoding(temp);
				
				if (crc_check(xor_neighbor_msgs[i], crcDegree, crc)){
					//std::cout << "We are at: " << i << "th neighbor" << std::endl;
					passingCodewords.push_back(xor_neighbor_cwds[i]);
					passingMessages.push_back(xor_neighbor_msgs[i]);
				}
				
			}
			else{
				continue;
			}
		}

		// std::cout << "  ++++++++++++ Passing Messages" << std::endl;
		// for (int i=0; i < passingMessages.size(); i++){
		// 	print_int_vector(passingMessages[i]);
		// }
		// Euclidean distance with received codeword of all passing codewords
		std::vector<double> euclidean_distance(passingCodewords.size(), 0);
		for (int i=0; i < passingCodewords.size(); i++){
			for (int j=0; j < passingCodewords[i].size(); j++){
				euclidean_distance[i] += std::pow(receivedMessage[j] - (double)passingCodewords[i][j], 2);
			}
		}

		// Find the index of the smallest hamming distance and use that to find the best message
		auto smallestIt = std::min_element(euclidean_distance.begin(), euclidean_distance.end());
		int index = std::distance(euclidean_distance.begin(), smallestIt);

		// Assign output.message value
		if (passingMessages.size() != 0){
			output.message = passingMessages[index];
			output.listSize = -3;
			return output;
		}
		output.listSizeExceeded = true;
		output.listSize = -4;
		return output;

	}
	output.listSizeExceeded = true;
	return output;
}

LowRateListDecoder::messageInformation LowRateListDecoder::TBGuaranteedListSearching(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	messageInformation output(this);

	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int curListSize = 0;
	while(curListSize < (this->listSize)){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);

		// TB Guarantee neighbors finder
		if (path[0] == path[lowrate_pathLength - 1]){
			
			for (int j=0; j<output.neighbor_cwds[0].size(); j++){
				output.neighbor_cwds[curListSize][j] = codeword[j];
			}
			for (int k=0; k<output.neighbor_msgs[0].size(); k++){
				output.neighbor_msgs[curListSize][k] = message[k];
			}
			output.path_ie[curListSize].push_back(path[0]);
			output.path_ie[curListSize].push_back(path[lowrate_pathLength-1]);
			
			output.message = message;
			output.path = path;
			output.listSize = numPathsSearched;
			output.metric = forwardPartialPathMetric;

			curListSize++;
		}
		numPathsSearched++;
	}
	output.listSizeExceeded = true;
	return output;
}


LowRateListDecoder::messageInformation LowRateListDecoder::TBGuaranteedLinearityDecoding(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	messageInformation output(this);

	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0; 
	while(numPathsSearched < (this->listSize)){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);

		// denoise process

		// if the closest TB-pass codeword also passes crc, we are done!
		// if (path[0] == path[lowrate_pathLength - 1] && crc_check(message, crcDegree, crc)){
		// 	// std::cout << "trivially correct - linear" << std::endl;
		// 	output.message = message;
		// 	output.path = path;
		// 	output.listSize = numPathsSearched + 1;
		// 	output.metric = forwardPartialPathMetric;
		// 	return output;
		// }

		if (path[0] == path[lowrate_pathLength - 1]){
			//std::cout << "Fouuuuund Ya!" << std::endl;
			std::vector<std::vector<int>> passingCodewords;
			std::vector<std::vector<int>> passingMessages;
			std::vector<std::vector<int>> xor_neighbor_cwds(this->neighboring_cwds.size(), std::vector<int>(this->neighboring_cwds[0].size(),0));
			std::vector<std::vector<int>> xor_neighbor_msgs(this->neighboring_msgs.size(), std::vector<int>(this->neighboring_msgs[0].size(),0));

			for (int i=0; i < this->neighboring_cwds.size(); i++){
					// XOR the codewords
					for (int j=0; j < this->neighboring_cwds[0].size(); j++){
						xor_neighbor_cwds[i][j] = this->neighboring_cwds[i][j] ^ codeword[j];
					}
					for (int k=0; k < this->neighboring_msgs[0].size(); k++){
						xor_neighbor_msgs[i][k] = this->neighboring_msgs[i][k] ^ message[k];
					}
					
					if (crc_check(xor_neighbor_msgs[i], crcDegree, crc)){
						// std::cout << "We are at: " << i << "th neighbor" << std::endl;
						passingCodewords.push_back(xor_neighbor_cwds[i]);
						passingMessages.push_back(xor_neighbor_msgs[i]);
					}
					
			}
			std::vector<double> euclidean_distance(passingCodewords.size(), 0);
			for (int i=0; i < passingCodewords.size(); i++){
				for (int j=0; j < passingCodewords[i].size(); j++){
					euclidean_distance[i] += std::pow(receivedMessage[j] - (double)passingCodewords[i][j], 2);
				}
			}

			// Find the index of the smallest hamming distance and use that to find the best message
			auto smallestIt = std::min_element(euclidean_distance.begin(), euclidean_distance.end());
			int index = std::distance(euclidean_distance.begin(), smallestIt);

			// Assign output.message value
			// std::cout << "Passing message length: " << passingMessages.size() << std::endl;
			if (passingMessages.size() != 0){
				output.message = passingMessages[index];
				// std::cout << "decoding success!" << std::endl;
				output.listSize = -3;
				return output;
			}
			output.listSizeExceeded = true;
			output.listSize = -4;
			return output;

		}
		numPathsSearched++;

	}
	output.listSizeExceeded = true;
	return output;
}

std::vector<LowRateListDecoder::messageInformation> LowRateListDecoder::fullListSearching(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	// std::vector<messageInformation> fullNeighborMI;
	
	for (int i=0; i < 256; i++){
		messageInformation temp = messageInformation();
		fullNeighborMI.push_back(temp);
	}

	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0; 
	while(MIListisCompleted(fullNeighborMI, 8) == 0){

		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);

		// generate a list with start and end condition of 0
		int id = path[0] ^ path[lowrate_pathLength - 1];
		std::vector<int> ie_states = {path[0], path[lowrate_pathLength - 1]};

		// if (path[0] == path[lowrate_pathLength -1]){
		// 	print_int_vector(message);
		// 	print_int_vector(ie_states);
		// 	print_int_vector(codeword);
		// }

		fullNeighborMI[id].neighbor_cwds.push_back(codeword);
		fullNeighborMI[id].neighbor_msgs.push_back(message);
		fullNeighborMI[id].path_ie.push_back(ie_states);
		// // if fullNeighborMI[id] is full, we continue the loop
		// if (fullNeighborMI[id].neighbor_cwds.size() == this->listSize){
		// 	numPathsSearched++;
		// 	continue;
		// }
		// else{
		// 	// if fullNeighborMI[id] is not full, we push_back the cwd, msg, and path_ie
		// 	fullNeighborMI[id].neighbor_cwds.push_back(codeword);
		// 	fullNeighborMI[id].neighbor_msgs.push_back(message);
		// 	fullNeighborMI[id].path_ie.push_back(ie_states);
		// }
		
		numPathsSearched++;
	}
	for (int i=0; i < this->fullNeighborMI.size(); i++){
		std::cout << "For " << i << "th list, there are " << fullNeighborMI[i].neighbor_cwds.size() << " neighbors" << std::endl;
	}
	std::cout << "Using " << numPathsSearched << " paths" << std::endl;
	std::cout << "Full List search completed!" << std::endl;
	return fullNeighborMI;
}

std::vector<std::vector<std::vector<LowRateListDecoder::cell>>> LowRateListDecoder::constructLowRateMultiTrellis(std::vector<double> receivedMessage){
	// trellisInfo is indexed [TB state][state][stage]
	std::vector<std::vector<std::vector<cell>>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	
	trellisInfo = std::vector<std::vector<std::vector<cell>>>(lowrate_numStates, std::vector<std::vector<cell>>(lowrate_numStates, std::vector<cell>(lowrate_pathLength)));
	//(numStates/2, std::vector<std::vector<cell>>(numStates, std::vector<cell>(pathLength)))

	// initializes all the valid starting states
	for(int i = 0; i < lowrate_numStates; i++){
		trellisInfo[i][i][0].pathMetric = 0;
		trellisInfo[i][i][0].init = true;
	}
	

	
	// building the trellis
	for(int TBState = 0; TBState < lowrate_numStates; TBState++){
		for(int stage = 0; stage < lowrate_pathLength - 1; stage++){
			for(int currentState = 0; currentState < lowrate_numStates; currentState++){
				// if the state / stage is invalid, we move on
				if(!trellisInfo[TBState][currentState][stage].init)
					continue;

				// otherwise, we compute the relevent information
				for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
					// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
					// beyond indexing the forward path

					int nextState = lowrate_nextStates[currentState][forwardPathIndex];
					
					// if the nextState is invalid, we move on
					if(nextState < 0)
						continue;
					
					double branchMetric = 0;
					std::vector<int> output_point = get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
					
					for(int i = 0; i < lowrate_symbolLength; i++){
						branchMetric += std::pow(receivedMessage[lowrate_symbolLength * stage + i] - (double)output_point[i], 2);
						// branchMetric += std::abs(receivedMessage[lowrate_symbolLength * stage + i] - (double)output_point[i]);
					}
					double totalPathMetric = branchMetric + trellisInfo[TBState][currentState][stage].pathMetric;
					
					// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
					if(!trellisInfo[TBState][nextState][stage + 1].init){
						trellisInfo[TBState][nextState][stage + 1].pathMetric = totalPathMetric;
						trellisInfo[TBState][nextState][stage + 1].optimalFatherState = currentState;
						trellisInfo[TBState][nextState][stage + 1].init = true;
					}
					else if(trellisInfo[TBState][nextState][stage + 1].pathMetric > totalPathMetric){
						trellisInfo[TBState][nextState][stage + 1].suboptimalPathMetric = trellisInfo[TBState][nextState][stage + 1].pathMetric;
						trellisInfo[TBState][nextState][stage + 1].suboptimalFatherState = trellisInfo[TBState][nextState][stage + 1].optimalFatherState;
						trellisInfo[TBState][nextState][stage + 1].pathMetric = totalPathMetric;
						trellisInfo[TBState][nextState][stage + 1].optimalFatherState = currentState;
					}
					else{
						trellisInfo[TBState][nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
						trellisInfo[TBState][nextState][stage + 1].suboptimalFatherState = currentState;
					}
				}

			}
		}
	}
	return trellisInfo;
}

std::vector<std::set<LowRateListDecoder::messageInformation>> LowRateListDecoder::multiTrellisFullListFinder(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<std::vector<cell>>> trellisInfo;
	trellisInfo = constructLowRateMultiTrellis(receivedMessage);

	// start search
	messageInformation output;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	std::vector<std::set<messageInformation>> result; 

	//for each ESD 
	for(int ESD = 0; ESD < lowrate_numStates; ESD++)
	{
		int worstPathMetric = 0; 
		for(int TBState = 0; TBState < lowrate_numStates; TBState++)
		{
			//run list decoding for each trellis for a specific start and end state pair
			//add the paths to a list and choose the best 32 paths by metric
			// create node for our only possible end state
				int endState = TBState^ESD; 
				DetourObject detour;
				detour.startingState = endState;
				detour.pathMetric = trellisInfo[endState][endState][lowrate_pathLength - 1].pathMetric;
				detourTree.insert(detour);
			int numPathsSearched = 0;
			while(numPathsSearched < listSize){
				DetourObject detour = detourTree.pop();
				std::vector<int> path(lowrate_pathLength);

				int newTracebackStage = lowrate_pathLength - 1;
				double forwardPartialPathMetric = 0;
				int currentState = detour.startingState;
				int TBState = detour.startingState;

				// if we are taking a detour from a previous path, we skip backwards to the point where we take the
				// detour from the previous path
				if(detour.originalPathIndex != -1){
					forwardPartialPathMetric = detour.forwardPathMetric;
					newTracebackStage = detour.detourStage;

					// while we only need to copy the path from the detour to the end, this simplifies things,
					// and we'll write over the earlier data in any case
					path = previousPaths[detour.originalPathIndex];
					currentState = path[newTracebackStage];

					double suboptimalPathMetric = trellisInfo[TBState][currentState][newTracebackStage].suboptimalPathMetric;

					currentState = trellisInfo[TBState][currentState][newTracebackStage].suboptimalFatherState;
					newTracebackStage--;
					
					double prevPathMetric = trellisInfo[TBState][currentState][newTracebackStage].pathMetric;

					forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
					
				}
				path[newTracebackStage] = currentState;

				// actually tracing back
				for(int stage = newTracebackStage; stage > 0; stage--){
					double suboptimalPathMetric = trellisInfo[TBState][currentState][stage].suboptimalPathMetric;
					double currPathMetric = trellisInfo[TBState][currentState][stage].pathMetric;

					// if there is a detour we add to the detourTree
					if(trellisInfo[TBState][currentState][stage].suboptimalFatherState != -1){
						DetourObject localDetour;
						localDetour.detourStage = stage;
						localDetour.originalPathIndex = numPathsSearched;
						localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
						localDetour.forwardPathMetric = forwardPartialPathMetric;
						localDetour.startingState = detour.startingState;
						detourTree.insert(localDetour);
					}
					currentState = trellisInfo[TBState][currentState][stage].optimalFatherState;
					double prevPathMetric = trellisInfo[TBState][currentState][stage - 1].pathMetric;
					forwardPartialPathMetric += currPathMetric - prevPathMetric;
					path[stage - 1] = currentState;
				}
				previousPaths.push_back(path);
				
				if(TBState == 0) //first, populate our list with the zero state case
				{
					std::vector<int> message = pathToMessage(path);
					output.message = message;
					output.path = path;
					output.listSize = numPathsSearched + 1;
					output.metric = forwardPartialPathMetric;
					result[ESD].insert(output);
					worstPathMetric = forwardPartialPathMetric; 
				}
				else if(forwardPartialPathMetric >= worstPathMetric) //if the metric is worse than our worst metric, there is no way we get better, so skip this stage
				{
					continue;
				}
				else //replace our worst path with the new better path
				{
					std::vector<int> message = pathToMessage(path);
					output.message = message;
					output.path = path;
					output.listSize = numPathsSearched + 1;
					output.metric = forwardPartialPathMetric;
					result[ESD].insert(output);
					result[ESD].erase(result[ESD].end());
					std::set<messageInformation>::iterator i = result[ESD].end();
					messageInformation worstMessage = *i; 
					worstPathMetric = worstMessage.metric;
				}
				numPathsSearched++;
			}
		}
	}
	return result;
}

LowRateListDecoder::messageInformation LowRateListDecoder::fullListDecoding(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	messageInformation output(this);

	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int curListSize = 0;
	while(curListSize < (this->listSize)){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);

		// if the message passes CRC and path is tailbiting, we return the MI
		if(path[0] == path[lowrate_pathLength - 1] && crc_check(message, crcDegree, crc)){
			output.message = message;
			output.path = path;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
		 	return output;
		}
		// elseif we find the right "bucket" using the correct id
		int id = path[0] ^ path[lowrate_pathLength - 1];
		// std::cout << "Decoding using full list" << std::endl;
		// std::cout << "Path start with: " << path[0] << ". End with: " << path[lowrate_pathLength - 1] << std::endl;
		// std::cout << "checking list number " << id << " with" << fullNeighborMI[id].neighbor_cwds.size() << std::endl;
		// std::cout << "List look like: ";
		// print_int_vector(this->fullNeighborMI[id].path_ie[0]);
		// std::cout << std::endl;

		// this->fullNeighborMI[id].
		std::vector<std::vector<int>> passingCodewords;
		std::vector<std::vector<int>> passingMessages;
		std::vector<std::vector<int>> xor_neighbor_cwds(this->fullNeighborMI[id].neighbor_cwds.size(), std::vector<int>(this->fullNeighborMI[id].neighbor_cwds[0].size(),0));
		std::vector<std::vector<int>> xor_neighbor_msgs(this->fullNeighborMI[id].neighbor_msgs.size(), std::vector<int>(this->fullNeighborMI[id].neighbor_msgs[0].size(),0));

		for (int i=0; i < this->fullNeighborMI[id].neighbor_cwds.size(); i++){
				// XOR the codewords
				for (int j=0; j < this->fullNeighborMI[id].neighbor_cwds[0].size(); j++){
					xor_neighbor_cwds[i][j] = this->fullNeighborMI[id].neighbor_cwds[i][j] ^ codeword[j];
				}
				for (int k=0; k < this->fullNeighborMI[id].neighbor_msgs[0].size(); k++){
					xor_neighbor_msgs[i][k] = this->fullNeighborMI[id].neighbor_msgs[i][k] ^ message[k];
				}
				
				if (crc_check(xor_neighbor_msgs[i], crcDegree, crc)){
					// std::cout << "We are at: " << i << "th neighbor" << std::endl;
					passingCodewords.push_back(xor_neighbor_cwds[i]);
					passingMessages.push_back(xor_neighbor_msgs[i]);
				}
				
		}
		std::vector<double> euclidean_distance(passingCodewords.size(), 0);
		for (int i=0; i < passingCodewords.size(); i++){
			for (int j=0; j < passingCodewords[i].size(); j++){
				euclidean_distance[i] += std::pow(receivedMessage[j] - (double)passingCodewords[i][j], 2);
			}
		}

		// Find the index of the smallest hamming distance and use that to find the best message
		auto smallestIt = std::min_element(euclidean_distance.begin(), euclidean_distance.end());
		int index = std::distance(euclidean_distance.begin(), smallestIt);

		// Assign output.message value
		// std::cout << "Passing message length: " << passingMessages.size() << std::endl;
		if (passingMessages.size() != 0){
			output.message = passingMessages[index];
			// std::cout << "decoding success!" << std::endl;
			output.listSize = -3;
			return output;
		}
		output.listSizeExceeded = true;
		output.listSize = -4;
		return output;

		curListSize++;
	}
	output.listSizeExceeded = true;
	return output;
}

LowRateListDecoder::messageInformation LowRateListDecoder::lowRateDecoding_listSize(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	messageInformation output;
	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;

	std::ofstream outputFile;
	std::string filename = "codewords.txt";
	outputFile.open(filename);
	while(numPathsSearched < listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		//std::cout << "path metric at " << numPathsSearched << ": " << forwardPartialPathMetric << std::endl;

		std::vector<int> codeword = pathToCodeword(path);
		std::vector<double> justPastMidpointCodeword;

		for(int i = 0; i < codeword.size(); i++){
			justPastMidpointCodeword.push_back(0.49*1 + 0.51*codeword[i]);
		}

		justPastMidpointCodeword[47] = 0;
		justPastMidpointCodeword[60] = 0;
		justPastMidpointCodeword[129] = 0;
		justPastMidpointCodeword[504] = 0;
		
		// std::cout << "List size is: "<< numPathsSearched << ", Distance metric is: " << forwardPartialPathMetric << std::endl;
		// print_int_vector(codeword);

		for(int i = 0; i < justPastMidpointCodeword.size() - 1; i++){
			outputFile << justPastMidpointCodeword[i] << " ";
		}
		outputFile << justPastMidpointCodeword[justPastMidpointCodeword.size() - 1];
		outputFile << std::endl;
		
		numPathsSearched++;
	}
	outputFile.close();
	return output;
}


LowRateListDecoder::messageInformation LowRateListDecoder::lowRateDecoding_recip(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	messageInformation output;
	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	while(numPathsSearched < listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);
		
		// if we reach the all-zero codeword
		std::vector<int> original;
		for (int i = 0; i < codeword.size(); i++){
			if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
				original.push_back(1);
			}
			else{
				original.push_back(0);
			}
		}
		if(codeword == original){
			output.message = message;
			output.path = path;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
		 	return output;
		}

		numPathsSearched++;
	}
	// output.listSizeExceeded = true;
	return output;
}


LowRateListDecoder::messageInformation LowRateListDecoder::lowRateReciprocityCheck(std::vector<double> originalPoint, std::vector<double> distantPoint){
	
	int localListSize = listSize;
	bool reciprocalDecoding = false;

	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;

	// in this case, we are taking the distant point and decoding until we return to the original point
	if(distantPoint.size() > 0){
		trellisInfo = constructLowRateTrellis(distantPoint);
		localListSize = 1e5;
		reciprocalDecoding = true;
	}
	else{
		trellisInfo = constructLowRateTrellis(originalPoint);
	}

	// start search
	messageInformation output;
	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	while(numPathsSearched < localListSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);
		
		std::vector<int> codeword = pathToCodeword(path);
		std::vector<double> punc_msg;
		for (int i = 0; i < codeword.size(); i++){
			if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
				punc_msg.push_back(codeword[i]);
			}
			else{
				punc_msg.push_back(0);
			}
		}

		std::vector<int> message = pathToMessage(path);
		
		numPathsSearched++;

		if(reciprocalDecoding){
			if(originalPoint == punc_msg){
				output.message = message;
				output.listSize = numPathsSearched;
				return output;
			}
		}
		else{
			if (numPathsSearched == listSize){
				return lowRateReciprocityCheck(originalPoint, punc_msg);
			}
		}
	}
	output.listSizeExceeded = true;
	return output;
}


std::vector<double> LowRateListDecoder::lowRateDistanceDistribution(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	messageInformation output;
	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	std::vector<double> distances;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	while(numPathsSearched < listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		numPathsSearched++;
		
		distances.push_back(forwardPartialPathMetric);
		
	}
	return distances;
}

void LowRateListDecoder::printRangeOfCodewords(std::vector<double> receivedMessage, int start, int end){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(receivedMessage);

	// start search
	messageInformation output;
	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	while(numPathsSearched < listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		numPathsSearched++;

		
		if(start < numPathsSearched && end > numPathsSearched){
			std::cout << "on the " << numPathsSearched << "th path: " << std::endl;
			print_int_vector(pathToCodeword(path));
		}
	}
}

//Determines whether size of each list is >= listSize. There are 2^v lists in total.
bool LowRateListDecoder::MIListisCompleted(std::vector<messageInformation>fullNeighborMI, int v){

    for(int i = 0; i < pow(2,v); i++)
        if(fullNeighborMI.at(i).neighbor_cwds.size() < this->listSize)
            return false;
    return true;
}

// converts a path through the tb trellis to the binary message it corresponds with
std::vector<int> LowRateListDecoder::pathToMessage(std::vector<int> path){
	std::vector<int> message;
	for(int pathIndex = 0; pathIndex < path.size() - 1; pathIndex++){
		for(int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++){
			if(lowrate_nextStates[path[pathIndex]][forwardPath] == path[pathIndex + 1])
				message.push_back(forwardPath);
		}
	}
	return message;
}

// converts a path through the tb trellis to the BPSK it corresponds with
// currently does NOT puncture the codeword
std::vector<int> LowRateListDecoder::pathToCodeword(std::vector<int> path){
	std::vector<int> nopunc_codeword;
	for(int pathIndex = 0; pathIndex < path.size() - 1; pathIndex++){
		for(int forwardPath = 0; forwardPath < numForwardPaths; forwardPath++){
			if(lowrate_nextStates[path[pathIndex]][forwardPath] == path[pathIndex + 1]){
				std::vector<int> output_bin;
				dec_to_binary(lowrate_outputs[path[pathIndex]][forwardPath], output_bin, lowrate_symbolLength);
				for (int outbit = 0; outbit < lowrate_symbolLength; outbit ++){
					nopunc_codeword.push_back(-2 * output_bin[outbit] + 1);
				}
			}
		}
	}
	/*
	// puncture 4 bits
	std::vector<int> codeword;
	for (int i = 0; i < nopunc_codeword.size(); i++){
		if ((i != 47) && (i != 60) && (i != 129) && (i != 504)){
			codeword.push_back(nopunc_codeword[i]);
		}
		else{
			codeword.push_back(0);
		}
	}
	*/
	return nopunc_codeword;
}


std::vector<std::vector<LowRateListDecoder::cell>> LowRateListDecoder::constructMinimumLikelihoodLowRateTrellis(std::vector<double> receivedMessage){
	std::vector<std::vector<cell>> trellisInfo;
	lowrate_pathLength = (receivedMessage.size() / lowrate_symbolLength) + 1;

	trellisInfo = std::vector<std::vector<cell>>(lowrate_numStates, std::vector<cell>(lowrate_pathLength));

	// initializes all the valid starting states
	for(int i = 0; i < lowrate_numStates; i++){
		trellisInfo[i][0].pathMetric = 0;
		trellisInfo[i][0].init = true;
	}
	
	// building the trellis
	for(int stage = 0; stage < lowrate_pathLength - 1; stage++){
		for(int currentState = 0; currentState < lowrate_numStates; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// otherwise, we compute the relevent information
			for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
				// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
				// beyond indexing the forward path

				int nextState = lowrate_nextStates[currentState][forwardPathIndex];
				
				// if the nextState is invalid, we move on
				if(nextState < 0)
					continue;
				
				double branchMetric = 0;
				std::vector<int> output_point = get_point(lowrate_outputs[currentState][forwardPathIndex], lowrate_symbolLength);
				
				for(int i = 0; i < lowrate_symbolLength; i++){
					branchMetric += std::pow(receivedMessage[lowrate_symbolLength * stage + i] - (double)output_point[i], 2);
					// branchMetric += std::abs(receivedMessage[lowrate_symbolLength * stage + i] - (double)output_point[i]);
				}
				double totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
				
				// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
				if(!trellisInfo[nextState][stage + 1].init){
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
					trellisInfo[nextState][stage + 1].init = true;
				}
				else if(trellisInfo[nextState][stage + 1].pathMetric < totalPathMetric){
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = trellisInfo[nextState][stage + 1].pathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = trellisInfo[nextState][stage + 1].optimalFatherState;
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				}
				else{
					trellisInfo[nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].suboptimalFatherState = currentState;
				}
			}

		}
	}
	return trellisInfo;
}

LowRateListDecoder::messageInformation LowRateListDecoder::minimumLikelihoodLowRateDecoding(std::vector<double> receivedMessage){
	std::ofstream outputFile;
	std::string filename = "farthest_codewords.txt";
	outputFile.open(filename);
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructMinimumLikelihoodLowRateTrellis(receivedMessage);

	// start search
	messageInformation output;
	//RBTree detourTree;
	maxheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;
	int TBCRCPathsSearched = 0;
	while(TBCRCPathsSearched < this->listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		bool isPresent = std::find(previousPaths.begin(), previousPaths.end(), path) != previousPaths.end();

		if (!isPresent){
			previousPaths.push_back(path);
		
			std::vector<int> message = pathToMessage(path);
			std::vector<int> codeword = pathToCodeword(path);

			if (path[0] == path[lowrate_pathLength - 1] && crc_check(message, crcDegree, crc) ){
				for(int i = 0; i < codeword.size() - 1; i++){
					outputFile << codeword[i] << " ";
				}
				outputFile << codeword[codeword.size() - 1] << std::endl;
				// outputFile << forwardPartialPathMetric << std::endl;
				TBCRCPathsSearched++;
				// std::cout << "metric is: " << forwardPartialPathMetric << std::endl;
			}
			
			numPathsSearched++;
		}
		else{
			std::cout << "repeated path" << std::endl;
		}
		
		// // one trellis decoding requires both a tb and crc check
		// if(path[0] == path[lowrate_pathLength - 1] && crc_check(message, crcDegree, crc)){
		// 	output.message = message;
		// 	output.codeword = codeword;
		// 	output.path = path;
		//  	output.listSize = numPathsSearched + 1;
		// 	output.metric = forwardPartialPathMetric;
		// 	output.TBListSize = TBPathsSearched + 1;
		//  	return output;
		// }

		// numPathsSearched++;

		// if(path[0] == path[lowrate_pathLength - 1])
		// 	TBPathsSearched++;
	}
	output.listSizeExceeded = true;
	outputFile.close();
	return output;
}

LowRateListDecoder::messageInformation LowRateListDecoder::multiTrellisLowRateDecoding(std::vector<double> receivedMessage){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<std::vector<cell>>> trellisInfo;
	trellisInfo = constructLowRateMultiTrellis(receivedMessage);

	// start search
	messageInformation output;
	//RBTree detourTree;
	minheap detourTree;
	std::vector<std::vector<int>> previousPaths;
	

	// create nodes for each valid ending state with no detours
	// std::cout<< "end path metrics:" <<std::endl;
	for(int i = 0; i < lowrate_numStates; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][i][lowrate_pathLength - 1].pathMetric;
		detourTree.insert(detour);
	}

	// std::ofstream outputFile_cwd;
	// std::string Cwdfilename = "multiTrellis_cwd_K64v8m8.txt";
	// outputFile_cwd.open(Cwdfilename, std::fstream::app);	//have appending

	int numPathsSearched = 0;
	while(numPathsSearched < listSize){
		DetourObject detour = detourTree.pop();
		std::vector<int> path(lowrate_pathLength);

		int newTracebackStage = lowrate_pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.startingState;
		int TBState = detour.startingState;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.originalPathIndex != -1){
			forwardPartialPathMetric = detour.forwardPathMetric;
			newTracebackStage = detour.detourStage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.originalPathIndex];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[TBState][currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[TBState][currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[TBState][currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[TBState][currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[TBState][currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[TBState][currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[TBState][currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[TBState][currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);
		std::vector<int> codeword = pathToCodeword(path);
		// for(int j = 0; j < codeword.size(); j++){
		// 	outputFile_cwd << std::to_string(codeword[j]) + " ";
		// }
		// outputFile_cwd << std::endl;
		
		// one trellis decoding requires both a tb and crc check
		if(crc_check(message, crcDegree, crc)){
			output.message = message;
			output.path = path;
		 	output.listSize = numPathsSearched + 1;
			output.metric = forwardPartialPathMetric;
			// outputFile_cwd << "/////////////" << std::endl;
			// outputFile_cwd.close();
		 	return output;
		}

		numPathsSearched++;
	}
	output.listSizeExceeded = true;
	// outputFile_cwd << "/////////////" << std::endl;
	// outputFile_cwd.close();
	return output;
}