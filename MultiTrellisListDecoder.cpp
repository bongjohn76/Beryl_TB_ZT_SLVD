#include "listdecoder.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <climits>

std::vector<std::vector<std::vector<ListDecoder::cell>>> ListDecoder::constructNTrellis(std::vector<double> receivedMessage){
	// trellisInfo is indexed [starting/ending state][state][stage]
    std::vector<std::vector<std::vector<cell>>> trellisInfo;
	trellisInfo = std::vector<std::vector<std::vector<cell>>>(numStates/2, std::vector<std::vector<cell>>(numStates, std::vector<cell>(pathLength)));

	// initializes all the valid starting states
	for(int i = 0; i < numStates/2; i++){
		trellisInfo[i][i][0].pathMetric = 0;
		trellisInfo[i][i][0].init = true;
	}

	// precomputing euclidean distance between the received signal and +/- 1
	std::vector<std::vector<double>> precomputedMetrics;
	precomputedMetrics = std::vector<std::vector<double>>(receivedMessage.size(), std::vector<double>(2));

	for(int stage = 0; stage < receivedMessage.size(); stage++){
		// precomputedMetrics[stage][0] = std::abs(receivedMessage[stage] - 1);
		// precomputedMetrics[stage][1] = std::abs(receivedMessage[stage] + 1);
		precomputedMetrics[stage][0] = std::pow(receivedMessage[stage] - 1,2);
		precomputedMetrics[stage][1] = std::pow(receivedMessage[stage] + 1,2);
	}

	// building the trellis
	for(int startingState = 0; startingState < numStates/2; startingState++){
		for(int stage = 0; stage < receivedMessage.size(); stage++){
			for(int currentState = 0; currentState < numStates; currentState++){
				// if the state / stage is invalid, we move on
				if(!trellisInfo[startingState][currentState][stage].init)
					continue;

				// otherwise, we compute the relevent information
				for(int forwardPathIndex = 0; forwardPathIndex < numForwardPaths; forwardPathIndex++){
					// note that the forwardPathIndex is also the bit that corresponds with the 
					// trellis transition

					int nextState = nextStates[currentState][stage%numTrellisSegLength][forwardPathIndex];

					// if the nextState is invalid, we move on
					if(nextState < 0)
						continue;

					double totalPathMetric = precomputedMetrics[stage][forwardPathIndex] + trellisInfo[startingState][currentState][stage].pathMetric;

					// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
					if(!trellisInfo[startingState][nextState][stage + 1].init){
						trellisInfo[startingState][nextState][stage + 1].pathMetric = totalPathMetric;
						trellisInfo[startingState][nextState][stage + 1].optimalFatherState = currentState;
						trellisInfo[startingState][nextState][stage + 1].init = true;
					}
					else if(trellisInfo[startingState][nextState][stage + 1].pathMetric > totalPathMetric){
						trellisInfo[startingState][nextState][stage + 1].suboptimalPathMetric = trellisInfo[startingState][nextState][stage + 1].pathMetric;
						trellisInfo[startingState][nextState][stage + 1].suboptimalFatherState = trellisInfo[startingState][nextState][stage + 1].optimalFatherState;
						trellisInfo[startingState][nextState][stage + 1].pathMetric = totalPathMetric;
						trellisInfo[startingState][nextState][stage + 1].optimalFatherState = currentState;
					}
					else{
						trellisInfo[startingState][nextState][stage + 1].suboptimalPathMetric = totalPathMetric;
						trellisInfo[startingState][nextState][stage + 1].suboptimalFatherState = currentState;
					}
				}

			}
		}
	}
	return trellisInfo;
}

ListDecoder::messageInformation ListDecoder::nTrellisDecoding(std::vector<double> receivedMessage){
	
	pathLength = receivedMessage.size() + 1;

	// trellisInfo is indexed [starting/ending state][state][stage]
    std::vector<std::vector<std::vector<cell>>> trellisInfo;
	trellisInfo = constructNTrellis(receivedMessage);

	// start search
	messageInformation output;
	//RBTree detourTree;
	minheap detourTree;

	std::vector<std::vector<int>> previousPaths;

	// create nodes for each valid ending state with no detours
	for(int i = 0; i < numStates / 2; i++){
		DetourObject detour;
		detour.startingState = i;
		detour.pathMetric = trellisInfo[i][i][pathLength - 1].pathMetric;
		//queue.push(detour);
		detourTree.insert(detour);
	}

	int numPathsSearched = 0;

	while(numPathsSearched < listSize){
		//DetourObject detour = queue.top();
		//queue.pop();
		DetourObject detour = detourTree.pop();
		std::vector<int> path(pathLength);

		int newTracebackStage = pathLength - 1;
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

			double suboptimalPathMetric = trellisInfo[detour.startingState][currentState][newTracebackStage].suboptimalPathMetric;

			currentState = trellisInfo[detour.startingState][currentState][newTracebackStage].suboptimalFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[detour.startingState][currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[detour.startingState][currentState][stage].suboptimalPathMetric;
			double currPathMetric = trellisInfo[detour.startingState][currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[detour.startingState][currentState][stage].suboptimalFatherState != -1){
				DetourObject localDetour;
				localDetour.detourStage = stage;
				localDetour.originalPathIndex = numPathsSearched;
				localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forwardPathMetric = forwardPartialPathMetric;
				localDetour.startingState = detour.startingState;
				//queue.push(localDetour);
				detourTree.insert(localDetour);
			}
			currentState = trellisInfo[detour.startingState][currentState][stage].optimalFatherState;
			double prevPathMetric = trellisInfo[detour.startingState][currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> message = pathToMessage(path);

		if(crc_check(message, crcDegree, crc)){
			output.message = message;
			output.path = path;
			output.listSize = numPathsSearched + 1;
			return output;
		}


		numPathsSearched++;
	}
	output.listSizeExceeded = true;
	return output;
}