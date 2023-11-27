#include "lowratelistdecoder.h"
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <climits>
#include <fstream>
#include <iterator>
#include <sstream>



void LowRateListDecoder::generateNeighborList(std::vector<double> allZerosMessage, std::string dir, int listSize){
	// trellisInfo is indexed [startingState][state][stage]
	std::vector<std::vector<std::vector<cell>>> trellisInfo = constructLowRateMultiTrellis(allZerosMessage);
	std::cout << "list size: " << listSize << std::endl;

	struct neighbor{
		int weight;
		std::vector<int> codeword;
		std::vector<int> message;
		bool operator<(const neighbor& obj) const{
			return weight < obj.weight;
		}
		bool operator>(const neighbor& obj) const{
			return weight > obj.weight;
		}
	};

	// we first iterate over the ending state differences, then over the trellises, and within
	// each trellis, we perform list decoding
	for(int ESD = 0; ESD < lowrate_numStates; ESD++){
		// int ESD = 0;
		std::cout << "starting with an ESD of " << std::to_string(ESD) << std::endl;
		std::multiset<neighbor> candidateCodewords;
		for(int trellisIndex = 0; trellisIndex < lowrate_numStates; trellisIndex++){
			//std::cout << "trellis index: " << trellisIndex << std::endl;
			std::vector<std::vector<cell>> currentTrellis = trellisInfo[trellisIndex];

			minheap detourTree;
			std::vector<std::vector<int>> previousPaths;
			
			// we only want the appropriate ESD to be a valid ending state
			DetourObject detour;
			//detour.startingState = std::abs(trellisIndex - ESD);
			//detour.pathMetric = currentTrellis[std::abs(trellisIndex - ESD)][lowrate_pathLength - 1].pathMetric;
			detour.startingState = trellisIndex^ESD;
			detour.pathMetric = currentTrellis[trellisIndex^ESD][lowrate_pathLength - 1].pathMetric;
			detourTree.insert(detour);

			int numPathsSearched = 0;
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

					double suboptimalPathMetric = currentTrellis[currentState][newTracebackStage].suboptimalPathMetric;

					currentState = currentTrellis[currentState][newTracebackStage].suboptimalFatherState;
					newTracebackStage--;
					
					double prevPathMetric = currentTrellis[currentState][newTracebackStage].pathMetric;

					forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
					
				}
				path[newTracebackStage] = currentState;

				// actually tracing back
				for(int stage = newTracebackStage; stage > 0; stage--){
					double suboptimalPathMetric = currentTrellis[currentState][stage].suboptimalPathMetric;
					double currPathMetric = currentTrellis[currentState][stage].pathMetric;

					// if there is a detour we add to the detourTree
					if(currentTrellis[currentState][stage].suboptimalFatherState != -1){
						DetourObject localDetour;
						localDetour.detourStage = stage;
						localDetour.originalPathIndex = numPathsSearched;
						localDetour.pathMetric = suboptimalPathMetric + forwardPartialPathMetric;
						localDetour.forwardPathMetric = forwardPartialPathMetric;
						localDetour.startingState = detour.startingState;
						detourTree.insert(localDetour);
					}
					currentState = currentTrellis[currentState][stage].optimalFatherState;
					double prevPathMetric = currentTrellis[currentState][stage - 1].pathMetric;
					forwardPartialPathMetric += currPathMetric - prevPathMetric;
					path[stage - 1] = currentState;
				}
				previousPaths.push_back(path);

				std::vector<int> codeword = pathToCodeword(path);
				std::vector<int> message = pathToMessage(path);
				int metric = std::count(codeword.begin(), codeword.end(), -1);
				
				if(candidateCodewords.size() >= listSize && metric > (*(candidateCodewords.end()--)).weight)
					break;
				else{
					neighbor currNeighbor;
					currNeighbor.codeword = codeword;
					currNeighbor.message = message;
					currNeighbor.weight = metric;

					candidateCodewords.insert(currNeighbor);

					if(candidateCodewords.size() > listSize)
						candidateCodewords.erase(std::prev(candidateCodewords.end()));
				}

				// ///// only store neighbors that pass CRC
				// else if(crc_check(message, crcDegree, crc)){
				// 	neighbor currNeighbor;
				// 	currNeighbor.codeword = codeword;
				// 	currNeighbor.message = message;
				// 	currNeighbor.weight = metric;

				// 	candidateCodewords.insert(currNeighbor);

				// 	if(candidateCodewords.size() > listSize)
				// 		candidateCodewords.erase(std::prev(candidateCodewords.end()));
				// }
				/*
				pair currPair;
				currPair.codeword = codeword;
				currPair.weight = metric;
				candidateCodewords.insert(currPair);
				*/
				numPathsSearched++;
				
			}
		}
		std::ofstream outputFile;
		std::string codewordFileName = "c_" + dir + std::to_string(ESD) + ".txt";
		outputFile.open(codewordFileName);	//no appending
		int count = 0;
		for(auto it = candidateCodewords.begin(); it != candidateCodewords.end(); it++){
			std::vector<int> codeword = (*it).codeword;
			for(int i = 0; i < codeword.size(); i++)
				outputFile << std::to_string(codeword[i]) + " ";
			outputFile << std::endl;
			count++;
			if(count == listSize)
				break;
		}
		outputFile.close();
		std::string messageFileName = "m_" + dir + std::to_string(ESD) + ".txt";
		outputFile.open(messageFileName);	//no appending
		count = 0;
		for(auto it = candidateCodewords.begin(); it != candidateCodewords.end(); it++){
			std::vector<int> message = (*it).message;
			for(int i = 0; i < message.size(); i++)
				outputFile << std::to_string(message[i]) + " ";
			outputFile << std::endl;
			count++;
			if(count == listSize)
				break;
		}
		outputFile.close();
	}
}



void LowRateListDecoder::readNeighborList(std::string path){
    codewordNeighborList = std::vector<std::vector<std::vector<int>>>(lowrate_numStates);
	messageNeighborList = std::vector<std::vector<std::vector<int>>>(lowrate_numStates);
	std::string codewordPath = "c_" + path;
	std::string messagePath = "m_" + path;
    for(int ESD = 0; ESD < lowrate_numStates; ESD++){
        std::string codewordFileName = codewordPath + std::to_string(ESD) + ".txt";
		std::string messageFileName = messagePath + std::to_string(ESD) + ".txt";

        std::ifstream inputFile;
        inputFile.open(codewordFileName);

        std::string codewordLine;
        // reads the file in line by line
        while(std::getline(inputFile, codewordLine)){
            std::vector<int> codeword;
            std::stringstream stream(codewordLine);
            std::string val;
            while(stream >> val)
                codeword.push_back(std::stoi(val));
            codewordNeighborList[ESD].push_back(codeword);

        }

	    inputFile.close();

        inputFile.open(messageFileName);

        std::string messageLine;
        // reads the file in line by line
        while(std::getline(inputFile, messageLine)){
            std::vector<int> message;
            std::stringstream stream(messageLine);
            std::string val;
            while(stream >> val){
                message.push_back(std::stoi(val));
			}
            messageNeighborList[ESD].push_back(message);

        }
	    inputFile.close();
    }
}

LowRateListDecoder::messageInformation LowRateListDecoder::linearityDecoder(std::vector<double> receivedMessage){
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
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				}
			}

		}
	}

    // perform the single traceback
	messageInformation output;

    double minMetric = INT_MAX;
    int endingState = -1;

    for(int i = 0; i < lowrate_numStates; i++){
        if(trellisInfo[i][lowrate_pathLength - 1].pathMetric < minMetric){
            minMetric = trellisInfo[i][lowrate_pathLength - 1].pathMetric;
            endingState = i;
        }
    }

    std::vector<int> path(lowrate_pathLength);
    path[lowrate_pathLength - 1] = endingState;
    int currentState = endingState;

    // traceback
    for(int stage = lowrate_pathLength - 1; stage > 0; stage--){
        currentState = trellisInfo[currentState][stage].optimalFatherState;
        path[stage - 1] = currentState;
    }

    std::vector<int> codeword = pathToCodeword(path);
    std::vector<int> message = pathToMessage(path);

    int ESD = path[0]^path[lowrate_pathLength - 1];
	// offset distance between the received word and the decoding center
	for(int i = 0; i < receivedMessage.size(); i++){
		output.offset += std::pow(receivedMessage[i] - codeword[i], 2);
	}
	// std::cout << ESD << std::endl;

	std::vector<int> indices;
	std::vector<std::vector<int>> candidateMessages;

    if(ESD == 0 && crc_check(message, crcDegree, crc)){
    // if(ESD == 0){
        output.message = message;
        output.path = path;
        output.listSize = 1;
		double metric = 0;
		for(int i = 0; i < receivedMessage.size(); i++)
			metric += std::pow(receivedMessage[i] - codeword[i], 2);
		output.metric = metric;
        return output;
    }
	else{
		// std::cout << messageNeighborList[ESD].size() << std::endl;
		// collate the messages that pass the crc
		for(int i = 0; i < messageNeighborList[ESD].size(); i++){
			std::vector<int> currMessage;
			for(int j = 0; j < message.size(); j++){
				currMessage.push_back(message[j]^messageNeighborList[ESD][i][j]);
			}
			// print_int_vector(currMessage);
			if(crc_check(currMessage, crcDegree, crc)){
				indices.push_back(i);
				candidateMessages.push_back(currMessage);
			}
		}
	}

	// std::cout << candidateMessages.size() << std::endl;
	if(candidateMessages.size() == 0){
		output.listSizeExceeded = true;
		return output;
	}
	else if(candidateMessages.size() == 1){
		output.message = candidateMessages[0];
		output.listSize = indices[0];
		std::vector<int> currCodeword = codewordNeighborList[ESD][indices[0]];
		// "XOR"ing the neighbors with the received nearest word
		for(int j = 0; j < currCodeword.size(); j++)
			currCodeword[j] = currCodeword[j]*codeword[j];
		output.path = currCodeword;
		double metric = 0;
		for(int i = 0; i < receivedMessage.size(); i++)
			metric += std::pow(receivedMessage[i] - currCodeword[i], 2);
		output.metric = metric;
		return output;
	}
	else{
		double currBestMetric = INT_MAX;
		double currBestIndex = -1;
		for(int i = 0; i < indices.size(); i++){
			std::vector<int> currCodeword = codewordNeighborList[ESD][indices[i]];
			// "XOR"ing the neighbors with the received nearest word
			for(int j = 0; j < currCodeword.size(); j++)
				currCodeword[j] = currCodeword[j]*codeword[j];
			// std::cout << indices[i] << std::endl;
			// print_int_vector(currCodeword);
			double currMetric = 0;
			for(int j = 0; j < currCodeword.size(); j++){
				currMetric += std::pow(currCodeword[j] - receivedMessage[j], 2);
			}
			if(currMetric < currBestMetric){
				output.path = currCodeword;
				currBestIndex = i;
				currBestMetric = currMetric;
			}
		}

		double metric = 0;
		for(int i = 0; i < receivedMessage.size(); i++)
			metric += std::pow(receivedMessage[i] - output.path[i], 2);
		output.metric = metric;

		output.message = candidateMessages[currBestIndex];
		output.metric = currBestMetric;
		output.listSize = indices[currBestIndex];
		return output;
	}

}



LowRateListDecoder::messageInformation LowRateListDecoder::multiLinearityDecoder(std::vector<double> receivedMessage, int NumBall){
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
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].optimalFatherState = currentState;
				}
			}
		}
	}

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
	std::vector<std::vector<int>> codewords;
	std::vector<std::vector<int>> messages;
	std::vector<int> path(lowrate_pathLength);
	std::vector<int> ESDs;

	// examine the neighbors of  NumBall codewords
	while(numPathsSearched <  NumBall){
		std::vector<int> path(lowrate_pathLength);
		DetourObject detour = detourTree.pop();

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

		// store multiple codeword information
		std::vector<int> cur_codeword = pathToCodeword(path);
		std::vector<int> cur_message = pathToMessage(path);
		int cur_ESD = path[0]^path[lowrate_pathLength - 1];
		
		// if the nearest codeword is TB & CRC
		if (cur_ESD == 0 &&  crc_check(cur_message, crcDegree, crc)){
			output.message = cur_message;
			output.codeword = cur_codeword;
			output.metric = forwardPartialPathMetric;
			return output;
		}

		codewords.push_back(cur_codeword);
		messages.push_back(cur_message);
		ESDs.push_back(cur_ESD);

		numPathsSearched ++;
	}	

	// collate the messages that pass the crc
	std::vector<std::vector<int>> candidateMessages;
	std::vector<std::vector<int>> candidateCodewords;

	for (int e = 0; e < ESDs.size(); e++){
		int cur_ESD = ESDs[e];
		std::vector<int> cur_msg = messages[e];
		std::vector<int> cur_cwd = codewords[e];

		for(int i = 0; i < messageNeighborList[cur_ESD].size(); i++){
			std::vector<int> xorMsg;
			std::vector<int> xorCwd;
			for(int j = 0; j < cur_msg.size(); j++){
				xorMsg.push_back(cur_msg[j]^messageNeighborList[cur_ESD][i][j]);
			}
			if(crc_check(xorMsg, crcDegree, crc)){
				candidateMessages.push_back(xorMsg);
				
				// "XOR"ing the neighbors with the received nearest word
				for(int j = 0; j < cur_cwd.size(); j++)
					xorCwd.push_back(cur_cwd[j]*codewordNeighborList[cur_ESD][i][j]);
				candidateCodewords.push_back(xorCwd);
			}
		}
	}

	if(candidateMessages.size() == 0){
		output.listSizeExceeded = true;
		return output;
	}
	else if(candidateMessages.size() == 1){
		output.message = candidateMessages[0]; 
		output.codeword = candidateCodewords[0];
		double metric = 0;
		for(int i = 0; i < receivedMessage.size(); i++)
			metric += std::pow(receivedMessage[i] -candidateCodewords[0][i], 2);
		output.metric = metric;
		return output;
	}
	else{
		double currBestMetric = INT_MAX;
		double currBestIndex = -1;
		for(int i = 0; i < candidateMessages.size(); i++){
			std::vector<int> currCodeword = candidateCodewords[i];
			double currMetric = 0;
			for(int j = 0; j < currCodeword.size(); j++){
				currMetric += std::pow(receivedMessage[j] - currCodeword[j], 2);
			}
			if(currMetric < currBestMetric){
				currBestIndex = i;
				currBestMetric = currMetric;
			}
		}
		output.metric = currBestMetric;
		output.message = candidateMessages[currBestIndex];
		output.codeword = candidateCodewords[currBestIndex];		
		
		return output;
	}

}

void LowRateListDecoder::findMiniBallNeighbors(std::vector<double> ESDcodeword, std::vector<int> ESDmessage, std::string dir, int ESD, int num_neighbors){
	// trellisInfo is indexed [state][stage]
	std::vector<std::vector<cell>> trellisInfo;
	trellisInfo = constructLowRateTrellis(ESDcodeword);

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
	int neighbor_count = 0;
	int dmin;
	bool terminate = false;
	bool firstTBpath = true;

	std::ofstream outputFile1;
	std::string codewordFileName = "miniball_Nneighbors_" + std::to_string(num_neighbors) + "_c_" + dir + std::to_string(ESD) + ".txt";
	outputFile1.open(codewordFileName);	//no appending

	std::ofstream outputFile2;
	std::string messageFileName = "miniball_Nneighbors_" + std::to_string(num_neighbors) + "_m_" + dir + std::to_string(ESD) + ".txt";
	outputFile2.open(messageFileName);

	while(!terminate){
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

		// tb and dmin check
		if(path[0] == path[lowrate_pathLength - 1]){
			if (firstTBpath){
				dmin = forwardPartialPathMetric;
				firstTBpath = false;
				// std::cout<< dmin<<std::endl;
			}
			if (forwardPartialPathMetric == dmin){
				// print to text files
				for(int i = 0; i < codeword.size(); i++)
					outputFile1 << std::to_string(codeword[i]*(int)ESDcodeword[i]) + " ";
				outputFile1 << std::endl;
				for(int i = 0; i < message.size(); i++)
					outputFile2 << std::to_string(message[i]^ESDmessage[i]) + " ";
				outputFile2 << std::endl;
				std::cout << "at ESD of " << ESD << ", metric is: " << forwardPartialPathMetric << std::endl;
			}
			// termination condition
			if(forwardPartialPathMetric > dmin){
				neighbor_count ++;
				std::cout << "larger metric is: " << forwardPartialPathMetric << std::endl;
				if(neighbor_count == num_neighbors)
					terminate = true;
				dmin = 	forwardPartialPathMetric;	
			}
		}

		numPathsSearched++;
	}
	outputFile1.close();
	outputFile2.close();
}


void LowRateListDecoder::readMiniBallNeighborList(std::string path, int num_neighbors){
    miniCodewordNeighborList = std::vector<std::vector<std::vector<int>>>(lowrate_numStates);
	miniMessageNeighborList = std::vector<std::vector<std::vector<int>>>(lowrate_numStates);
	std::string codewordPath = "miniball_Nneighbors_" +std::to_string(num_neighbors)+ "_c_" + path;
	std::string messagePath = "miniball_Nneighbors_" +std::to_string(num_neighbors)+ "_m_" + path;
    for(int ESD = 0; ESD < lowrate_numStates; ESD++){
        std::string codewordFileName = codewordPath + std::to_string(ESD) + ".txt";
		std::string messageFileName = messagePath + std::to_string(ESD) + ".txt";

        std::ifstream inputFile;
        inputFile.open(codewordFileName);

        std::string codewordLine;
        // reads the file in line by line
        while(std::getline(inputFile, codewordLine)){
            std::vector<int> codeword;
            std::stringstream stream(codewordLine);
            std::string val;
            while(stream >> val)
                codeword.push_back(std::stoi(val));
            miniCodewordNeighborList[ESD].push_back(codeword);
        }

	    inputFile.close();

        inputFile.open(messageFileName);

        std::string messageLine;
        // reads the file in line by line
        while(std::getline(inputFile, messageLine)){
            std::vector<int> message;
            std::stringstream stream(messageLine);
            std::string val;
            while(stream >> val){
                message.push_back(std::stoi(val));
			}
            miniMessageNeighborList[ESD].push_back(message);

        }
	    inputFile.close();
    }
}

LowRateListDecoder::messageInformation LowRateListDecoder::miniBallDecoder(std::vector<double> receivedMessage, double threshold){
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

	// std::ofstream outputFile;
	// std::string ESDfilename = "miniball_ESD_K64v8m8.txt";
	// outputFile.open(ESDfilename, std::fstream::app);	//have appending

	// std::ofstream outputFile_cwd;
	// std::string Cwdfilename = "miniball_cwd_K64v8m8.txt";
	// outputFile_cwd.open(Cwdfilename, std::fstream::app);	//have appending

	// keep track of the TB CRC codeword w minimal metric, even if it fails threshold	
	std::vector<int> codewordMinMetirc;
	std::vector<int> messageMinMetirc;
	double minMetric = 1e5; // initialize to a large number
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

		int ESD = path[0]^path[lowrate_pathLength - 1];
		// outputFile << std::to_string(ESD) << std::endl;  // track ESD values

		// // record all TB codewords at this ESD
		// if (ESD == 0){
		// 	for(int j = 0; j < codeword.size(); j++){
		// 		outputFile_cwd << std::to_string(codeword[j]) + " ";
		// 	}
		// 	outputFile_cwd << std::endl;
		// 	outputFile_cwd << "/////////////" << std::endl;
		// }
		// else{
		// 	for (int i=0; i< miniCodewordNeighborList[ESD].size(); i++){
		// 		std::vector<int> currCodeword = miniCodewordNeighborList[ESD][i];
				// for(int j = 0; j < currCodeword.size(); j++){
				// 	outputFile_cwd << std::to_string(currCodeword[j]) + " ";
				// }
				// outputFile_cwd << std::endl;
		// 	}
		// 	outputFile_cwd << "/////////////" << std::endl;
		// }

		// miniball decoder: for each path, apply the smallest ball and see if we get a TB & CRC passing cwd
		// std::vector<int> indices;

		std::vector<int> candidateMetrics;
		std::vector<std::vector<int>> candidateMessages;
		std::vector<std::vector<int>> candidateCodewords;

		if(ESD == 0 && crc_check(message, crcDegree, crc)){
			// outputFile << "/////////////" << std::endl;
			// outputFile_cwd << "/////////////" << std::endl;
			// outputFile.close();
			// outputFile_cwd.close();

			// before returning output, compare current best w codeword_minmetric (if it exits/list rank > 0)
			if (numPathsSearched > 0){
				if (minMetric > forwardPartialPathMetric){
					output.message = message;
					output.codeword = codeword;
					output.listSize = numPathsSearched + 1;
					output.metric = forwardPartialPathMetric;
				}
				else{
					output.message = messageMinMetirc;
					output.codeword = codewordMinMetirc;
					output.listSize = numPathsSearched + 1;
					output.metric = minMetric;
				}
			}
			return output;
		}
		else{
			// collate the messages that pass the crc
			for(int i = 0; i < miniMessageNeighborList[ESD].size(); i++){
				std::vector<int> currMessage;
				for(int j = 0; j < message.size(); j++){
					currMessage.push_back(message[j]^miniMessageNeighborList[ESD][i][j]);
				}

				if(crc_check(currMessage, crcDegree, crc)){
					std::vector<int> currCodeword = miniCodewordNeighborList[ESD][i];
					// "XOR"ing the neighbors with the received nearest word
					for(int j = 0; j < currCodeword.size(); j++){
						currCodeword[j] = currCodeword[j]*codeword[j];
					}
					double currMetric = 0; // metric of the neighbor codeword
					for(int j = 0; j < currCodeword.size(); j++){
						currMetric += std::pow(currCodeword[j] - receivedMessage[j], 2);
					}
					
					double metric_diff =  currMetric - forwardPartialPathMetric;
					if ( metric_diff < threshold && 0 < metric_diff){
						// indices.push_back(i);
						candidateMessages.push_back(currMessage);
						candidateCodewords.push_back(currCodeword);
						candidateMetrics.push_back(currMetric);
					}
					// store the codeword & message with minimal metric, regardless of it not passing threshold
					else if (currMetric < minMetric){
						minMetric = currMetric;
						for(int j = 0; j < currCodeword.size(); j++){
							codewordMinMetirc.push_back(currCodeword[j]);
						}
						for(int j = 0; j < currMessage.size(); j++){
							messageMinMetirc.push_back(currMessage[j]);
						}
					}
					
				}
			}
		}

		if(candidateMessages.size() == 1){
			// before returning output, compare current best w codeword_minmetric (if it exits/list rank > 0)
			if (numPathsSearched > 0){
				if (minMetric > candidateMetrics[0]){
					output.message = candidateMessages[0];
					output.codeword = candidateCodewords[0];
					output.listSize = numPathsSearched + 1;
					output.metric = candidateMetrics[0];
				}
				else{
					output.message = messageMinMetirc;
					output.codeword = codewordMinMetirc;
					output.listSize = numPathsSearched + 1;
					output.metric = minMetric;
				}
			}
			return output;
			// outputFile << "/////////////" << std::endl;
			// outputFile_cwd << "/////////////" << std::endl;
			// outputFile.close();
			// outputFile_cwd.close();
			
		}
		else if (candidateMessages.size() > 1){
			// std::cout << "multiple nearest CRC-passing codewords" << std::endl;
			double currBestMetric = INT_MAX;
			double currBestIndex = -1;
			//  compare all the candidate metrics
			for(int i = 0; i < candidateMetrics.size(); i++){
				double currMetric = candidateMetrics[i];
				if(currMetric < currBestMetric){
					currBestIndex = i;
					currBestMetric = currMetric;
				}
			}

			// outputFile << "/////////////" << std::endl;
			// outputFile_cwd << "/////////////" << std::endl;
			// outputFile.close();
			// outputFile_cwd.close();
			// before returning output, compare current best w codeword_minmetric (if it exits/list rank > 0)
			if (numPathsSearched > 0){
				if (minMetric > currBestMetric){
					output.message = candidateMessages[currBestIndex];
					output.codeword = candidateCodewords[currBestIndex];
					output.metric = currBestMetric;
					output.listSize = numPathsSearched + 1;
				}
				else{
					output.message = messageMinMetirc;
					output.codeword = codewordMinMetirc;
					output.listSize = numPathsSearched + 1;
					output.metric = minMetric;
				}
			}
			return output;
		}
		numPathsSearched++;
		
	}
	// outputFile << "/////////////" << std::endl;
	// outputFile_cwd << "/////////////" << std::endl;
	// outputFile.close();
	// outputFile_cwd.close();
	output.listSizeExceeded = true;
	return output;
}


void LowRateListDecoder::miniBallMetrics(std::vector<int> originalCwd,std::vector<double> receivedMessage,std::vector<int> MLcwd){
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

	std::ofstream output_metrics;
	std::string filename = "miniball_K64v8m8_3metrics_EbNo1.txt";
	output_metrics.open(filename, std::fstream::app);	//have appending

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
		int ESD = path[0]^path[lowrate_pathLength - 1];

		// record the current codeword's metric
		output_metrics << "current codeword's metric: " << forwardPartialPathMetric << std::endl;

		// miniball decoder: for each path, apply the smallest ball and see if we get a TB & CRC passing cwd
		std::vector<int> indices;
		std::vector<std::vector<int>> candidateMessages;

		if(ESD == 0 && crc_check(message, crcDegree, crc)){
			std::cout << ESD << std::endl;
			if(codeword == originalCwd){
				output_metrics << "correct codeword found" << std::endl;
			}
			if (codeword == MLcwd){
				output_metrics << "ML codeword found" << std::endl;
				output_metrics << "/////////////" << std::endl;
				output_metrics.close();
				return;
			}
			else{
				output_metrics << "non-ML TB CRC codeword found, continue searching" << std::endl;
			}
		}
		else{
			// collate the messages that pass the crc
			for(int i = 0; i < miniMessageNeighborList[ESD].size(); i++){
				std::vector<int> currMessage;
				for(int j = 0; j < message.size(); j++){
					currMessage.push_back(message[j]^miniMessageNeighborList[ESD][i][j]);
				}
				if(crc_check(currMessage, crcDegree, crc)){
					indices.push_back(i);
					candidateMessages.push_back(currMessage);
				}
				
			}
		}

		if(candidateMessages.size() == 1){
			std::vector<int> currCodeword = miniCodewordNeighborList[ESD][indices[0]];
			// "XOR"ing the neighbors with the received nearest word
			for(int j = 0; j < currCodeword.size(); j++){
				currCodeword[j] = currCodeword[j]*codeword[j];
			}

			double currMetric = 0;
			for(int j = 0; j < currCodeword.size(); j++){
				currMetric += std::pow(currCodeword[j] - receivedMessage[j], 2);
			}

			output_metrics << "Only neighbor TB codeword metric is: " << currMetric << std::endl;
			if (currCodeword == originalCwd){
				output_metrics << "correct codeword found" << std::endl;
			}
			if (currCodeword == MLcwd){
				output_metrics << "ML codeword found" << std::endl;
				output_metrics << "/////////////" << std::endl;
				output_metrics.close();
				return;
			}
			else{
				output_metrics << "non-ML TB CRC codeword found, continue searching" << std::endl;
			}

		}
		else if (candidateMessages.size() > 1){
			std::cout << "multiple nearest CRC-passing codewords" << std::endl;
			double currBestMetric = INT_MAX;
			double currBestIndex = -1;
			std::vector<int> bestCodeword; 
			for(int i = 0; i < indices.size(); i++){
				std::vector<int> currCodeword = miniCodewordNeighborList[ESD][indices[i]];
				// "XOR"ing the neighbors with the received nearest word
				for(int j = 0; j < currCodeword.size(); j++)
					currCodeword[j] = currCodeword[j]*codeword[j];
				double currMetric = 0;
				for(int j = 0; j < currCodeword.size(); j++){
					currMetric += std::pow(currCodeword[j] - receivedMessage[j], 2);
				}
				if(currMetric < currBestMetric){
					currBestIndex = i;
					currBestMetric = currMetric;
					bestCodeword = currCodeword;
				}
			}

			output_metrics << "Closest neighbor TB codeword metric is: " << currBestMetric << std::endl;
			if(bestCodeword == originalCwd){
				output_metrics << "correct codeword found" << std::endl;
			}
			if (bestCodeword == MLcwd){
				output_metrics << "ML codeword found" << std::endl;
				output_metrics << "/////////////" << std::endl;
				output_metrics.close();
				return;
			}
			else{
				output_metrics << "non-ML TB CRC codeword found, continue searching" << std::endl;
			}
		}
		numPathsSearched++;
		
	}
	output_metrics << "list size exceeded" << std::endl;
	output_metrics << "/////////////" << std::endl;
	output_metrics.close();
	output.listSizeExceeded = true;
	return;
}
