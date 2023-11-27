#include "maxheap.h"

maxheap::maxheap(){
    // constructor if necessary
}

int maxheap::parentIndex(int index){return (index - 1)/2;}
// int maxheap::leftChildIndex(int index){return (2 * index + 1);}
// int maxheap::rightChildIndex(int index){return (2 * index + 2);}
int maxheap::leftChildIndex(int index){return (2 * index + 2);}
int maxheap::rightChildIndex(int index){return (2 * index + 1);}

void maxheap::insert(DetourObject detour){
    detourList.push_back(detour);
    int index = detourList.size() - 1;
    while(index > 0 && detourList[parentIndex(index)] < detourList[index]){
        std::swap(detourList[parentIndex(index)], detourList[index]);
        index = parentIndex(index);
    }
}

DetourObject maxheap::pop(){
    DetourObject detour = detourList[0];
    detourList[0] = detourList[detourList.size() - 1];
    detourList.pop_back();
    reHeap(0);

    return detour;
}

void maxheap::reHeap(int index){
    int leftIndex = leftChildIndex(index);
    int rightIndex = rightChildIndex(index);
    int maxDetourIndex = index;
    if(leftIndex < detourList.size() && detourList[leftIndex] > detourList[maxDetourIndex])
        maxDetourIndex = leftIndex;
    if(rightIndex < detourList.size() && detourList[rightIndex] > detourList[maxDetourIndex])
        maxDetourIndex = rightIndex;
    
    if(maxDetourIndex != index){
        std::swap(detourList[index], detourList[maxDetourIndex]);
        reHeap(maxDetourIndex);
    }
}

int maxheap::size(){
    return detourList.size();
}