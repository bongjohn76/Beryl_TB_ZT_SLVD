#ifndef MAX_HEAP_H
#define MAX_HEAP_H

#include <vector>
#include "MinHeap.h"

class maxheap{
public:
    maxheap();
    void insert(DetourObject);
    DetourObject pop();
    int size();
private:
    std::vector<DetourObject> detourList;
    void reHeap(int index);
    int parentIndex(int index);
    int rightChildIndex(int index);
    int leftChildIndex(int index);
};


#endif