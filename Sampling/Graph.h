#ifndef GRAPH_H_
#define GRAPH_H_
#include "Motif.h"
#include <map>

class Graph{
public:
	int NodeSize;
	int EdgeSize;
    vector<dynamic_bitset<> > adjacencyList;
    int **adjacencyMatrix;
    vector<double> edgeProbabilities;
    map<int,string> edgeNameMap;

public:
	Graph(int nodeNum, char* fileName): NodeSize(nodeNum), adjacencyList(nodeNum, dynamic_bitset<>(nodeNum)){
        adjacencyMatrix = new int*[nodeNum];
        for(int i=0;i<nodeNum;i++)
            adjacencyMatrix[i] = new int[nodeNum]();
        readNetworkFile(fileName);
	}
    //create a subgraph from the original graph
    Graph(dynamic_bitset<>& subgraph, int nodeNum, Graph* g0): NodeSize(nodeNum), adjacencyList(nodeNum, dynamic_bitset<>(nodeNum)){
        adjacencyMatrix = new int*[nodeNum];
        for(int i=0;i<nodeNum;i++)
            adjacencyMatrix[i] = new int[nodeNum]();
        int* mapNode = new int[nodeNum]();
        int j=0;
        for(int i=0;i<subgraph.size();i++){
            if(subgraph[i]==1){
                 mapNode[j]=i;
                 j++;
            }
                
        }
    
        j = 0;
        for(int i=0;i<nodeNum-1;i++){
            int iIndex = mapNode[i];
            for(int k=i+1;k<nodeNum;k++){
                int index = mapNode[k];
                if(g0->adjacencyList[iIndex][index]==1){
                    adjacencyList[i][k] = 1;
                    adjacencyList[k][i] = 1;
                    adjacencyMatrix[i][k] = j;
                    adjacencyMatrix[k][i] = j;
                    int edgeIndex = g0->adjacencyMatrix[iIndex][index];
                    edgeProbabilities.push_back(g0->edgeProbabilities[edgeIndex]);
                    edgeNameMap.insert(pair<int,string>(j,g0->edgeNameMap[edgeIndex]));
                    j++;
                }
            }
        }
        delete[] mapNode;

        EdgeSize = j;
    }

	~Graph(){
        for(int i=0;i<NodeSize;i++)
        	delete[] adjacencyMatrix[i];
        delete[] adjacencyMatrix;
        edgeProbabilities.clear();
        adjacencyList.clear();
        edgeNameMap.clear();
	}

	void findMotifs(vector<Motif*>& pattern, int caseN);

private:
	void readNetworkFile(char* fileName);	
};

#endif