#include "Graph.h"
#include <fstream>
#include <boost/algorithm/string.hpp>

void Graph::readNetworkFile(char* fileName, double threshold){
	int sourceNode, destNode, mappingCounter = 0;
	ifstream in_stream;
	string sourceName, destName, edgeProb,str1;	
    int edgeIndex=0;
    in_stream.open(fileName);
	vector<string> strs;
	map<string,int> nodeNameToNoMap;

	while(!in_stream.eof()) {
        getline(in_stream, str1);
	    if(str1.size()>0){
	        strs.clear();
	        boost::split(strs,str1,boost::is_any_of("\t"));
	        sourceName=strs[0];
	        destName=strs[1];
	        edgeProb=strs[2];
	        if(edgeProb.compare("-")==0)	            	
                continue;	            
			 
			if(nodeNameToNoMap.find(sourceName) == nodeNameToNoMap.end()) {
				nodeNameToNoMap.insert(pair<string,int>(sourceName, mappingCounter));
				++mappingCounter;
			}
			    
			if(nodeNameToNoMap.find(destName) == nodeNameToNoMap.end()) {
				nodeNameToNoMap.insert(pair<string,int>(destName, mappingCounter));
				++mappingCounter;
			}

			sourceNode = nodeNameToNoMap.find(sourceName)->second;
			destNode = nodeNameToNoMap.find(destName)->second;

			double prob=strtof(edgeProb.c_str(), 0);

			if(adjacencyList[sourceNode][destNode]==0 && prob > threshold){
				adjacencyMatrix[sourceNode][destNode] = edgeIndex;
				adjacencyMatrix[destNode][sourceNode] = edgeIndex;
				adjacencyList[sourceNode][destNode] = 1;
				adjacencyList[destNode][sourceNode] = 1;
				edgeProbabilities.push_back(prob);
				edgeNameMap.insert(pair<int,string>(edgeIndex,sourceName+" "+destName));
				edgeIndex++;
			}
		}
	}

	in_stream.close();
    EdgeSize=edgeIndex;
    strs.clear();
    nodeNameToNoMap.clear();
}

void Graph::findMotifs(vector<Motif*>& pattern, int caseN){
    dynamic_bitset<>::size_type it,itt,iit;

	int index1,index2,index3;
	index1=index2=index3=-1;

	int count = 0;
	if(caseN==1){
	    for(int i=0;i<NodeSize;i++){
		    it=adjacencyList[i].find_first();
		    while(it!=dynamic_bitset<>::npos){
		        itt=adjacencyList[i].find_next(it);
		        while(itt!=dynamic_bitset<>::npos){
		    	    //pattern1		    	
		            Motif* t=new Motif(count);
		            index1=adjacencyMatrix[i][it];
		            index2=adjacencyMatrix[i][itt];		           
		            t->edges.push_back(index1);
		            t->edges.push_back(index2);	
		            t->A = edgeProbabilities[index1] * edgeProbabilities[index2];            	
		            pattern.push_back(t);
		            count++;		    	
		    	    itt=adjacencyList[i].find_next(itt);
		        }
		        it=adjacencyList[i].find_next(it);
		    }		
	    }
	}

	else if(caseN==2){
	    for(int i=0;i<NodeSize;i++){
		    it=adjacencyList[i].find_next(i);
		    while(it!=dynamic_bitset<>::npos){
		        itt=adjacencyList[i].find_next(it);
		        while(itt!=dynamic_bitset<>::npos){
		    		if(adjacencyList[it][itt]==1){
		        		Motif* t=new Motif(count);
		        		index1=adjacencyMatrix[i][it];
		        		index2=adjacencyMatrix[i][itt];
		        		index3=adjacencyMatrix[it][itt];
		                t->edges.push_back(index1);
		                t->edges.push_back(index2);
						t->edges.push_back(index3);
						t->A = edgeProbabilities[index1] * edgeProbabilities[index2] * edgeProbabilities[index3];
		        		pattern.push_back(t);
                        count++;		        		
		    	    }
		    	    itt=adjacencyList[i].find_next(itt);
		        }
		        it=adjacencyList[i].find_next(it);
		    }		
	    }
	}

	else if(caseN==3){
	    for(int i=0;i<NodeSize;i++){
		    it=adjacencyList[i].find_first();
		    while(it!=dynamic_bitset<>::npos){
		        itt=adjacencyList[i].find_next(it);
		        while(itt!=dynamic_bitset<>::npos){
		    		if(adjacencyList[it].count()>1){
		    			iit=adjacencyList[it].find_next(itt);
		    			while(iit!=dynamic_bitset<>::npos){
		    				if(iit!=i){
		    					Motif* t=new Motif(count);
		    					index1=adjacencyMatrix[i][it];
		    					index2=adjacencyMatrix[i][itt];
		    					index3=adjacencyMatrix[it][iit];
		                        t->edges.push_back(index1);
		                        t->edges.push_back(index2);
						        t->edges.push_back(index3);
						        t->A = edgeProbabilities[index1] * edgeProbabilities[index2] * edgeProbabilities[index3];
                                count++;
		        		        pattern.push_back(t);		  
		    				}		    				
		    				iit=adjacencyList[it].find_next(iit);
		    			}
		    		}
		    			
		    		if(adjacencyList[itt].count()>1){
		    			iit=adjacencyList[itt].find_next(it);
		    			while(iit!=dynamic_bitset<>::npos){
		    				if(iit!=i){
		    					Motif* t=new Motif(count);
		    					index1=adjacencyMatrix[i][it];
		    					index2=adjacencyMatrix[i][itt];
		    					index3=adjacencyMatrix[itt][iit];
		                        t->edges.push_back(index1);
		                        t->edges.push_back(index2);
						        t->edges.push_back(index3);
						        t->A = edgeProbabilities[index1] * edgeProbabilities[index2] * edgeProbabilities[index3];
                                count++;		        		
		        		        pattern.push_back(t);		  
		    				}		    				
		    				iit=adjacencyList[itt].find_next(iit);
		    			}
		    		}
		    	    itt=adjacencyList[i].find_next(itt);
		        }
		        it=adjacencyList[i].find_next(it);
		    }
	    }
	}
    else if(caseN==4){
        for(int i=0;i<NodeSize;i++){
		    it=adjacencyList[i].find_first();
		    while(it!=dynamic_bitset<>::npos){
		        itt=adjacencyList[i].find_next(it);
		        while(itt!=dynamic_bitset<>::npos){
		    		iit=adjacencyList[i].find_next(itt);
		    		while(iit!=dynamic_bitset<>::npos){
		    			Motif* t=new Motif(count);
		    			index1=adjacencyMatrix[i][it];
		    			index2=adjacencyMatrix[i][itt];
		    			index3=adjacencyMatrix[i][iit];
		                t->edges.push_back(index1);
		                t->edges.push_back(index2);
						t->edges.push_back(index3);
						t->A = edgeProbabilities[index1] * edgeProbabilities[index2] * edgeProbabilities[index3];
                        count++;		        		
		        		pattern.push_back(t);
		        		iit=adjacencyList[i].find_next(iit);
		    		}
		    	    itt=adjacencyList[i].find_next(itt);
		        }
		        it=adjacencyList[i].find_next(it);
		    }		
	    }
    }
}

