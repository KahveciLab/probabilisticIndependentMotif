#include <iostream>
#include "Graph.h"
#include <set>
#include <boost/heap/fibonacci_heap.hpp>
#include <unordered_map>
#include <iomanip>
#include <new>
#include <armadillo>
#include <random>
#include <cmath>

#define NUM_SIZE 2097152

using namespace arma;

Graph* g;
Graph* g0;

//max heap, sort the priority of the motifs, from max to min
struct CompareNode : public std::binary_function<Motif*, Motif*, bool>
{
	bool operator()(const Motif* lhs, const Motif* rhs) const{
	    return lhs->priority < rhs->priority;
	}
};

typedef heap::fibonacci_heap<Motif*,heap::compare<CompareNode> >::handle_type cn;
typedef Mat<unsigned char> bmat;

void getUpdateEdge(vector<Motif*>& pattern, int motifIndex, dynamic_bitset<>& neighbors, dynamic_bitset<>& updateEdge){
	 updateEdge.reset();
     dynamic_bitset<>::size_type it;
     it = neighbors.find_first();
     Motif* a;
     while(it!=dynamic_bitset<>::npos){
     	a = pattern[it];
     	for(int i=0;i<a->edges.size();i++)
     		updateEdge[a->edges[i]]=1;
     	it=neighbors.find_next(it);
     }

     a = pattern[motifIndex];
     for(int i=0;i<a->edges.size();i++)
        updateEdge[a->edges[i]]=0;   
}

double calculateLoss(bmat& temp, bmat& edgeM, mat& edgeP, int edgeNum){
    bmat tt = temp; //this is used to determine which motifs will be completed
    double res = 0;
    bmat* result = new bmat(temp.n_elem,1);
    result->zeros();
    mat* rcoeff = new mat(1,1,arma::fill::ones);
    imat* y = new imat(1,1);
    y->zeros();
    bmat* nresult;
    mat* nrcoeff;
    imat* ny;
    uvec degree;

    int i = edgeNum-1;
    for(;i>=0;i--){    
        nresult = new bmat(*result);                    
        nrcoeff = new mat(*rcoeff);       
        ny = new imat(*y);

        *rcoeff *= edgeP(i);
        *nrcoeff *= (1-edgeP(i));
        
        degree = arma::find(edgeM.col(i)==1);
        tt.cols(degree) -= 1;

        int numComple = 0;
        for(int index=0; index < degree.n_elem; index++ ){
            int ii = degree(index);
            result->row(ii) += 1;

            if(tt(ii)==0)//complete
            {
                uvec q1 = arma::find(result->row(ii)==temp(ii));
                y->elem(q1) += 1;
                numComple++;
            }
        }
        
        result->insert_cols(result->n_cols, *nresult);
        rcoeff->insert_cols(rcoeff->n_cols, *nrcoeff);
        y->insert_rows(y->n_rows, *ny);

        delete nresult;
        delete nrcoeff;
        delete ny; 
          
    }

    mat mulres(*rcoeff * *y);
    res = mulres(0);

    delete result;
    delete rcoeff;
    delete y;
    tt.clear();

    return res;
}

void calculatePriority(vector<Motif*>& pattern, int motifIndex, dynamic_bitset<>& neighbors,vector<dynamic_bitset<> >& madjacencyList){
    Motif* t=pattern[motifIndex];
    int numNeighbors=neighbors.count();
    bmat temp(1,numNeighbors); //template
    temp.zeros();

    int j=0;//the index of edges
	int k=0;//the index of motifs
	map<int,int> edgeIndex;
	vector<int> diffEdge; //store the different edges' index
    bmat edgeM; //edge matrix
    bmat erow(numNeighbors,1); //size numneighbors X  numedges
    erow.zeros();
    mat edgeP; //probability of edges, size 1 X numedges
    mat prow(1,1);

	dynamic_bitset<>::size_type it;
	it=neighbors.find_first();
	while(it!=dynamic_bitset<>::npos){
	    Motif* a=pattern[it];
		
		for(int i=0;i<a->edges.size();i++){
			int eIndex = a->edges[i];
			if(madjacencyList[eIndex][motifIndex]==0)
				diffEdge.push_back(eIndex);
		}
	   
	    temp(k) = diffEdge.size(); 
		while(!diffEdge.empty()){ 
			int ii = diffEdge.back();
			diffEdge.pop_back();
		
			if(edgeIndex.find(ii)==edgeIndex.end()){
	        	edgeIndex.insert(pair<int, int>(ii,j));
                prow(0) = g->edgeProbabilities[ii];
                edgeP.insert_cols(j,prow);
                edgeM.insert_cols(j,erow);
                edgeM(k,j) = 1;
	            j++;
	        }
	        else{
	        	edgeM(k,edgeIndex[ii]) = 1;
	        }
		}
	    k++;	    
	    it=neighbors.find_next(it);
	}
	diffEdge.clear();
	erow.clear();
	edgeIndex.clear();  
	prow.clear(); 
	int edgeNum = j; //the number of edges
	
    double res= calculateLoss(temp,edgeM,edgeP,edgeNum);    
    t->priority = t->A/res;   
    temp.clear();    
    edgeP.clear();
    edgeM.clear();
}

/*
void calculatePriority(vector<Motif*>& pattern, int currIndex, dynamic_bitset<>& neighbors,vector<dynamic_bitset<> >& madjacencyList){
    Motif* t=pattern[currIndex];
    double mean=0;

    int numNeighbors = neighbors.count();
    unordered_map<int,int> edgeIndexMap; // originalEdgeIndex, currentIndex
    unordered_map<int,int> motifIndexMap; // currentIndex, originalMotifId
    vector<dynamic_bitset<> > edgeMatrix;
    vector<double> tempEdgeProbability;


    dynamic_bitset<>::size_type it;
    it=neighbors.find_first();
    int edgeIndex=0;
    int motifIndex=0;
    while(it!=dynamic_bitset<>::npos){
        Motif* a=pattern[it];
        
        double existProb=1;
        a->tempEdges.clear();
        for(int i=0;i<a->edges.size();i++){
            int eIndex = a->edges[i];
            if(madjacencyList[eIndex][currIndex]==0)
            {
                existProb *= g->edgeProbabilities[eIndex];
                
                if(edgeIndexMap.find(eIndex)==edgeIndexMap.end()){
                    dynamic_bitset<> erow(numNeighbors);
                    erow[motifIndex]=1;
                    edgeIndexMap.insert(pair<int,int>(eIndex,edgeIndex));
                    tempEdgeProbability.push_back(g->edgeProbabilities[eIndex]);
                    edgeIndex++;
                    edgeMatrix.push_back(erow);
                }
                else edgeMatrix[edgeIndexMap[eIndex]][motifIndex]=1;

                a->tempEdges.push_back(edgeIndexMap[eIndex]);               
            }                    
        }
        a->tempA=existProb;
        mean+=existProb;
        
        motifIndexMap.insert(pair<int,int>(motifIndex,a->id));
        motifIndex++;
        it=neighbors.find_next(it);
    }
    t->mean=mean;
    t->priority=t->A/t->mean;

    //calculate the variance
    double variance=0;
    it=neighbors.find_first();
    dynamic_bitset<>::size_type nit;
    motifIndex=0;
    while(it!=dynamic_bitset<>::npos){
        Motif* a=pattern[it];
        
        dynamic_bitset<> tdb(numNeighbors); //find motifs that have edges in common with motif a
        for(int i=0;i<a->tempEdges.size();i++){
            int eIndex=a->tempEdges[i];
            tdb |= edgeMatrix[eIndex];
        }

        double sendSum=0; //second summation of variance equation of Andrei's paper //counting motifs in probabilistic biological networks
        nit=tdb.find_first();
        while(nit!=dynamic_bitset<>::npos){
            Motif* b=pattern[motifIndexMap[nit]];
            double same=1;
            double diff=1;
            for(int i=0;i<b->tempEdges.size();i++){
                int eIndex=b->tempEdges[i];
                if(edgeMatrix[eIndex][motifIndex]==1)
                    same*=tempEdgeProbability[eIndex];
                else diff*=tempEdgeProbability[eIndex];
            }
            same = 1- same;
            sendSum+=same*diff;
            nit=tdb.find_next(nit);
        }

        variance+=a->tempA*sendSum;
        motifIndex++;
        it=neighbors.find_next(it);
    }
    t->std=sqrt(variance);

    cout<<currIndex<<"\t"<<t->priority<<"\t"<<t->mean<<"\t"<<t->std<<endl;
}
*/

double calculateF2(vector<Motif*>& pattern, int& countF2){
	int size=pattern.size();
	if(size==0)
		return 0;	
	double F2=0;
	
	//build overlap graph
	vector<dynamic_bitset<> > madjacencyList(g->EdgeSize, dynamic_bitset<>(size));
	dynamic_bitset<> neighbors(size);
	heap::fibonacci_heap<Motif*,heap::compare<CompareNode> > pq;//priority queue to store the priority of motifs
	map<int, cn > tool;

	for(int i=0;i<size;i++){
		Motif* a=pattern[i];
		for(int j=0;j<a->edges.size();j++)		
            madjacencyList[a->edges[j]][i] =1;
    }
    	
	for(int i=0;i<size;i++){
		Motif* a=pattern[i];
		neighbors.reset();
        for(int j=0;j<a->edges.size();j++){
			int index = a->edges[j];
			neighbors |= madjacencyList[index];
		}
		neighbors[i]=0;		
            			
		if(neighbors.none()){
			countF2++;
			F2 += a->A;
			delete a;
		}
		else{
            calculatePriority(pattern,i,neighbors,madjacencyList);
            tool.insert(pair<int, cn >(i,pq.push(a)));
        }
            					      
	}

    dynamic_bitset<> updateEdge(g->EdgeSize);
	dynamic_bitset<> nn(size);//record which embeddings are to update
	dynamic_bitset<>::size_type it;
	dynamic_bitset<> tdb(size);
		
	while(!pq.empty()){		
		Motif* a =pq.top();
        F2 += a->A;
	    countF2++;

		//delete a and its neighbor
        neighbors.reset();
        for(int j=0;j<a->edges.size();j++){
        	int index = a->edges[j];
            neighbors |= madjacencyList[index];
        }
        getUpdateEdge(pattern,a->id,neighbors,updateEdge);

	    it=neighbors.find_first();
	    while(it!=dynamic_bitset<>::npos){
            pq.erase(tool[it]); //remove motif from the priority queue	    				            
	        delete pattern[it];
	        it=neighbors.find_next(it);
	    }
          
        nn.reset();
        it=updateEdge.find_first();
        while(it!=dynamic_bitset<>::npos){
        	tdb = madjacencyList[it];
        	tdb &= neighbors;
        	madjacencyList[it] ^= tdb;
            nn |= madjacencyList[it];
            it = updateEdge.find_next(it);
        }

	    it=nn.find_first();
	    while(it!=dynamic_bitset<>::npos){
	        Motif* b=pattern[it];

	        neighbors.reset();
            for(int j=0;j<b->edges.size();j++){
			    int index = b->edges[j];
			    neighbors |= madjacencyList[index];
		    }		
            neighbors[it]=0;

		    if(neighbors.none()){
			    countF2++;
			    F2 += b->A;
			    pq.erase(tool[it]);
			    delete b;
		    }

		    else{			    
                calculatePriority(pattern,it,neighbors,madjacencyList);
                pq.update(tool[it]);
		    }			      
	        	
	        it=nn.find_next(it);
	    }
	}
    
    updateEdge.clear();
    nn.clear();
    tdb.clear();
    neighbors.clear();
	pattern.clear();
	madjacencyList.clear();
	tool.clear();

	return F2;
}

void DFS(dynamic_bitset<>& marked, dynamic_bitset<>& subgraph, int i){
    marked[i]=1;
    subgraph[i]=1;
    for(int j=0;j<g0->NodeSize;j++){
        if(g0->adjacencyList[i][j]==1 && marked[j] == 0)
            DFS(marked,subgraph,j);
    }
}

void DepthFirstSearch(){
    double pf1 = 0;
	double pf2 = 0;
    double pf3 = 0;
    double pf4 = 0;

	dynamic_bitset<> marked(g0->NodeSize);
	dynamic_bitset<> subgraph(g0->NodeSize);
	dynamic_bitset<>::size_type it;

    vector<Motif*> pattern;
    int count1 = 0;
    int count2 = 0;
    int count3 = 0;
    int count4 = 0;

    for(int i=0;i<g0->NodeSize;i++){
        if(marked[i]==0){
            subgraph.reset();
        	DFS(marked,subgraph,i);
        	int nodeNum = subgraph.count();

        	//then work on the subgraph with nodes in set res
        	if(nodeNum>=3){ 
        		if(nodeNum != g0->NodeSize)
        		    g = new Graph(subgraph, nodeNum, g0);
        		else g = g0;

                g->findMotifs(pattern, 1);            
                pf1 += calculateF2(pattern, count1);

                g->findMotifs(pattern, 2);
                pf2 += calculateF2(pattern, count2); 
                               
                g->findMotifs(pattern, 3);
                pf3 += calculateF2(pattern, count3);
             
                g->findMotifs(pattern, 4);
                pf4 += calculateF2(pattern, count4);
                
                if(nodeNum != g0->NodeSize)
                    delete g;
               
        	}
        	
        }
    }
    subgraph.clear();
    marked.clear();

/*
    cout<<"F2 count"<<endl;
    cout<<count1<<"\t"<<count2<<"\t"<<count3<<"\t"<<count4<<endl;
    cout<<"Probability sum"<<endl;*/
    cout<<pf1<<"\t"<<pf2<<"\t"<<pf3<<"\t"<<pf4<<"\t";
}

int main(int argc, char *argv[]){
	int numberNodes = atoi(argv[1]);
	char fileName[100];

    strcpy(fileName,argv[2]);
    g0 = new Graph(numberNodes, fileName);
    const clock_t begin_time=clock();
    DepthFirstSearch(); 
    cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl;
    delete g0;
    /*
	char buffer[5];
	for(int i=1;i<=1;i++){
		strcpy(fileName,argv[2]);
		sprintf(buffer,"%d",i);
		strcat(fileName,buffer);
		g0 = new Graph(numberNodes, fileName);
	    const clock_t begin_time=clock();
        DepthFirstSearch();	
	    cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl;
	    delete g0;
	}
	*/

	return 0;
}