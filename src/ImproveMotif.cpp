#include <iostream>
#include "Graph.h"
#include "Edge.h"
#include <set>
#include <boost/heap/fibonacci_heap.hpp>
#include <unordered_map>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#define NUM_LIMIT 18
#define TIME_LIMIT 50

Graph* g;
Graph* g0;
int randTime; //when the graph is too dense, number of times we randomly select 
//max heap, sort the priority of the motifs, from max to min
struct CompareNode : public std::binary_function<Motif*, Motif*, bool>
{
	bool operator()(const Motif* lhs, const Motif* rhs) const{
	    return lhs->priority < rhs->priority;
	}
};

//min heap, for the undetermined embeddings
struct CompareNode2 : public std::binary_function<Motif*, Motif*, bool>
{
	bool operator()(const Motif* lhs, const Motif* rhs) const{
	    return lhs->priority > rhs->priority;
	}
};
//sort sortedPattern from max to min
struct CompareDegree : public std::binary_function<Motif*, Motif*, bool>
{
	bool operator()(const Motif* lhs, const Motif* rhs) const{
	  	if(lhs->pri == rhs->pri) 	 	
	  	 	return lhs->id < rhs->id; 	 	
	    else return lhs->pri > rhs->pri;
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

void calculatePriorityForSizeTwo(vector<Motif*>& pattern, int motifIndex, dynamic_bitset<>& neighbors){
	Motif* t=pattern[motifIndex];
	int numNeighbors=neighbors.count();
	
	int firstEdge = t->edges[0];
	int secondEdge = t->edges[1];
	dynamic_bitset<> progress(numNeighbors+1);
	dynamic_bitset<> cprog(numNeighbors+1);

	map<int,int> edges; //edge index and how many motifs connected to it
    dynamic_bitset<>::size_type it;
    it = neighbors.find_first();
	while(it!=dynamic_bitset<>::npos){
	    Motif* a=pattern[it];

	    int temp = a->edges[0];
        if(!(temp != firstEdge && temp!= secondEdge))
			temp = a->edges[1];
		
	    if(edges.find(temp)==edges.end())
            edges.insert(pair<int,int>(temp,1));	    
	    else edges[temp] = edges[temp] + 1; 
	        
	    it=neighbors.find_next(it);
	}
    
    map<int,int>::iterator edgeit;
    rowvec result = zeros(1,numNeighbors+1);
    result(0) = 1;
    rowvec nresult = zeros(1,numNeighbors+1);

    progress[0]=1;

    for(edgeit = edges.begin(); edgeit != edges.end(); ++edgeit){
        int index = edgeit->first;
        int tnumber = edgeit->second;

        cprog.reset();
        it = progress.find_first();
        while(it!=dynamic_bitset<>::npos){
        	double tempr = result(it);        	
        	result(it) = tempr * (1-g->edgeProbabilities[index]);
 
            int y = (int)it + tnumber; 
            nresult(y) = nresult(y) + tempr * g->edgeProbabilities[index];

            cprog[y] = 1;
            it = progress.find_next(it);
        }
        
        result = result + nresult;
        nresult.fill(0);
        progress |= cprog;
    }

    it = progress.find_first();
    double res=0;
    while(it!=dynamic_bitset<>::npos){
        res += (int)it * result(it);
        it = progress.find_next(it);
    }

    result.clear();
    nresult.clear();
    progress.clear();
    cprog.clear();
    edges.clear();

    t->priority = t->A/res;
}

int chooseEdge(mat& edgeM, bmat& tt){
    double pri = -1;
	double sePri = -1;
	int bIndex=-1;
    int k = edgeM.n_rows-2;
    uvec degree;
    uvec q1;

	for(int j=0;j<edgeM.n_cols;j++){
		int ndegree = 0;
		int mcomplete = 0;
		if(k>=1){
            degree= arma::find(edgeM.col(j).rows(0,k-1)==1); //find which motifs need this edge
            ndegree = degree.n_elem;
            q1 = (arma::find(tt.cols(degree)==1));
            mcomplete = q1.n_elem;
		}
	    double a= edgeM(k+1,j) * (mcomplete+edgeM(k,j)); 
	    double b= edgeM(k+1,j) * (ndegree+edgeM(k,j));    
	    if(a > pri){
	        pri = a;
	        sePri = b;
	        bIndex = j;
	    }
	    else if(a == pri){
	        if(b > sePri){
	            sePri = b;
	            bIndex = j;
	        }
	    }
    }
    return bIndex;
}
bool calculatePriority(vector<Motif*>& pattern, int motifIndex, dynamic_bitset<>& neighbors, double best, vector<dynamic_bitset<> >& madjacencyList){
    Motif* t=pattern[motifIndex];
    
    //build Edge matrix
    int numNeighbors=neighbors.count();
    bmat temp(1,numNeighbors);
    temp.zeros();

	int j=0;//the index of edges
	int k=0;//the index of motifs
	int i=0;
	map<int,int> edgeIndex;
	vector<int> diffEdge; //store the different edges' index

    mat edgeM; //edge matrix
    mat erow(numNeighbors+2,1); //size numneighbors+1(t value) X  numedges
    erow.zeros();

	dynamic_bitset<>::size_type it;
	it=neighbors.find_first();
	while(it!=dynamic_bitset<>::npos){
	    Motif* a=pattern[it];
		
		for(i=0;i<a->edges.size();i++){
			int eIndex = a->edges[i];
			if(madjacencyList[eIndex][motifIndex]==0)
				diffEdge.push_back(eIndex);
		}

	    if(diffEdge.size() == 1){
			i = diffEdge.back();
			diffEdge.pop_back();
	        if(edgeIndex.find(i)==edgeIndex.end()){
	           	edgeIndex.insert(pair<int, int>(i,j));
                edgeM.insert_cols(j,erow);
                edgeM(numNeighbors,j)=1;
                edgeM(numNeighbors+1,j)=g->edgeProbabilities[i];
	           	j++;
	        }	        
	        else edgeM(numNeighbors,edgeIndex[i]) = edgeM(numNeighbors,edgeIndex[i]) + 1;
	    }
	    
	    else{
	    	temp(k) = diffEdge.size();
			while(!diffEdge.empty()){ 
				i = diffEdge.back();
			    diffEdge.pop_back();
				
				if(edgeIndex.find(i)==edgeIndex.end()){
	        		edgeIndex.insert(pair<int, int>(i,j));	        		
                    edgeM.insert_cols(j,erow);
                    edgeM(k,j)=1;
                    edgeM(numNeighbors+1,j)=g->edgeProbabilities[i];
	        		j++;
	        	}
	        	else{
	        		edgeM(k,edgeIndex[i]) = 1;
	        	}
			}
	        k++;
	    }

	    it=neighbors.find_next(it);
	}
	diffEdge.clear();
	erow.clear();
	edgeIndex.clear();

	if(k<numNeighbors){
		edgeM.shed_rows(k,numNeighbors-1);
		temp.shed_cols(k,numNeighbors-1);
		numNeighbors = k;
	}

	bmat tt = temp; //this is used to determine which motifs will be completed
	double res = 0;
	bmat* result = new bmat(temp.n_elem,1);
	result->zeros();
    mat* rcoeff = new mat(1,1,fill::ones);
    imat* y=new imat(1,1);
    y->zeros();
    bmat* nresult;
	mat* nrcoeff;
    imat* ny;
    uvec degree;
    double res1, res2, calRes;
    calRes = 0;
	double remainRes =0; //the loss value sum for the remaining matrices
    vector<bmat*> remain; //record which matrices have not been calculated
	vector<mat*> remainCoeff; //corresponding coefficient value 
	vector<imat*> remainy;//corresponding y value
	vector<double> remainValue; //the current loss value for the remaining matrices 
	vector<int> stackProgress;//record which part that has not been calculated
	vector<Edge> midSate; //record which edges that are needed to be recorded 
	mat stackBit;//record which edges have not been dealt with
  	 
  	//for the previous NUM_LIMIT edges, no collapse and no recursive strategy    
	for(i=0;i<j;i++){
        nresult = new bmat(*result);
        nrcoeff = new mat(*rcoeff);
        ny = new imat(*y);

	    int eIndex = chooseEdge(edgeM,tt);
	    *rcoeff *= edgeM(numNeighbors+1,eIndex);
	    *nrcoeff *= (1-edgeM(numNeighbors+1,eIndex));
	    mat mulres(*nrcoeff * *ny);
	    res1 = mulres(0);

	    degree = arma::find(edgeM.col(eIndex).head(numNeighbors)==1);
	    tt.cols(degree) -= 1;
	    *y += (int)(edgeM(numNeighbors,eIndex));
	    	
	    for(int index =0; index < degree.n_elem; index++){
	    	int ii = degree(index);
	    	result->row(ii) += 1;

	    	if(tt(ii)==0){
	    		uvec q1 = arma::find(result->row(ii)==temp(ii));
	    		y->elem(q1) += 1;
	    	}
	    }
        mulres = (*rcoeff * *y);
        res2 = mulres(0);
        res = res1+res2+remainRes;

        t->priority = t->A/res;
        if(t->priority <= best){
            delete nresult;
            delete nrcoeff;
            delete ny;
            break;
        }
        else if(i<=NUM_LIMIT){
        	if(i==NUM_LIMIT){
                degree = arma::find(tt != 0);
                if(degree.n_elem < tt.n_elem){   
                    tt = tt.cols(degree);
                    temp = temp.cols(degree);
                    *result = result->rows(degree);
                    *nresult = nresult->rows(degree);
                    degree.resize(degree.n_elem+2);
                    degree(degree.n_elem-2) = numNeighbors;
                    degree(degree.n_elem-1) = numNeighbors+1;
                    edgeM = edgeM.rows(degree);
                    numNeighbors = degree.n_elem-2;
                }
            }
            result->insert_cols(result->n_cols, *nresult);
            rcoeff->insert_cols(rcoeff->n_cols, *nrcoeff);
            y->insert_rows(y->n_rows, *ny);

            delete nresult;
            delete nrcoeff;
            delete ny;        	
        }
        else{
            Edge e;
            e.degree = edgeM.col(eIndex);
            e.temp=zeros(1,numNeighbors);
            uvec q1 = arma::find(tt.cols(degree)==0);
            q1 = degree.elem(q1);
            e.temp.cols(q1) = (conv_to<mat>::from(temp)).cols(q1);
            midSate.push_back(e);
            if(i==j-1){//last edge, no need to store
               calRes = res1+res2;
               delete nresult;
               delete nrcoeff;
               delete ny;
            }
            else{
            	stackBit.insert_cols(stackBit.n_cols, edgeM.col(eIndex).head(numNeighbors));
            	remain.push_back(nresult);
            	remainCoeff.push_back(nrcoeff);
                remainy.push_back(ny);
                remainValue.push_back(res1);
                remainRes += res1;
                stackProgress.push_back(i-NUM_LIMIT);
            }         
        }

        edgeM.shed_col(eIndex);       
	}

    delete result;
    delete rcoeff;
    delete y;
    temp.clear();
    tt.clear();
    edgeM.clear();
    //stop because priority value is smaller than the current best
	if(i<j){        
        midSate.clear();
        stackBit.clear();
        remainValue.clear();
        bool flag = false;
        if(stackProgress.empty())
        	flag = true;
        stackProgress.clear();

        for(int m=0;m<remain.size();m++){
        	delete remain[m];
        	delete remainCoeff[m];
        	delete remainy[m];
        }
        remain.clear();
        remainCoeff.clear();
        remainy.clear();


        if(i==j-1 && flag)
        	return true;
        else return false;
	}
	else{
		if(stackProgress.empty())
			return true;
		else{
			int time = 0;
			int edgeNum = midSate.size();
			while(!stackProgress.empty()){
				i = stackProgress.back();
                stackProgress.pop_back();

                result = remain.back();
                remain.pop_back();
                rcoeff = remainCoeff.back();
                remainCoeff.pop_back();
                y = remainy.back();
                remainy.pop_back();
                colvec bitcol = stackBit.col(stackBit.n_cols-1);
                stackBit.shed_col(stackBit.n_cols-1);
                remainRes -= remainValue.back();
                remainValue.pop_back();

                for(;i<edgeNum;i++){
                	Edge e = midSate[i];                	
                	degree = arma::find(e.degree.head(numNeighbors)==1);
                	uvec q1 = arma::find(bitcol.elem(degree)==0);
                	if(q1.n_elem==0 && e.degree(numNeighbors)==0){
                		//no need to expand                	    
                	    mat mulres(*rcoeff * *y);
	                    res = mulres(0)+remainRes+calRes;
                	}
                	else{
                		nresult = new bmat(*result);
                        nrcoeff = new mat(*rcoeff);
                        ny = new imat(*y);

                        *rcoeff *= e.degree(numNeighbors+1);
	                    *nrcoeff *= (1-e.degree(numNeighbors+1));
	                    mat mulres(*nrcoeff * *ny);
	                    res1 = mulres(0);
	                    colvec tempBit(bitcol);

                        *y += (int)(e.degree(numNeighbors));
                		degree = degree.elem(q1);
                		for(int index =0;index < degree.n_elem;index++){
                            int ii = degree(index);
                            result->row(ii) += 1;
                            tempBit(ii) = 1;

                            if(e.temp(ii)>0){
                            	q1 = arma::find(result->row(ii)==e.temp(ii));
                            	y->elem(q1) +=1;
                            }                            
                		}

                		mulres = (*rcoeff * *y);
                		res2 = mulres(0);
                        res = res1+res2+remainRes+calRes;

                        if(i!=edgeNum-1){
                        	stackBit.insert_cols(stackBit.n_cols, tempBit);
                        	remain.push_back(nresult);
            	            remainCoeff.push_back(nrcoeff);
                            remainy.push_back(ny);
                            remainValue.push_back(res1);
                            remainRes += res1;
                            stackProgress.push_back(i+1);
                        }
                        else{
                        	calRes += res1 + res2;
                        	delete nresult;
                        	delete nrcoeff;
                        	delete ny;
                        }
                        time ++;

                	}

	                t->priority = t->A/res;
                    if(t->priority <= best)           
                        break;  

                    if(time == TIME_LIMIT)
                        break;      
                }
                delete result;
                delete rcoeff;
                delete y;
                //stop if less than the current best or the count is too much 
                if(i<edgeNum){
                    midSate.clear();
                    stackBit.clear();
                    remainValue.clear();
                    for(int m=0;m<remain.size();m++){
        	            delete remain[m];
        	            delete remainCoeff[m];
        	            delete remainy[m];
                    }
                    remain.clear();
                    remainCoeff.clear();
                    remainy.clear();

                    if(t->priority <= best){
                        if(i==edgeNum-1 && stackProgress.empty())
                        	return true;
                        else return false;
                    }
                    else return false;
                }
			}
			midSate.clear();
            return true;
		}
	}
	           
}

//deal with the specified motif pattern with three nodes and two edges
double calculateF2ForSizeTwo(vector<Motif*>& pattern, int& countF2){
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
            calculatePriorityForSizeTwo(pattern,i,neighbors);
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
                calculatePriorityForSizeTwo(pattern,it,neighbors);
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

//deal with motif patterns with four nodes and three edges 
double calculateF2(vector<Motif*>& pattern, int& countF2){
	int size=pattern.size();
	if(size==0)
		return 0;
	if(size==1){
		countF2++;
        double res = pattern[0]->A;
        delete pattern[0];
        pattern.clear();
	    return res;	
	}

	double F2=0;

	//build overlap graph
	vector<dynamic_bitset<> > madjacencyList(g->EdgeSize, dynamic_bitset<>(size));
	for(int i=0;i<size;i++){
		Motif* a=pattern[i];
		for(int j=0;j<a->edges.size();j++)		
            madjacencyList[a->edges[j]][i] =1;
    }


    dynamic_bitset<> neighbors(size);
	heap::fibonacci_heap<Motif*,heap::compare<CompareNode> > pq;//priority queue to store the priority of motifs
	heap::fibonacci_heap<Motif*,heap::compare<CompareNode2> > pq2;//priority queue for undetermined embeddings
	set<Motif*, CompareDegree> sortPQ; 
	map<int, cn > tool;
	double best=-1;
	
	
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
            a->pri = a->A / neighbors.count();
            sortPQ.insert(a);
		}			
	}

    if(sortPQ.empty()){
    	neighbors.clear();
        madjacencyList.clear();
    	pattern.clear();
    	return F2;
    }

    //calculate the priority value of top k embeddings
    int k = ( 100 > sortPQ.size()/10 ) ? 100 : sortPQ.size()/10;
    int count = 0;
	set<Motif*,CompareDegree>::iterator mit;
	for(mit = sortPQ.begin(); mit != sortPQ.end() && count < k; ++mit, ++count){
        Motif* a = *mit;
        int id= a->id;

    	neighbors.reset();
    	for(int j=0;j<a->edges.size();j++){
			int index = a->edges[j];
			neighbors |= madjacencyList[index];
		}
        neighbors[id]=0;

        if(calculatePriority(pattern,id,neighbors,best,madjacencyList)){
        	if(a->priority > best){
                best = a->priority;
                while(!pq2.empty()&&pq2.top()->priority <= best){
                	pq2.pop();
                }
        	}        	
        	tool.insert(pair<int, cn>(id,pq.push(a)));  
        }
        else{
            if(a->priority > best){
            	pq2.push(a);
            }
        }       
	}

    dynamic_bitset<> updateEdge(g->EdgeSize);
	dynamic_bitset<> nn(size);//record which embeddings are to update
	dynamic_bitset<>::size_type it;
	map<int, cn>::iterator pqit;
	dynamic_bitset<> tdb(size);
	srand(time(NULL));
	int randNumber;

    while(!sortPQ.empty() || !pq.empty()){
    	Motif* a;
        if(!pq2.empty()){
            double maxA = -1;
            if(!pq.empty()){                
                a = pq.top();
                maxA = a->A;
            }
           	              
            while(!pq2.empty()){
            	Motif* teM = pq2.top();
            	pq2.pop();
            	if(maxA < teM->A){
            		maxA = teM->A;
            		a = teM;
            	}
            }
            randTime++;
        }                
        else
           a = pq.top();           
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
	    	sortPQ.erase(pattern[it]);
            if( tool.find(it) != tool.end() ){
                pq.erase(tool[it]);
                tool.erase(it);
            }            	            				            
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
        
        //first update priority queue
        for(pqit = tool.begin(); pqit != tool.end(); ){
            Motif* b = pattern[pqit->first];
            int id = b->id;

            if(nn[id] == 0)//no update
            {
            	++pqit;
            	continue;
            }
            nn[id] = 0;

            neighbors.reset();
            for(int j=0;j<b->edges.size();j++){
			    int index = b->edges[j];
			    neighbors |= madjacencyList[index];
		    }
		    neighbors[id]=0;

            if(neighbors.none()){
		        countF2++;
            	F2 += b->A;
                pq.erase(pqit->second);
                pqit = tool.erase(pqit);
                sortPQ.erase(pattern[id]);
            	delete b;
            }
            else{
            	//update sortPQ           	
            	sortPQ.erase(pattern[id]);
            	b->pri = b->A/neighbors.count();
            	sortPQ.insert(pattern[id]);

            	calculatePriority(pattern,id,neighbors,-1,madjacencyList);
            	pq.update(pqit->second);
            	++pqit;
            }       	
        }

        if(pq.empty())
        	best = -1;
        else best = pq.top()->priority;

        //then update nn's pri
        it = nn.find_first();
        while(it!= dynamic_bitset<>::npos){
        	Motif* b= pattern[it];
        	int id = b->id;

        	neighbors.reset();
            for(int j=0;j<b->edges.size();j++){
			    int index = b->edges[j];
			    neighbors |= madjacencyList[index];
		    }
		    neighbors[id]=0;

		    if(neighbors.none()){
		        countF2++;
            	F2 += b->A;
                sortPQ.erase(pattern[id]);
            	delete b;
            }
            else{
            	//update sortPQ            	
            	sortPQ.erase(pattern[id]);
            	b->pri = b->A/neighbors.count();
            	sortPQ.insert(pattern[id]);
            } 

            it = nn.find_next(it);
        }

        //redertimine the k and refill the sortPQ
        k = ( 100 > sortPQ.size()/10 ) ? 100 : sortPQ.size()/10;
	    for(mit = sortPQ.begin(),count=0; mit != sortPQ.end() && count < k; ++mit, ++count){
            Motif* b = *mit;
            int id= b->id;
            
            if(tool.find(id) == tool.end()){
            	if(nn[id]==0 && b->priority <= best)
            		continue;

            	neighbors.reset();
    	        for(int j=0;j<b->edges.size();j++){
			        int index = b->edges[j];
			        neighbors |= madjacencyList[index];
		        }
                neighbors[id]=0;

                if(calculatePriority(pattern,id,neighbors,best,madjacencyList)){
                	if(b->priority > best){
                		best = b->priority;
                		while(!pq2.empty()&&pq2.top()->priority <= best){
                	        pq2.pop();
                        }
                	}               	
                	tool.insert(pair<int,cn>(id,pq.push(b)));
                } 
                else{
                	if(b->priority > best)
                		pq2.push(b);
                }              
            }
	    }

    }
    neighbors.clear();
    updateEdge.clear();
    nn.clear();
    tdb.clear();
    tool.clear();
    pq.clear();
    sortPQ.clear();
	pattern.clear();
	madjacencyList.clear();

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
    double pf1 = 0; //to record the expected number of independent embeddings of motif 1
	double pf2 = 0;
    double pf3 = 0;
    double pf4 = 0;

	dynamic_bitset<> marked(g0->NodeSize);
	dynamic_bitset<> subgraph(g0->NodeSize);

    vector<Motif*> pattern;
    int count1 = 0; //to record how many independent embeddings (this is an integer)
    int count2 = 0;
    int count3 = 0;
    int count4 = 0;

    for(int i=0;i<g0->NodeSize;i++){
        if(marked[i]==0){
            subgraph.reset();
        	DFS(marked,subgraph,i);
        	int nodeNum = subgraph.count();

            //for each disconnected component, calculate its f2 count 
        	//then work on the subgraph with nodes in set res
        	if(nodeNum>=3){ 
        		if(nodeNum != g0->NodeSize)
        		    g = new Graph(subgraph, nodeNum, g0);
        		else g = g0;
        	                  
//        		cout<<"Subgraph with node size: "<<subgraph.count()<<endl;
//        		cout<<"Starting calculating pattern 1"<<endl;
                g->findMotifs(pattern, 1);
                pf1 += calculateF2ForSizeTwo(pattern, count1);

//              cout<<"Starting calculating pattern 2"<<endl;
                g->findMotifs(pattern, 2);
                pf2 += calculateF2(pattern, count2); 
                       
//              cout<<"Starting calculating pattern 3"<<endl;               
                g->findMotifs(pattern, 3);
                pf3 += calculateF2(pattern, count3);
                                      
//              cout<<"Starting calculating pattern 4"<<endl;
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
    cout<<"Probability sum"<<endl;
    */
    cout<<pf1<<"\t"<<pf2<<"\t"<<pf3<<"\t"<<pf4<<"\t";
//    cout<<"random selection times: "<<randTime<<endl;
}

int main(int argc, char *argv[]){
	int numberNodes = atoi(argv[1]);
	char fileName[100];
    strcpy(fileName,argv[2]);
    g0 = new Graph(numberNodes, fileName);
    const clock_t begin_time=clock();
    DepthFirstSearch(); 
    cout <<float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl;
    delete g0;
	return 0;
}