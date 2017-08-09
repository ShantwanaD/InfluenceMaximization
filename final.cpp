#include <iostream>
#include <fstream>
#include <cstdlib>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>
#include <time.h>
#include <sys/stat.h>
#include <string>
#include "table.cpp"

using namespace std;


//structure for the priority queue in Influence Maximization
struct tup
{
    int u;
    float mg;
    int flag;
    bool operator<(const tup& rhs) const
    {
        return mg < rhs.mg;
    }
};

//structure for the priority queue in Page Rank
struct tup_pr
{
    int u;
    double r;
    bool operator<(const tup_pr& rhs) const
    {
        return r < rhs.r;
    }
};

//Returns the union of two vectors
void set_union(vector<int> *s1, vector<int> *s2)
{
     vector<int>::iterator it;
     for(it = (*s2).begin(); it != (*s2).end(); ++it){
      if(find((*s1).begin(), (*s1).end(), *it) == (*s1).end())
        (*s1).push_back(*it);
     }
}

//Returns the set difference of two vectors
void set_difference(vector<int> *s1, vector<int> *s2)
{
     vector<int>::iterator it;
     for(it = (*s2).begin(); it != (*s2).end(); ++it){
      if(find((*s1).begin(), (*s1).end(), *it) != (*s1).end())
        (*s1).erase((find((*s1).begin(), (*s1).end(), *it)));
     }
}

//Generates 'sim' graph instances with edge weight 'prob' 
void generateGraphs(float prob, int sim){
	ifstream myfile; 
	ofstream outfile;
	int i, a, b; float r;
	string OutputFolder = to_string(prob);
	mkdir(OutputFolder.c_str(), S_IRWXU);
	string filename = to_string(prob) + "/g";
	for(i = 0; i < sim; i++){
	 filename += to_string(i); 
	 outfile.open(filename, ios::out);
	 myfile.open("data.txt", ios::in);
	 while((myfile >> a) && (myfile >> b)){
	  r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
      if(r <= prob){
        outfile << a << " " << b << "\n";
      }
	 }
	 filename = to_string(prob) + "/g";
	 outfile.close(); myfile.close();
	}
}


//Makes the adjacency list of all the graph instances & holds it in an array(of size sim) of maps. 
void graphs(int sim, float prob, map<int, vector<int> > *g){
  int i, a, b;
  string filename = to_string(prob) + "/g";;
  for (i = 0; i < sim; i++){
      ifstream myfile;
      filename += to_string(i); 
    myfile.open(filename, ios::in);
    while((myfile >> a) && (myfile >> b)){
      g[i][a].push_back(b);
    }
    myfile.close();
    filename = to_string(prob) + "/g";
     }
   }

//Finds the influence of given node in all the MC simulations and holds them in a vector of size 'sim'
void  findInfluence(map<int, vector<int> > *graphs, vector<vector<vector<int> > > *infset, int node, int sim, vector<vector<int> > *max_infset){
  int i, n;
  queue<int> q; 
  vector<int> s;
  vector<int>::const_iterator iter;
  map<int,bool>::const_iterator it;
  for(i = 0; i < sim; i++){
    map<int,bool> visited; 
    q.push(node);
    visited[node] = true;
    while(!q.empty()){
      n = q.front(); q.pop();
      if(graphs[i].find(n) != graphs[i].end()){
        s = graphs[i].find(n) -> second;
        for(iter = s.begin(); iter != s.end(); iter++){
          if(!visited[*iter]){
            q.push(*iter);
            visited[*iter] = true;
          }
        }
      }
    }
    for(it = visited.begin(); it != visited.end(); ++it){
    if(it -> second)
      (*infset)[i][node].push_back(it -> first);
  }
  set_difference(&((*infset)[i][node]), &((*max_infset)[i]));
 }
}

//Greedy Algorithm for influence maximization
vector<int> greedy(map<int, vector<int> >* graphs, int k, int sim){
  vector<int> A, tmp; 
  vector<vector<int> > max_infset(sim);
  float max_influence, inf;
  int i, j, x, s, node;
  for(i = 0; i < k; i++){
  	vector<vector<int> > tmpset(sim);
  	vector<vector <vector <int> > > infset(sim);
	for(x = 0; x < sim; x++){
		infset[x].resize(15229);
	}
    max_influence = -1;
    for(j = 0; j < 15229; j++){
      findInfluence(graphs, &infset, j, sim, &max_infset);
      inf = 0;
      for(x = 0; x < sim; x++){
      	inf += infset[x][j].size();
      }
      inf /= sim;
      if(inf > max_influence){
        max_influence = inf;
        for(x = 0; x < sim; x++){
      		tmpset[x] = infset[x][j];
      	}
        node = j;
      }
    }
    for(x = 0; x < sim; x++){
    	set_union(&max_infset[x], &tmpset[x]);
	}
    A.push_back(node);
    cout << "Node: " << node << "\tGain: " << max_influence << "\n";
    for(j = 0; j < sim; j++){
      graphs[j].erase(node);
    }
  } 
  return A;

}

//Optimised greedy using CELF algorithm
vector<int> celf(map<int, vector<int> >* graphs, int k, int sim){
  vector<int> A, tmp; 
  vector<vector<int> > max_infset(sim);
  float  inf, sum = 0;
  int i, j, x, s, node;
   vector<vector <vector <int> > > infset(sim);
  for(x = 0; x < sim; x++){
    infset[x].resize(15229);
  }
  priority_queue<tup> pq; 
  tup v;
  for(j = 0; j < 15229; j++){
    findInfluence(graphs, &infset, j, sim, &max_infset);
      inf = 0;
      for(x = 0; x < sim; x++){
        inf += infset[x][j].size();
      }
      inf /= sim;
    tup t = {j, inf, 0};
    pq.push(t);
  }
  cout << "k\tNode\tGain\tSpread\n";
  while(A.size() < k){
    v = pq.top(); pq.pop();
    if(v.flag == A.size()){
      A.push_back(v.u);
      sum += v.mg;
      cout << A.size() << "\t" << v.u << "\t" << v.mg << "\t" << sum << "\n";
      for(x = 0; x < sim; x++){
         set_union(&max_infset[x], &infset[x][v.u]);
      } 
    }
    else{
      for(x = 0; x < sim; x++){
         infset[x][v.u].clear();
      }
      findInfluence(graphs, &infset, v.u, sim, &max_infset);
      inf = 0;
      for(x = 0; x < sim; x++){
        inf += infset[x][v.u].size();
      }
      inf /= sim;
      v.mg = inf;
      v.flag = A.size();
      pq.push(v);
    }
  }
  return A;
}

// Returns the influence of k highest pagerank nodes 
vector<float> pagerank(map<int, vector<int> >* graphs, float prob, int sim, int k){
	Table t;
	int i, j; double sum;
	vector<int>::iterator it;
	vector<int> v; vector<float> spread;
	vector<vector<double> > p_ranks(sim);
	for(i = 0; i < sim; i++){
      p_ranks[i].resize(15229);
	}
    string filename = to_string(prob) + "/g";
    for (i = 0; i < sim; i++){
      filename += to_string(i); 
      t.reset();
      t.set_delim(" ");
	t.set_numeric(true);
      t.read_file(filename);
      t.pagerank();
	  for(j = 0; j < 15229; j++){
	  	p_ranks[i][j] = t.pr[j];
	  }
	  filename = to_string(prob) + "/g";
	 }
	 priority_queue<tup_pr> pq;
	 double ranks[15229]; tup_pr t1;
	 for(i = 0; i < 15229; i++){
	 	sum = 0;
	 	for(j = 0; j < sim; j++){
	 		sum += p_ranks[j][i];
	 	}
	 	ranks[i] = sum / sim;
	 	t1 = {i, ranks[i]};
	 	pq.push(t1);
	 }
	 for(i = 0; i < k; i++){
	 	t1 = pq.top(); pq.pop();
	 	v.push_back(t1.u);
	 }
	 cout << "\nk highest pagerank nodes: ";
	 for(it = v.begin(); it != v.end(); it++){
	 	cout << *it << " ";
	 }
	 vector<vector<int> > max_infset(sim);
	 vector<vector <vector <int> > > infset(sim);
	  for(i = 0; i < sim; i++){
	    infset[i].resize(15229);
	  }
	 for(it = v.begin(); it != v.end(); it++){
	 	findInfluence(graphs, &infset, *it, sim, &max_infset);
	 	sum = 0;
	      for(i = 0; i < sim; i++){
	        sum += infset[i][*it].size();
	      }
	      sum /= sim;
	      spread.push_back(sum);
	      for(i = 0; i < sim; i++){
	         set_union(&max_infset[i], &infset[i][*it]);
	      } 
	 }
	 return spread;
}

int main()
{
 clock_t start = clock(); int i;
 vector<int>::const_iterator iter;
 vector<float>::const_iterator iter1;
 int sim = 10, k = 200;
 float prob = 0.1, sum = 0;
 // generateGraphs(prob, sim);
 map<int, vector<int> > g[sim];
 graphs(sim, prob, g);
 cout << "Graph Instances Generated!\n";
 vector<int> A = celf(g, k, sim);
 
 cout << "\n\nTime elapsed: " << (double)(clock() - start)/CLOCKS_PER_SEC << " seconds.\n\n";

 vector<float> spread = pagerank(g, prob, sim, k);
 cout << "\nFinal Seed Set using Influence Maximisation: ";
 for(iter = A.begin(); iter != A.end(); ++iter){
  cout << *iter << " ";
 }
 cout << "\nSpread of k highest pagerank nodes:\n k\tSpread\n";
 i = 0;
 for(iter1 = spread.begin(); iter1 != spread.end(); ++iter1){
 	i++;
 	sum += *iter1;
  cout << i << "\t" << sum << "\n";
 }
 cout  << "\n";
 return 0;           
}
