#include <stdlib.h>
#include <cstdio>
#include <stdbool.h>
#include <string.h>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <windows.h>

#define NLINKS 100000000
using namespace std;

class edge
{
public:
	int s;
	int t;
};



class Graph{
public:
	void readedgelist(string edgelist);
	void coreDecomposition();
	void mkGraph();
	int outLargeClique();
	bool isEdge(int, int);
	Graph mksub(Graph &, vector<int>);
	//bool isEdge(int, int);

	int n;
	int e;
	int maxDeg;
	vector<edge> edges;

	int* deg;
	int* cd;
	int* adj;
	int *coreRank;			//increasing core number order
	int *coreNum;			//coreNum[i] is the core number of node i.
	int *bin;
};

inline int max3(int a, int b, int c) {
	a = (a>b) ? a : b;
	return (a>c) ? a : c;
}



void Graph::readedgelist(string edgelist) {
	
	int e1 = NLINKS;
	n = 0;
	e = 0;
	edges.resize(e1);
	ifstream file;
	file.open(edgelist);

	while (file >> edges[e].s >> edges[e].t)
	{
		n = max3(n, edges[e].s, edges[e].t);
		e++;
		if (e == e1) {
			e1 += NLINKS;
			edges.resize(e1);
		}
	}
	file.close();
	n++;
	edges.resize(e);
}

void Graph::mkGraph()
{
	deg = new int[n];
	cd = new int[n + 1];
	adj = new int[2 * e];
	maxDeg = 0;
	for (int i = 0; i < e; i++)
	{
		deg[edges[i].s]++;
		deg[edges[i].t]++;
		maxDeg = max3(maxDeg, deg[edges[i].s], deg[edges[i].t]);
	}

	cd[0] = 0;
	for (int i = 1; i < n + 1; i++) {
		cd[i] = cd[i - 1] + deg[i - 1];
		deg[i - 1] = 0;
	}
	for (int i = 0; i < e; i++) {
		adj[cd[edges[i].s] + deg[edges[i].s]++] = edges[i].t;
		adj[cd[edges[i].t] + deg[edges[i].t]++] = edges[i].s;
	}
	for (int i = 0; i < n; i++) sort(adj + cd[i], adj + cd[i] + deg[i]);
}

bool cmp(const pair<int, int> &a, const pair<int, int> &b)
{
	return a.second > b.second;
}

bool Graph::isEdge(int a, int b)
{
	if (deg[a] > deg[b]) a = a^b, b = a^b, a = a^b;
	for (int i = cd[a]; i < cd[a] + deg[a]; i++)
		if (adj[i] == b) return true;
	return false;
}
int Graph::outLargeClique()
{
	int CSize = 0;
	for (int i = n-1; i >= 0; i--)
	{
		int id = coreRank[i];
		if (coreNum[id] >= CSize)
		{
			pair<int, int> *SCore = new pair<int, int>[deg[id]];
			//int *S = new int[deg[id]], cnt = 0;
			int cnt = 0, ind = 0;
			for (int j = cd[id]; j < cd[id] + deg[id]; j++)
				if (coreNum[adj[j]] >= CSize)
				{
					SCore[cnt].first = adj[j];
					SCore[cnt].second = coreNum[adj[j]];
					cnt++;
				}
					//S[cnt++] = adj[j];
			sort(SCore, SCore + cnt, cmp);
			//sort(S, S + cnt, cmp);
			int *C = new int[deg[id]];
			for (int j = 0; j < cnt; j++)
			{
				int flag = 1;
				for (int k = 0; k < ind; k++)
				{
					if (isEdge(SCore[j].first, C[k]) == false)
					{
						flag = 0;
						break;
					}
				}
				if (flag) C[ind++] = SCore[j].first;
			}
			ind++;	//node "id" ?
			if (ind > CSize) CSize = ind;	
		}
	}
	return CSize;
}

void Graph::coreDecomposition()
{
	bin = new int[maxDeg + 2]();
	//int *bin = (int *)calloc(g->maxDeg + 2, sizeof(int));


	for (int i = 0; i < n; i++) bin[deg[i]]++;
	int lastBin = bin[0], nowBin;
	bin[0] = 0;
	for (int i = 1; i <= maxDeg; i++)
	{
		nowBin = lastBin + bin[i - 1];
		lastBin = bin[i];
		bin[i] = nowBin;
	}
	int *vert = new int[n], *pos = new int[n], *tmpDeg = new int[n];
	//int *vert = (int *)malloc(g->n * sizeof(int)), *pos = (int *)malloc(g->n * sizeof(int)), *tmpDeg = (int *)malloc(g->n * sizeof(int));
	for (int i = 0; i < n; i++)
	{
		pos[i] = bin[deg[i]];
		vert[bin[deg[i]]++] = i;
		tmpDeg[i] = deg[i];
	}

	bin[0] = 0;
	for (int i = maxDeg; i >= 1; i--)
	{
		bin[i] = bin[i - 1];
	}

	//int core = 0;
	int *cNum = new int[n];
	//int *cNum = (int *)malloc(g->n * sizeof(int));
	for (int i = 0; i < n; i++)
	{
		int id = vert[i], nbr, binFrontId;
		//if (i == bin[core + 1]) ++core;
		cNum[id] = tmpDeg[id];
		for (int i = cd[id]; i < cd[id] + deg[id]; i++)
		{
			nbr = adj[i];


			if (tmpDeg[nbr] > tmpDeg[id])
			{
				binFrontId = vert[bin[tmpDeg[nbr]]];
				if (binFrontId != nbr)
				{
					pos[binFrontId] = pos[nbr];
					pos[nbr] = bin[tmpDeg[nbr]];
					vert[bin[tmpDeg[nbr]]] = nbr;
					vert[pos[binFrontId]] = binFrontId;  
				}
				bin[tmpDeg[nbr]]++;
				tmpDeg[nbr]--;
				
			}
			
		}

	}
	coreNum = cNum;
	coreRank = vert;
	delete[] pos;
}

Graph Graph::mksub(Graph &sg, vector<int> nodes)
{
	sg.n = nodes.size();

	for (int i = 0; i < sg.n; i++)
	{
		int ind = cd[nodes[i]];
		for (int j = i+1; j < sg.n; j++)
		{
			while (ind < cd[nodes[i]] + deg[nodes[i]])
			{
				if (adj[ind] == nodes[j])
				{
					sg.e++;
					ind++;
					break;
				}
				else if (adj[ind] > nodes[j]) break;
				ind++;
			}
		}
	}
	sg.edges.resize(sg.e);
	sg.e = 0;
	int *lab = new int[sg.n], cnt = -1;
	for (int i = 0; i < sg.n; i++) lab[nodes[i]] = -1;

	for (int i = 0; i < sg.n; i++)
	{
		int ind = cd[nodes[i]];
		for (int j = i + 1; j < sg.n; j++)
		{
			while (ind < cd[nodes[i]] + deg[nodes[i]])
			{
				if (adj[ind] == nodes[j])
				{
					//cout << "nodes[i] = " << nodes[i] << " nodes[j] = " << nodes[j] << endl;
					lab[nodes[i]] = (lab[nodes[i]] == -1 ? (++cnt) : lab[nodes[i]]);
					lab[adj[ind]] = (lab[adj[ind]] == -1 ? (++cnt) : lab[adj[ind]]);
					sg.edges[sg.e].s = lab[nodes[i]];
					sg.edges[sg.e].t = lab[adj[ind]];
					sg.e++;
					ind++;
					break;
				}
				else if (adj[ind] > nodes[j]) break;
				ind++;
			}
		}
	}

	sg.mkGraph();

	return sg;
}

long long combination(int n, int m)
{
	if (n < m) return 0;
	if (n == m) return 1;
	long long res = 1;
	for (int i = n, j = 1; i >= n - m + 1; i--, j++) res *= i, res /= j;
	return res;
}

int main(int argc, char** argv) {

	Graph g;
	cout << "Reading edgelist from file " << argv[2] << endl;
	g.readedgelist(argv[2]);
	cout << "Reading edgelist finished!" << endl;
	g.mkGraph();

	//printf("LCSize[[]] \n");
	g.coreDecomposition();

	int k = atoi(argv[1]);
	unsigned long long n;

	/*
	ifstream infile("C:\\Users\\Gawssin\\Source\\Repos\\kClist\\kClist\\testOut.txt");

	int tmp;
	for (int i = 0; i < g.n; i++)
	{
		infile >> tmp;
		//fscanf(file, "%d", &tmp);
		if (tmp != g.coreNum[i])
		{
			cout << "error i = " << i << " tmp = " << tmp << " coreNum = " << g.coreNum[i] << endl;
			//printf("error i = %d tmp = %d coreNum = %d!\n", i, tmp, g.coreNum[i]);
			return 0;
		}
	}
	*/
	long long *lowerBound = new long long[g.n];
	for (int i = 0; i < g.n; i++)
	{
		vector<int> iNbr(g.adj+g.cd[i], g.adj + g.cd[i]+g.deg[i]);
		
		Graph sg;
		g.mksub(sg, iNbr);
		sg.coreDecomposition();
		int LCSize = sg.outLargeClique();
		printf("LCSize[%d] = %d\n", i, LCSize);
		lowerBound[i] = combination(LCSize, k - 1);

		printf("lowerBound[%d] = %lld\n\n", i, lowerBound[i]);

		//cout << "lowerBound = "

		/*
		cout << "k = " << k << " sg.n = " << sg.n << " deg = " << g.deg[k] << endl;
		for (int i = 0; i < sg.n; i++)
		{
			cout << "node " << i << " cd " << sg.cd[i] << " deg  " << sg.deg[i] << endl;
			for (int j = sg.cd[i]; j < sg.cd[i]+sg.deg[i]; j++)
			{
				cout << sg.adj[j] << " ";
			}
			cout << endl << endl;
		}

		*/
	}
	

	return 0;

}
