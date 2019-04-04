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

#define NLINKS 1000000
using namespace std;

class edge
{
public:
	int s;
	int t;
};

class iddeg
{
public:
	int id;
	int degree;
};

class Graph{
public:
	Graph();
	Graph(const Graph &obj);
	~Graph();
	void readedgelist(string edgelist);
	void coreDecomposition();
	void mkGraph();
	int outLargeClique();
	bool isEdge(int, int);
	void mksub(Graph &, int *, int);
	int color(int *);
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

Graph::Graph(void){}
Graph::~Graph(void)
{
	if (deg != NULL) delete[] deg;
	if (cd != NULL) delete[] cd;
	if (adj != NULL) delete[] adj;
	if (coreRank != NULL) delete[] coreRank;
	if (coreNum != NULL) delete[] coreNum;
	if (bin != NULL) delete[] bin;
}
Graph::Graph(const Graph &obj)
{
	n = obj.n, e = obj.e, maxDeg = obj.maxDeg, edges = obj.edges;
	edges = obj.edges;

	if (deg != NULL) delete[] deg;
	if (obj.deg != NULL) deg = new int[n], memcpy(deg, obj.deg, n * sizeof(int));

	if (cd != NULL) delete[] cd;
	if (obj.cd != NULL) cd = new int[n], memcpy(cd, obj.cd, n * sizeof(int));

	if (adj != NULL) delete[] adj;
	if (obj.adj != NULL) adj = new int[2*e], memcpy(adj, obj.adj, 2*e* sizeof(int));

	if (coreRank != NULL) delete[] coreRank;
	if (obj.coreRank != NULL) coreRank = new int[n], memcpy(coreRank, obj.coreRank, n * sizeof(int));

	if (coreNum != NULL) delete[] coreNum;
	if (obj.coreNum != NULL) coreNum = new int[n], memcpy(coreNum, obj.coreNum, n * sizeof(int));

	if (bin != NULL) delete[] bin;
	if (obj.bin != NULL) bin = new int[maxDeg + 2], memcpy(bin, obj.bin, (maxDeg + 2) * sizeof(int));
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
	deg = new int[n]();
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
bool IGCmp(const iddeg &a, const iddeg &b)
{
	return a.degree == b.degree ? (a.id < b.id) : (a.degree > b.degree);
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

	for (int i = 0; i < n; i++) 
		bin[deg[i]]++;

	int lastBin = bin[0], nowBin;
	bin[0] = 0;
	for (int i = 1; i <= maxDeg; i++)
	{
		nowBin = lastBin + bin[i - 1];
		lastBin = bin[i];
		bin[i] = nowBin;
	}
	int *vert = new int[n](), *pos = new int[n](), *tmpDeg = new int[n]();
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

	delete[] tmpDeg;
	delete[] pos;
}

void Graph::mksub(Graph &sg, int *nodes, int NodeNum)
{
	sg.n = NodeNum, sg.e = 0;
	int * newFg = new int[n]();

	for (int i = 0; i < NodeNum; i++) newFg[nodes[i]] = 1;

	for (int i = 0; i < e; i++)
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1) sg.e++;

	//printf("sg.e = %d\n", sg.e);
	sg.edges.resize(sg.e);



	//sort(nodes, nodes+NodeNum);
	

	/*
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
	printf("sg.e = %d\n", sg.e);
	sg.edges.resize(sg.e);
	*/
	
	sg.e = 0;
	int *lab = new int[n], cnt = 0;
	for (int i = 0; i < sg.n; i++) lab[nodes[i]] = -1;


	for (int i = 0; i < e; i++)
	{
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1)
		{
			lab[edges[i].s] = (lab[edges[i].s] == -1 ? (cnt++) : lab[edges[i].s]);
			lab[edges[i].t] = (lab[edges[i].t] == -1 ? (cnt++) : lab[edges[i].t]);
			sg.edges[sg.e].s = lab[edges[i].s];
			sg.edges[sg.e].t = lab[edges[i].t];
			sg.e++;
		}
	}


	/*
	
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
	
	*/
	//printf("sg labeled\n");

	delete[] newFg;
	delete[] lab;
	sg.mkGraph();
}

int Graph::color(int *color)
{
	iddeg *ig = new iddeg[n];
	for (int i = 0; i < n; i++)
	{
		ig[i].id = i;
		ig[i].degree = deg[i];
	}

	sort(ig,ig+n,IGCmp);

	//color = new int[n];
	memset(color, -1, sizeof(int)*n);
	int *C = new int[(ig[0].degree + 1)]();

	color[ig[0].id] = 0;
	int colorNum = 1;

	for (int i = 1; i < n; i++)
	{
		int tmpDeg = ig[i].degree, tmpid = ig[i].id;
		for (int j = 0; j < tmpDeg; j++)
		{
			int now = adj[cd[tmpid] + j];
			if (color[now] != -1)
				C[color[now]] = 1;
		}
		for (int j = 0; j < ig[0].degree + 1; j++)
			if (C[j] == 0)
			{
				color[ig[i].id] = j;
				colorNum = j > colorNum ? j : colorNum;
				break;
			}

		for (int j = 0; j < tmpDeg; j++)
		{
			int now = adj[cd[tmpid] + j];
			if (color[now] != -1)
				C[color[now]] = 0;
		}

	}
	//printf("color number = %d\n", colorNum);
	delete[] ig;
	delete[] C;
	return colorNum;
}


long long combination(int n, int m)
{
	if (n < m) return 0;
	if (n == m) return 1;
	long long res = 1;
	for (int i = n, j = 1; i >= n - m + 1; i--, j++) res *= i, res /= j;
	return res;
}

int main(int argc, char** argv) 
{
	Graph g;
	int k = atoi(argv[1]);
	cout << "Reading edgelist from file " << argv[2] << endl;
	g.readedgelist(argv[2]);
	cout << "Reading edgelist finished!" << endl;
	g.mkGraph();
	cout << "mkGraph finished!" << endl;

	g.coreDecomposition();
	cout << "coreDecomposition finished!" << endl;


	if (g.coreNum[g.coreRank[g.n-1]] < k - 1)
	{
		cout << "No " << k << "-clique core!" << endl;
		return 0;
	}


	Graph cg;
	int wCoreN = g.outLargeClique();
	for (int Core_ind = 0; Core_ind < g.n; Core_ind++)	//first node with k-clique core number >= k-1
	{
		if (g.coreNum[g.coreRank[Core_ind]] >= wCoreN - 1)
		{

			int *CGNodes = new int[g.n - Core_ind];
			printf("g.n - Core_ind = %d!\n", g.n - Core_ind);
			memcpy(CGNodes, g.coreRank + Core_ind, sizeof(int) * (g.n - Core_ind));

			g.mksub(cg, CGNodes, g.n - Core_ind);
			delete[] CGNodes;
			break;
		}
	}

	//long long n;


	long long *lowerBound = new long long[cg.n];
	printf("\nstart cg.n = %d!\n\n", cg.n);



	// lowerBound
	/*
	for (int i = 0; i < cg.n; i++)
	{

		Graph sg;
		
		int *iNbr = new int[cg.deg[i]];
		memcpy(iNbr, cg.adj + cg.cd[i], sizeof(int) * cg.deg[i]);

		cg.mksub(sg, iNbr, cg.deg[i]);

		sg.coreDecomposition();
		int LCSize = sg.outLargeClique();
		//printf("LCSize[%d] = %d\n", i, LCSize);
		lowerBound[i] = combination(LCSize, k - 1);

		//printf("lowerBound[%d] = %lld\n\n", i, lowerBound[i]);

	}
	*/


	
	//Graph oG = cg;
	Graph *newG, *oldG = new Graph(cg);

	//*oldG = cg;

	oldG->coreDecomposition();
	int w = oldG->outLargeClique();
	int lwb = combination(w - 1, k - 1);



	double **combin = new double *[k];
	int maxD = oldG->maxDeg;
	for (int i = 0; i < k; i++)
		combin[i] = new double[maxD]();
	for (int i = 0; i < maxD; i++) combin[0][i] = 1.0;
	for (int i = 0; i < k; i++)	combin[i][i] = 1.0;

	for (int i = 2; i < maxD; i++)
		for (int j = 1; j < i && j < k; j++)
			combin[j][i] = combin[j][i - 1] + combin[j - 1][i - 1];


	printf("combin[k - 1][maxD - 1] = %lf\n", combin[k - 1][maxD - 1]);

	while (true)
	{	
		
		
		int *col = new int[oldG->n], newN = 0;
		int colNum = oldG->color(col);
		int *newND = new int[oldG->n];
		//delete[] col;
		int *C = new int[oldG->maxDeg + 1]();
		
		int t1 = 0, t2 = 0;
		for (int i = 0; i < oldG->n; i++)
		{
			int cv = 0, maxCN = 0, secCN = 0, tmpCN, wCol = -1;
			for (int j = oldG->cd[i]; j < oldG->cd[i] + oldG->deg[i]; j++)
			{
				if (C[col[oldG->adj[j]]] == 0) cv++;
				C[col[oldG->adj[j]]]++;
				
				tmpCN = C[col[oldG->adj[j]]];

				if (tmpCN > maxCN)
				{
					if (col[oldG->adj[j]] != wCol)
					{
						secCN = maxCN;
						wCol = col[oldG->adj[j]];
					}
					maxCN = tmpCN;
					continue;
				}

				secCN = secCN > tmpCN ? secCN : tmpCN;
				
			}



			double uLimit = 0;
			if (cv == 1)
				if (oldG->deg[i] >= w - 1)
					newND[newN++] = i;
			if (cv == 2)
				uLimit = maxCN * combin[k - 2][oldG->deg[i] - maxCN] + combin[k - 1][oldG->deg[i] - maxCN];

			if (cv >= 3)
				uLimit = maxCN * secCN * combin[k - 3][oldG->deg[i] - maxCN - secCN]
				+ (maxCN + secCN) * combin[k - 2][oldG->deg[i] - maxCN - secCN]
				+ combin[k - 1][oldG->deg[i] - maxCN - secCN];

			if (uLimit < 0) 
				printf("over flow!\n");


			if (uLimit >= lwb) newND[newN++] = i;



			for (int j = oldG->cd[i]; j < oldG->cd[i] + oldG->deg[i]; j++)
				C[col[oldG->adj[j]]] = 0;
		}

		//printf("t1 = %d t2 = %d\n naive upperBound there newN = %d oldG.n = %d\n", t1, t2, newN, oldG->n);


		if (newN == oldG->n)
		{
			printf("start dp there\n");
			/*
			ofstream outFile;
			outFile.open("C:\\Users\\Gawssin\\Source\\Repos\\Code\\kCliqueListing\\outCoreIn.txt");
			for (int i = 0; i < oldG->e; i++)
			{
				outFile << oldG->edges[i].s << " " << oldG->edges[i].t << endl;

			}
			outFile.close();
			*/

			newN = 0;
			
			for (int i = 0; i < oldG->n; i++)
			{
				//printf("i = %d\n", i);
				int colN = 0;
				for (int j = oldG->cd[i]; j < oldG->cd[i] + oldG->deg[i]; j++)
				{
					if (C[col[oldG->adj[j]]] == 0) colN++;
					C[col[oldG->adj[j]]]++;
				}

				if (colN >= k - 1)
				//if (colN == k - 1) upperBound = 1;
				{
					int *nbrCol = new int[colN];
					colN = 0;

					for (int j = 0; j < oldG->maxDeg + 1; j++)
						if(C[j]) nbrCol[colN++] = C[j];


					double **dpCol = new double *[colN + 1];
					for (int j = 0; j < colN + 1; j++)
						dpCol[j] = new double[k];

					for (int j = 0; j < colN + 1; j++)	dpCol[j][0] = 1;
					for (int j = 1; j < k; j++)			dpCol[j][j] = nbrCol[j - 1] * dpCol[j - 1][j - 1];
					for (int j = 2; j < colN + 1; j++)
						for (int p = 1; p < j && p < k; p++)
							dpCol[j][p] = nbrCol[j - 1] * dpCol[j - 1][p - 1] + dpCol[j - 1][p];

					
					if (dpCol[colN][k - 1] < 0)
						printf("dpCol over flow\n");


					if (dpCol[colN][k - 1] >= lwb*1.0)
						newND[newN++] = i;

					/*
					if ((colN + k - 2 < w - 1) && dpCol[colN][k - 1] >= lwb)
					{
						printf("pause!\n");
					}
					*/



					delete[] nbrCol;

					for (int j = 0; j < colN + 1; j++)
						delete[] dpCol[j];
					delete[] dpCol;
				}
				for (int j = oldG->cd[i]; j < oldG->cd[i] + oldG->deg[i]; j++)
					C[col[oldG->adj[j]]] = 0;

			}
			
		}
			
		delete[] C;
		if (newN < oldG->n)
		{
			newG = new Graph;
			oldG->mksub(*newG,newND,newN);
			delete oldG;
			oldG = newG;
		}
		else break;
	}













	int minDeg = 1000000;
	for (int i = 0; i < oldG->n; i++) minDeg = minDeg < oldG->deg[i] ? minDeg : oldG->deg[i];


	/*
	ofstream outFile;
	outFile.open("C:\\Users\\Gawssin\\Source\\Repos\\Code\\kCliqueListing\\outCoreInW.txt");
	
	for (int i = 0; i < oldG->e; i++)
	{
		outFile << oldG->edges[i].s << " " << oldG->edges[i].t << endl;

	}

	outFile.close();
	*/


	printf("\nNow G.n = %d minDeg = %d \n",oldG->n, minDeg);



	return 0;

}
