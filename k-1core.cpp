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
#include <time.h>
#include <queue>
#include <cstdarg>
#include <chrono>

#define NLINKS 1000000
using namespace std;
using namespace chrono;

int *label, *updateMark, *deleteArray, *col, leafN, leafM;


int* uadj, *uadjt, *Merge, *mark, *oloc;

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

class Graph
{
public:
	Graph();
	Graph(const Graph& obj);
	~Graph();
	void readedgelist(string edgelist);
	void coreDecomposition();
	void mkGraph();
	int outLargeClique();
	bool isEdge(int, int);
	Graph* mksub(int, ...);
	//Graph* mksubMark(int*, int, int*);
	int color(int*);
	void kClique(int, long long*, long long*, int*, int);
	void kCliqueCount(int, long long*, long long*);
	//bool isEdge(int, int);

	int n;
	int e;
	int maxDeg;
	edge* edges;

	int* deg;
	int* cd;
	int* adj;
	int* coreRank;			//increasing core number order
	int* coreNum;			//coreNum[i] is the core number of node i.
	int* bin;
};

inline int max3(int a, int b, int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

Graph::Graph(void) {}
Graph::~Graph(void)
{
	if (deg != NULL) delete[] deg;
	if (cd != NULL) delete[] cd;
	if (adj != NULL) delete[] adj;
	if (coreRank != NULL) delete[] coreRank;
	if (coreNum != NULL) delete[] coreNum;
	if (bin != NULL) delete[] bin;
}
Graph::Graph(const Graph& obj)
{
	n = obj.n, e = obj.e, maxDeg = obj.maxDeg, edges = obj.edges;
	edges = obj.edges;

	if (deg != NULL) delete[] deg;
	if (obj.deg != NULL) deg = new int[n], memcpy(deg, obj.deg, n * sizeof(int));

	if (cd != NULL) delete[] cd;
	if (obj.cd != NULL) cd = new int[n + 1], memcpy(cd, obj.cd, (n + 1) * sizeof(int));

	if (adj != NULL) delete[] adj;
	if (obj.adj != NULL) adj = new int[2 * e], memcpy(adj, obj.adj, 2 * e * sizeof(int));

	if (coreRank != NULL) delete[] coreRank;
	if (obj.coreRank != NULL) coreRank = new int[n], memcpy(coreRank, obj.coreRank, n * sizeof(int));

	if (coreNum != NULL) delete[] coreNum;
	if (obj.coreNum != NULL) coreNum = new int[n], memcpy(coreNum, obj.coreNum, n * sizeof(int));

	if (bin != NULL) delete[] bin;
	if (obj.bin != NULL) bin = new int[maxDeg + 2], memcpy(bin, obj.bin, (maxDeg + 2) * sizeof(int));
}

void Graph::readedgelist(string edgelist)
{
	ifstream file;
	file.open(edgelist);
	file >> n >> e;
	edges = new edge[e];
	e = 0;
	while (file >> edges[e].s >> edges[e].t) e++;
	file.close();
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
	for (int i = 1; i < n + 1; i++)
	{
		cd[i] = cd[i - 1] + deg[i - 1];
		deg[i - 1] = 0;
	}

	for (int i = 0; i < e; i++)
	{
		adj[cd[edges[i].s] + deg[edges[i].s]++] = edges[i].t;
		adj[cd[edges[i].t] + deg[edges[i].t]++] = edges[i].s;
	}

	for (int i = 0; i < n; i++) sort(adj + cd[i], adj + cd[i] + deg[i]);
}

bool cmp(const pair<int, int>& a, const pair<int, int>& b)
{
	return a.second > b.second;
}
bool IGCmp(const iddeg& a, const iddeg& b)
{
	return a.degree == b.degree ? (a.id < b.id) : (a.degree > b.degree);
}

bool Graph::isEdge(int a, int b)
{
	if (deg[a] > deg[b]) a = a ^ b, b = a ^ b, a = a ^ b;
	for (int i = cd[a]; i < cd[a] + deg[a]; i++)
		if (adj[i] == b) return true;
	return false;
}
int Graph::outLargeClique()
{
	int CSize = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		int id = coreRank[i];
		if (coreNum[id] >= CSize)
		{
			pair<int, int>* SCore = new pair<int, int>[deg[id]];
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
			int* C = new int[deg[id]];
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
	int* vert = new int[n](), *pos = new int[n](), *tmpDeg = new int[n]();
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
	int* cNum = new int[n];
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


/*
Graph* Graph::mksubMark(int* nodes, int NodeNum, int *mark)
{
	Graph* sg = new Graph;
	sg->n = NodeNum, sg->e = 0;
	int* newFg = new int[n]();

	for (int i = 0; i < NodeNum; i++) newFg[nodes[i]] = 1;

	for (int i = 0; i < e; i++)
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1) sg->e++;
	sg->edges.resize(sg->e);

	sg->e = 0;
	int* lab = new int[n], cnt = 0;
	//for (int i = 0; i < sg.n; i++) lab[nodes[i]] = -1;

	for (int i = 0; i < sg->n; i++)
	{
		lab[nodes[i]] = cnt;
		mark[cnt++] = nodes[i];
	}

	for (int i = 0; i < e; i++)
	{
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1)
		{
			//lab[edges[i].s] = (lab[edges[i].s] == -1 ? (cnt++) : lab[edges[i].s]);
			//lab[edges[i].t] = (lab[edges[i].t] == -1 ? (cnt++) : lab[edges[i].t]);
			//mark[lab[edges[i].s]] = edges[i].s;
			//mark[lab[edges[i].t]] = edges[i].t;
			sg->edges[sg->e].s = lab[edges[i].s];
			sg->edges[sg->e].t = lab[edges[i].t];
			sg->e++;
		}
	}

	delete[] newFg;
	delete[] lab;
	sg->mkGraph();
	return sg;
}
*/




int Graph::color(int* color)
{
	iddeg* ig = new iddeg[n];
	for (int i = 0; i < n; i++)
	{
		ig[i].id = i;
		ig[i].degree = deg[i];
	}

	sort(ig, ig + n, IGCmp);

	//color = new int[n];
	memset(color, -1, sizeof(int) * n);
	int* C = new int[(ig[0].degree + 1)]();

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


double combination(int n, int m)
{
	if (n < m) return 0;
	if (n == m) return 1;
	double res = 1;
	for (int i = n, j = 1; i >= n - m + 1; i--, j++) res *= i, res /= j;
	return res;
}

int outList1(int* list1, int* i, int* j, int s, int s1)
{
	if (*i == s - 1)
		return list1[++ * j];
	if (*j == s1 - 1)
		return list1[++ * i];
	if (list1[*i + 1] < list1[*j + 1])
		return list1[++ * i];
	return list1[++ * j];
}
int merging(int s, int* list1, int s1, int* list2, int s2, int* list3)
{
	int i = 0, j = 0, p = -1, q = s - 1, s3 = 0;
	int x = outList1(list1, &p, &q, s, s1), y = list2[0];
	while (i < s1 && j < s2)
	{
		if (x < y)
		{
			x = outList1(list1, &p, &q, s, s1);
			++i;
			//x = list1[++i];
			continue;
		}
		if (y < x)
		{
			y = list2[++j];
			continue;
		}
		list3[s3++] = x;
		x = outList1(list1, &p, &q, s, s1);
		++i;
		//x = list1[++i];
		y = list2[++j];
	}
	return s3;
}



int *lab, *cdv, *adjv, *ns, ** degS, ** subS, *ver, *dTMP, K;


void init(Graph *g, int k)
{
	leafN = g->n;
	leafM = g->e;
	K = k;
	
	label = new int[g->n]();
	updateMark = new int[g->n]();
	col = new int[g->n];

	dTMP = new int[g->n];
	cdv = new int[g->n + 1];
	lab = new int[g->n]();
	adjv = new int[g->e * 2];
	ns = new int[k + 1];
	degS = new int*[k + 1], subS = new int*[k + 1];
	for (int i = 2; i <= k; i++)
	{
		degS[i] = new int[g->n];
		subS[i] = new int[g->n];
	}
	ver = new int[k + 1];
}
void deleteSubG(int k)
{
	delete[] lab, delete[] cdv, delete[] adjv, delete[] ns;
	for (int i = 2; i <= k; i++)
	{
		delete[] degS[i], delete[] subS[i];
	}
	delete[] degS, delete[] subS;
}

void Graph::kCliqueCount(int l, long long* tol,long long* cnt)
{
	int u, v, w, end;
	if (l == 2)
	{
		//printf("l = 2\n");
		int k = K;
		for (int i = 0; i < ns[2]; i++)
		{//list all edges
			u = subS[2][i];
			ver[2] = u;
			//(*n)+=g->d[2][u];
			end = cdv[u] + degS[2][u];
			for (int p = 2; p <= k; p++)
				cnt[ver[p]] += degS[2][u];

			(*tol) += degS[2][u];

			for (int j = cdv[u]; j < end; j++)
			{
				//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
				cnt[adjv[j]]++;
			}

		}
		return;
	}

	for (int i = 0; i < ns[l]; i++)
	{
		u = subS[l][i];
		ver[l] = u;
		//printf("%d %d\n",i,u);
		ns[l - 1] = 0;
		end = cdv[u] + degS[l][u];
		//printf("dd degS[l][u] = %d\n", degS[l][u]);
		for (int j = cdv[u]; j < end; j++) //relabeling nodes and forming U'.
		{
			
			v = adjv[j];
			//printf("ddv = %d \n", v);
			if (lab[v] == l) {
				lab[v] = l - 1;
				subS[l - 1][ns[l - 1]++] = v;
				degS[l - 1][v] = 0;//new degrees
			}
		}
		//printf("dd %d %d\n", i, u);
		for (int j = 0; j < ns[l - 1]; j++) //reodering adjacency list and computing new degrees
		{
			v = subS[l - 1][j];
			end = cdv[v] + degS[l][v];
			for (int k = cdv[v]; k < end; k++)
			{
				w = adjv[k];
				if (lab[w] == l - 1)
					degS[l - 1][v]++;
				else
				{
					adjv[k--] = adjv[--end];
					adjv[end] = w;
				}
			}
		}
		//printf("cc %d %d\n", i, u);
		kCliqueCount(l - 1, tol, cnt);

		for (int j = 0; j < ns[l - 1]; j++) {//restoring labels
			v = subS[l - 1][j];
			lab[v] = l;
		}

	}


}




void Graph::kClique(int k, long long* tol, long long* cnt, int *subg, int size)
{
	int cdTmp = 0, nsg = 0;
	for (int i = 0; i < size; i++) lab[subg[i]] = k;
	for (int i = 0; i < size; i++)
	{
		int u = subg[i];
		degS[k][u] = 0;
		for (int j = cd[u]; j < cd[u] + deg[u]; j++)
		{
			int v = adj[j];
			if (lab[v] == k) degS[k][u]++;
		}
	}
	//printf("fff\n");
	for (int i = 0; i < size; i++)
	{
		int u = subg[i];
		dTMP[u] = 0;
		cdv[u] = cd[u];
		subS[k][nsg++] = u;
		for (int j = cd[u]; j < cd[u] + deg[u]; j++)
		{
			int v = adj[j];
			if (lab[v] == k)
			{
				if (degS[k][u] < degS[k][v] || (degS[k][u] == degS[k][v] && u < v))
				{
					adjv[cdv[u] + dTMP[u]++] = v;
				}
			}
		}
	}
	//printf("ggg\n");
	for (int i = 0; i < size; i++) degS[k][subg[i]] = dTMP[subg[i]];

	//for (int i = 0; i < size; i++) printf("deg = %d\n", d[subg[i]]);

	ns[k] = nsg;
	//degS[k] = d;
	//subS[k] = sub;
	//printf("nsg = %d\n", nsg);


	

	kCliqueCount(k, tol, cnt);



	

	//printf("eee\n");

	for (int i = 0; i < size; i++) lab[subg[i]] = 0;


	/*

	// ord_core
	coreDecomposition();

	// relabel
	int sCore, tCore;
	for (int i = 0; i < e; i++)  
	{
		if (coreNum[edges[i].s] > coreNum[edges[i].t])
		{
			edges[i].s = edges[i].s ^ edges[i].t;
			edges[i].t = edges[i].s ^ edges[i].t;
			edges[i].s = edges[i].s ^ edges[i].t;
		}
	}

	// mkspecial
	
	//init d!!!!

	int nsg, maxDv;
	
	for (int i = 0; i < e; i++)	d[edges[i].s]++;

	
	nsg = 0, maxDv = 0, cdv[0] = 0;
	

	for (int i = 1; i < n + 1; i++)
	{
		cdv[i] = cdv[i - 1] + d[i - 1];
		maxDv = (maxDv > d[i - 1]) ? maxDv : d[i - 1];
		sub[nsg++] = i - 1;
		d[i - 1] = 0;
		lab[i - 1] = k;
	}

	
	for (int i = 0; i < e; i++)
		adjv[cdv[edges[i].s] + d[edges[i].s]++] = edges[i].t;

	
	ns[k] = nsg;
	degS[k] = d;
	subS[k] = sub;
	kCliqueCount(k, tol, ver, lab, cdv, adjv, ns, degS, subS, cnt);
	*/
}





long long Binary(double left, double right, Graph& g, double* uBd)
{
	//int left = 0, right = n;      //解的范围初始为(0,n],不包含0                             
	while (left + 1 < right)
	{
		long long mid = (right + left) / 2;
		if (1)//check()
		{
			right = mid;          //修正解的范围(left,right]
		}
		else
			left = mid;
	}
	return right;                //最后left + 1 = right
}




void deleteNodes(Graph *g, int *delArray, int size)
{
	leafN -= size;

	if (size == 1)
	{
		int u = delArray[0];
		leafM -= g->deg[u];
		for (int j = g->cd[u]; j < g->cd[u] + g->deg[u]; j++)
		{
			int v = g->adj[j];
			for (int k = g->cd[v]; k < g->cd[v] + g->deg[v]; k++)
			{
				int w = g->adj[k];
				if (w == u)
				{
					g->adj[k] = g->adj[g->cd[v] + g->deg[v] - 1];
					g->adj[g->cd[v] + g->deg[v] - 1] = w;
					g->deg[v]--;
					break;
				}
			}
		}
		g->deg[u] = 0;
		label[u] = -1;
		return;
	}




	for (int i = 0; i < size; i++)
	{
		int u = delArray[i];
		label[u] = -1;
		for (int j = g->cd[u]; j < g->cd[u] + g->deg[u]; j++) updateMark[g->adj[j]] = 1;
	}

	g->maxDeg = 0;
	for (int i = 0; i < g->n; i++)
	{
		if (label[i] == -1) g->deg[i] = 0;
		else
		{
			if (updateMark[i])
			{
				for (int j = g->cd[i]; j < g->cd[i] + g->deg[i]; j++)
				{
					int u = g->adj[j];
					if (label[u] == -1)
					{
						g->adj[j] = g->adj[g->cd[i] + g->deg[i] - 1];
						g->adj[g->cd[i] + g->deg[i] - 1] = u;
						g->deg[i]--;
						j--;
					}
				}

			}
			updateMark[i] = 0;
			g->maxDeg = max(g->deg[i], g->maxDeg);

			//g->maxDeg = max(g->deg[i], g->maxDeg);
		}
	}

	/*
	for (int i = 0; i < size; i++)
	{
		int u = g->coreRank[i];
		for (int j = g->cd[u]; j < g->cd[u] + g->deg[u]; j++)
		{
			int v = g->adj[j];
			for (int k = g->cd[v]; k < g->cd[v] + g->deg[v]; k++)
			{
				int w = g->adj[k];
				if (w == u)
				{
					g->adj[k] = g->adj[g->cd[v] + g->deg[v] - 1];
					g->adj[g->cd[v] + g->deg[v] - 1] = w;
					g->deg[v]--;
					break;
				}
			}
		}
		g->deg[u] = 0;
	}
	*/



}

void dpDecp(Graph *g, int k, double lwb)
{

	int *delNodes = new int[g->n], delN = 0;
	double minDPval;
	while (true)
	{
		int tolcol = g->color(col);

		int *colCnt = new int[tolcol]();
		double **ck;
		ck = new double*[tolcol + 1];

		for (int i = 0; i < tolcol + 1; i++)
		{
			ck[i] = new double[k + 1];
			ck[i][0] = 1;
		}

		for (int i = 1; i < k + 1; i++) ck[0][i] = 0;

		minDPval = 1e300;
		for (int i = 0; i < g->n; i++)
		{
			if (label[i] != -1)
			{
				for (int j = g->cd[i]; j < g->cd[i] + g->deg[i]; j++)
				{
					int u = g->adj[j];
					colCnt[col[u]]++;
				}

				for (int q = 1; q < tolcol + 1; q++)
					for (int p = 1; p < k + 1; p++)
						ck[q][p] = ck[q - 1][p - 1] * colCnt[q - 1] + ck[q - 1][p];

				
				minDPval = ck[tolcol][k] < minDPval ? ck[tolcol][k] : minDPval;

				if (ck[tolcol][k] < lwb) delNodes[delN++] = i;

				for (int j = g->cd[i]; j < g->cd[i] + g->deg[i]; j++) colCnt[col[g->adj[j]]] = 0;

			}

		}

		delete[] colCnt;

		for (int i = 0; i < tolcol + 1; i++) delete[] ck[i];
		delete[] ck;

		if (delN == 0) break;
		else deleteNodes(g, delNodes, delN);
		delN = 0;
	}

	leafM = 0;
	for (int i = 0; i < g->n; i++)
		if (label[i] != -1) leafM += g->deg[i];

	leafM /= 2;

	cout << fixed;
	cout << "minDPval: " << minDPval / leafN << endl;
	cout << "colorful k-star core density: " << 1.0 * leafM / leafN << endl;

	delete[] delNodes;
}

void peeling(Graph *g, int k)
{
	long long tol = 0, tolCliques = 0, *cnt = new long long[g->n](), *subcnt = new long long[g->n]();
	int *subG = new int[leafN], subN = 0, lessN;
	for (int i = 0; i < g->n; i++)
		if (label[i] != -1) subG[subN++] = i;

	lessN = subN;
	g->kClique(k, &tol, cnt, subG, subN);
	tolCliques = tol;

	//for (int i = 0; i < g->n; i++)
	//	if (label[i] != -1) printf("i = %d cliques = %lld\n",i, cnt[i]);


	priority_queue< pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>> > q;

	for (int i = 0; i < subN; i++) q.push(make_pair(cnt[subG[i]], subG[i]));
	
	
	int counter = 0;
	long long maxKCC = 0;
	double kCliqueDensity = 0, density = 0;
	while (!q.empty())
	{
		int c = q.top().second;
		long long coreNum = q.top().first;

		if (coreNum >= maxKCC)
		{
			maxKCC = coreNum;
			if (1.0*tolCliques / lessN > kCliqueDensity)
			{
				kCliqueDensity = 1.0*tolCliques / lessN;
				density = 1.0 * leafM / lessN;
			}
		}

		q.pop();
		if (label[c] == -1) continue;
		
		//printf("coreNum = %lld tolCliques = %lld\n", coreNum, tolCliques);

		tolCliques -= coreNum;



		

		//printf("c = %d\n", c);
		for (int i = g->cd[c]; i < g->cd[c] + g->deg[c]; i++) label[g->adj[i]] = 1;
		//printf("cc = %d\n", g->deg[c]);


		for (int i = g->cd[c]; i < g->cd[c] + g->deg[c]; i++)
		{
			int v = g->adj[i];
			if (label[v] == -1) continue; //

			subN = 0, tol = 0;
			for (int j = g->cd[v]; j < g->cd[v] + g->deg[v]; j++)
			{
				int w = g->adj[j];
				if (label[w] == 1) subG[subN++] = w;
			}
			//printf("v = %d subN = %d\n", v, subN);


			//printf("a******a\n");



			g->kClique(k - 2, &tol, subcnt, subG, subN);





			//printf("tol = %lld\n", tol);
			cnt[v] -= tol;
			q.push(make_pair(cnt[v], v));
		}

		

		for (int i = g->cd[c]; i < g->cd[c] + g->deg[c]; i++) label[g->adj[i]] = 0;

		deleteNodes(g, &c, 1);
		lessN--;
		//printf("lessN = %d\n", lessN);
		if (lessN == 0) break;

		
	}
	cout << fixed;
	cout << "kCliqueDensity: " << kCliqueDensity << endl;
	cout << "density: " << density << endl;
}





int main(int argc, char** argv)
{
	Graph *g = new Graph;
	int k = atoi(argv[1]);
	cout << "Reading edgelist from file " << argv[2] << endl;
	g->readedgelist(argv[2]);
	cout << "Reading edgelist finished!" << endl;

	cout << "N: " << g->n << endl << "M: " << g->e << endl;


	high_resolution_clock::time_point kccDP_s = high_resolution_clock::now();

	g->mkGraph();
	cout << "mkGraph finished!" << endl;

	init(g, k);
	
	g->coreDecomposition();
	for (int Core_ind = 0; Core_ind < g->n; Core_ind++)
	{
		if (g->coreNum[g->coreRank[Core_ind]] >= k - 1)
		{
			deleteNodes(g, g->coreRank, Core_ind);
			printf("Core_ind = %d!\n",Core_ind);
			break;
		}
	}

	//double lwb = combination(wCoreN-1, k-1);
	cout << "k-1 core: " << leafN << endl;

	//dpDecp(g, k, lwb);

	//cout << "colorful k-star core: " << leafN << endl;



	leafM = 0;
	for (int i = 0; i < g->n; i++)
		if (label[i] != -1) leafM += g->deg[i];
	leafM /= 2;

	peeling(g, k);


	deleteSubG(k);


	high_resolution_clock::time_point kccDP_e = high_resolution_clock::now();
	auto kccDP = std::chrono::duration_cast<std::chrono::microseconds>(kccDP_e - kccDP_s).count();

	cout << "toltal time: " << kccDP/1e6 << endl;


	return 0;

}
