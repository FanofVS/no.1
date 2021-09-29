/**
 * 问题本质为在有向图中寻找负权环
 **/

#include<iostream>
#include<string>
#include<vector>
#include<stack>
#include<queue>
#include<cstring>
#include<cmath>
#include<cstdlib>

#define inf 1e9
#define N 100005

using namespace std;


struct edge {
	int dest;
	double w;
	edge() {}
	edge(int ee, double ww) :dest(ee), w(ww) {}
};


//m, n
int m, n;

// dist (every round)
double dist[N];

// update times (every vertice)
int updates[N];

// in the queue
int inqueue[N];

// graph
vector<edge> V[N];


// divide the graph into several full connected parts
int mark[N];
vector<int> st_ver;
vector<int> part_size;

/*
bool bellman_ford(int v) {
	for (int i = 1; i <= m; i++) {
		dist[i] = inf;
	}
	dist[v] = 0;
}*/

// spfa: average T(n) = O(E)
bool spfa(int v, int con_size) {
	for (int i = 1; i <= m; i++)
		dist[i] = inf;
	dist[v] = 0;

	memset(updates, 0, sizeof(updates));
	memset(inqueue, 0, sizeof(inqueue));

	queue<int> Q;
	Q.push(v);
	inqueue[v] = 1;

	while (!Q.empty()) {
		int v = Q.front();
		Q.pop();
		inqueue[v] = 0;


		int len = V[v].size();

		for (int i = 0; i < len; i++) {
			int u = V[v][i].dest;

			if (dist[v] + V[v][i].w < dist[u]) {
				dist[u] = dist[v] + V[v][i].w;
				updates[u]++;

				// >= num of vertices
				if (updates[u] >= con_size)
					return 1;

				if (inqueue[u] == 0) {
					Q.push(u);
					inqueue[u] = 1;
				}
			}
		}
	}

	return 0;

}

// connected just like non...
int dfs(int v) {
	int len = V[v].size();
	int res = 1;
	mark[v] = 1;
	for (int i = 0; i < len; i++) {
		int dest = V[v][i].dest;
		if (!mark[dest])
			res += dfs(dest);
	}

	return res;
}

void preproc() {
	memset(mark, 0, sizeof(mark));
	for (int i = 1; i <= m; i++) {
		if (mark[i] == 0) {
			st_ver.push_back(i);
			part_size.push_back(dfs(i));
		}
	}
	return;
}

bool check() {
	int sum = 0;
	for (int i = 0; i < part_size.size(); i++)
		sum += part_size[i];
	return sum == m;
}

int main() {
	int x, y, z;
	double val;
	int flg = 0;

	scanf("%d%d", &m, &n);
	for (int i = 0; i < n; i++) {
		scanf("%d%d%d", &x, &y, &z);
		val = log((double)z);

		V[x].push_back(edge(y, val));
		V[y].push_back(edge(x, -val));
	}

	preproc();

	if (!check()) {
		exit(1);
	}

	int len = st_ver.size();

	for (int i = 0; i < len; i++) {
		int temp = st_ver[i];
		int tsiz = part_size[i];
		flg = (flg || spfa(temp, tsiz));

		if (flg) {
			printf("Yes\n");
			return 0;
		}
	}
	printf("No\n");
	return 0;
}