// 杂项代码文件.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

//某题动规的一维实现
/*
#include<iostream>
#include<algorithm>
#include<cmath>
#include<cstring>


using namespace std;

int f[1005];

int main(){
	char s[1005], t[1005];
	int sum;
	int temp;
	cin >> sum;
	for (int i = 0; i < sum; i++){
		cin >> s >> t;
		memset(f, 0, sizeof(f));

		int x = strlen(s), y = strlen(t);
		int p = (x>y ? x : y);

		while (p){
			f[p] = p;
			p--;
		}

		for (int j = 1; j <= x; j++){
			f[0] = temp = j - 1;
			for (int k = 1; k <= y; k++){
				if (s[j - 1] == t[k - 1]){
					int c = temp;
					temp = f[k];
					f[k] = c;
				}
				else{
					int c = temp < f[k - 1] ? temp : f[k - 1];
					int g = f[k];
					f[k] = f[k] < c ? f[k] : c;
					f[k]++;
					temp = g;
				}
				
			}
		}
		cout << f[y] << endl;
	}
	return 0;

}
*/

//popular cows
/*
#include<iostream>
#include<algorithm>
#include<vector>
#include<stack>
#include<queue>
#include<set>
#include<memory.h>
#define N 10005


using namespace std;

struct node{
	vector<int> likes;
	int visited;
	int instack;
	int dfn;
	int low;
	int myc;
}cow[N];

struct color{
	int vis ;
	int sum ;
	set<int> outs;
}ele[N];


stack<int> branch;
stack<int> forvis;
stack<int> forser;

int m, n;
int ind = 0;
int bloc = 1;
int plac = 0;


void tar(int x){
	int i, j, k;
	cow[x].dfn = cow[x].low = ind++;
	cow[x].visited =cow[x].instack = 1;

	//用于判断图的连通

	branch.push(x);
	int siz = cow[x].likes.size();
	for (i = 0; i < siz; i++){
		j = cow[x].likes[i];
		if (!cow[j].visited)
		{
			tar(j);
			if (cow[x].low>cow[j].low)
				cow[x].low = cow[j].low;
		}
		if (cow[j].instack){
			if (cow[x].low>cow[j].low)
				cow[x].low = cow[j].low;
		}
	}

	if (cow[x].dfn == cow[x].low){
		k = 0;
		ele[k].sum = ele[k].vis = 0;
		forser.push(bloc);
		while (k != x){
			k = branch.top();
			branch.pop();
			cow[k].instack = 0;

			//块号
			cow[k].myc = bloc;

			//块大小
			ele[bloc].sum++;

			//块出边
			int pl = cow[k].likes.size();
			for (i = 0; i < pl; i++)
			{
				int j = cow[k].likes[i];
		
				ele[bloc].outs.insert(j);
			}
		}
		bloc++;
	}
}

void ser(int x){
	int j, k = 0;
	ele[x].vis = 1;
	set<int>::iterator i;
	for (i = ele[x].outs.begin(); i != ele[x].outs.end(); i++){
		j = cow[*i].myc;
		if (j != x)
			k = 1;
		if (!ele[j].vis)
			ser(j);
	}
	if (!k){
		if (!plac)
			plac = x;
		else
			plac = N;
	}
}


int main(){
	int i, j, k, p, q;
	memset(cow, 0, sizeof(cow));
	scanf("%d%d", &n, &m);

	for (i = 1; i <= n; i++)
		forvis.push(i);

	for (i = 0; i < m; i++){

		scanf("%d%d", &p, &q);
		cow[p].likes.push_back(q);

	}

	while (!forvis.empty())
	{
		p = forvis.top();
		forvis.pop();
		if (!cow[p].visited)
			tar(p);
	}

	while (!forser.empty()){
		p = forser.top();
		forser.pop();
		if (!ele[p].vis)
			ser(p);
	}


	if (plac <= n)
		p = ele[plac].sum;
	else
		p = 0;
	printf("%d\n", p);
	return 0;



}

*/

//going from u to v or from v to u(照搬上一题)
/*
#include<iostream>
#include<algorithm>
#include<vector>
#include<stack>
#include<queue>
#include<set>
#include<memory.h>
#define N 1005


using namespace std;

struct node{
	vector<int> likes;
	int visited;
	int instack;
	int dfn;
	int low;

}cow[N];



stack<int> branch;



int m, n;
int ind = 0;

int plac = 0;


void tar(int x){
	int i, j, k;
	cow[x].dfn = cow[x].low = ind++;
	cow[x].visited = cow[x].instack = 1;

	//用于判断图的连通

	branch.push(x);
	int siz = cow[x].likes.size();
	for (i = 0; i < siz; i++){
		j = cow[x].likes[i];
		if (!cow[j].visited)
		{
			tar(j);
			if (cow[x].low>cow[j].low)
				cow[x].low = cow[j].low;
		}
		if (cow[j].instack){
			if (cow[x].low>cow[j].low)
				cow[x].low = cow[j].low;
		}
	}

	if (cow[x].dfn == cow[x].low){
		k = 0;
		
		while (k != x){
			k = branch.top();
			branch.pop();
			cow[k].instack = 0;
		}
		if (x == 1)
			plac = 1;
	}
}



int main(){
	int i, j, k, p, q, t;
	scanf("%d", &t);
	for (j = 0; j < t; j++){
		plac = 0;
		memset(cow, 0, sizeof(cow));
		scanf("%d%d", &n, &m);

		for (i = 0; i < m; i++){

			scanf("%d%d", &p, &q);
			cow[p].likes.push_back(q);

		}

		tar(1);

		if (plac)
			printf("Yes\n");
		else
			printf("No\n");

	}

	
	return 0;



}
*/

//Caocao's bridges
//不连通输出0
/*
#include<iostream>
#include<algorithm>
#include<stack>
#include<queue>
#include<memory.h>
#include<vector>
#define N 1005
#define MA 100000

using namespace std;


struct node{
	vector<int> lik;
	int dfn;
	int low;
	int visited;
	int par;
}isl[N];

int n, m, u, v, w;

int cou;
int bc[N][N];
int gua[N][N];
stack<int> fos;



void tar(int x){
	int i, j, k;
	isl[x].visited = 1;

	isl[x].dfn = cou;
	isl[x].low = cou++;
	
	k = isl[x].lik.size();
	for (i = 0; i < k; i++){
		j = isl[x].lik[i];
		
		
		if (!isl[j].visited)
		{
			isl[j].par = x;
			tar(j);
			if (isl[x].low>isl[j].low)
				isl[x].low = isl[j].low;
		}
		else if (isl[x].par!=j)
			if (isl[x].low > isl[j].dfn)
				isl[x].low = isl[j].dfn;
	}
	

}

int getit(){
	int res = MA;
	int i, j, k;
	for (i = 1; i <= n; i++){
		j = isl[i].par;
		if (j)
			if (bc[j][i] == 1)
				if (isl[j].dfn < isl[i].low){
					
					if (gua[j][i] < res)
						res = gua[j][i];
				}
	}
	return res;
}

int main(){
	int i, j, k;
	int flag;
	while (scanf("%d%d", &n, &m)){
		if (n == 0)
			break;
		flag = 1;
		for (i = 1; i <= n; i++)
			fos.push(i);
		

		memset(isl, 0, sizeof(isl));
		memset(bc, 0, sizeof(bc));
		memset(gua, 0, sizeof(gua));
	
		cou = 1;
		for (i = 0; i < m; i++){
			scanf("%d%d%d", &u, &v, &w);
			isl[u].lik.push_back(v);
			isl[v].lik.push_back(u);
			gua[u][v] = gua[v][u] = w;
			bc[u][v]++;
			bc[v][u]++;
		}
		tar(1);
		while (!fos.empty()){
			j = fos.top();
			fos.pop();
			if (!isl[j].visited)
				flag = 0;
		}
		if (flag){
			k = getit();
			if (k == MA)
				k = -1;
			else
				k = k>0 ? k : 1;
		}
		else
			k = 0;
		printf("%d\n", k);

	}
	return 0;

}

*/

//the rabbits && the stary sky
/*
#include<iostream>
#include<algorithm>
#include<stack>
#include<memory.h>
#include<queue>
#define N 50
using namespace std;


struct edge{
	int dest;
	int weigh;
	bool operator<(const edge &e) const {
		return weigh > e.weigh;
	}
	edge &operator=(const edge& e) {
		dest = e.dest;
		weigh = e.weigh;
		return *this;
	}
	edge(){}
	edge(int i, int j) :dest(i), weigh(j){}
};

priority_queue<edge> q;

struct node{
	bool visited;
	vector<edge> edg;
}nd[N];


int mst(){
	int sum = 0;
	edge p;
	q.push(edge(0,0));

	while (!q.empty()){
		p = q.top();
		q.pop();
		if (!nd[p.dest].visited){
			nd[p.dest].visited = 1;
			sum += p.weigh;
			int siz = nd[p.dest].edg.size();
			for (int i = 0; i < siz; i++){
				if (!nd[nd[p.dest].edg[i].dest].visited)
					q.push(nd[p.dest].edg[i]);
			}
		}
	}
	return sum;
}


int main(){
	int n,m;
	char c,e;
	int i, j, k;
	cin >> n;
	memset(nd, 0, sizeof(nd));
	for (i = 0; i < n - 1; i++){
		cin >> c >> m;
		for (j = 0; j < m; j++){
			cin >> e >> k;
			nd[c - 'A'].edg.push_back(edge(e - 'A', k));
			nd[e - 'A'].edg.push_back(edge(c - 'A', k));
		}
	}
	printf("%d\n", mst());
	return 0;
}
*/

//after the earthquake(it's wrong)
/*
#include<iostream>
#include<algorithm>
#include<queue>
#include<memory.h>
#include<math.h>
#include<vector>
#include<set>
#define N 150

using namespace std;

struct way{
	int src;
	int dest;
	double len;
	bool operator<(const way &w) const{
		return len < w.len;
	}
	way &operator=(const way &w) {
		src = w.src;
		dest = w.dest;
		len = w.len;
		return *this;
	}
	way(){}
	way(int a, int b, double j) :src(a),dest(b), len(j){}
};

struct vil{
	double x;
	double y;
	bool visited;
	int indegre;
	vector<way> w;
}v[N];

int m, n;
set<way> q;
double sum;
int vc;

void mst(){
	way p;
	q.insert(way(-1, 0, 0));
	while (!q.empty()){
		set<way>::iterator it = q.begin();
		
		while (!q.empty() && v[(*it).dest].indegre)
			q.erase(it++);

		if (q.empty())
			break;

		p = *it;
		q.erase(it);

		//新加入的nk
		int nk = v[p.dest].visited ? p.src : p.dest;
		v[nk].visited = 1;
		
		v[p.dest].indegre = 1;

		sum += p.len;

		vc++;
		printf("add %d %d\n", p.src, p.dest);

		//通过nk加边
		int siz = v[nk].w.size();
		for (int i = 0; i < siz; i++){
			if (v[nk].w[i].src == nk&&!v[v[nk].w[i].dest].visited
				|| v[nk].w[i].dest == nk&&!v[v[nk].w[i].src].visited)
				q.insert(v[p.dest].w[i]);
		}


	}

}

int main(){
	int i, j, k;
	int p, q;
	double xi, yi;
	while (scanf("%d%d", &n, &m) != EOF){

		sum = 0;
		vc = 0;
		for (i = 0; i < n; i++){
			scanf("%lf%lf", &xi, &yi);
			v[i].x = xi;
			v[i].y = yi;
			v[i].w.clear();
			v[i].visited = 0;
			v[i].indegre = 0;
		}
		for (i = 0; i < m; i++){
			scanf("%d%d", &p, &q);
			double dis = sqrt(pow(v[p - 1].x - v[q - 1].x, 2) + pow(v[p - 1].y - v[q - 1].y, 2));
			v[p - 1].w.push_back(way(p - 1, q - 1, dis));
			v[q - 1].w.push_back(way(p - 1, q - 1, dis));

		}
		mst();

		if (vc < n)
			printf("NO\n");
		else
			printf("%.2lf\n", sum);

	}
	return 0;
}
*/

//表达式树(incomplete)
/*
#include<iostream>
#include<stack>
#include<algorithm>
#include<memory.h>
#include<string.h>
#define N 100

using namespace std;

char mid[N];
char rev[N];
stack<char> s;
stack<int> cal;

int n,a;
char c;


int mapping[26] = { 0 };

int main(){
	int k = 0;
	int x, y;
	int res;
	cin >> mid;
	cin >> n;
	for (int i = 0; i < n; i++){
		cin >> c >> a;
		mapping[c - 'a'] = a;
	}

	a = strlen(mid);

	
	for (int i = 0; i < a; i++){
		if (mid[i] >= 'a'&&mid[i] <= 'z')
			rev[k++] = mid[i];
		else if (mid[i] == ')'){
			while (!s.empty() && s.top() != '('){
				c = s.top();
				s.pop();
				rev[k++] = c;
			}
			s.pop();
		}
		else if (mid[i] == '(')
			s.push(mid[i]);
		else{
			while (!s.empty() && s.top() != '('){
				if (mid[i] == '*' || mid[i] == '/')
					if (s.top() == '+' || s.top() == '-')
						break;
				c = s.top();
				s.pop();
				rev[k++] = c;
			}
			s.push(mid[i]);
		}
	}

	while (!s.empty()){
		c = s.top();
		s.pop();
		rev[k++] = c;
	}

	rev[k] = 0;

	for (int i = 0; i < k; i++){
		if (rev[i] >= 'a'&&rev[i] <= 'z')
			cal.push(mapping[rev[i] - 'a']);
		else {
			x = cal.top();
			cal.pop();
			y = cal.top();
			cal.pop();
			switch (rev[i]){
			case '+':cal.push(y + x); break;
			case '-':cal.push(y - x); break;
			case '*':cal.push(y*x); break; 
			case '/':cal.push(y / x);
			}
		}
	}
	
	res = cal.top();

	printf("%s\n", rev);

	



	printf("%d\n", res);
}
*/

//currency exchange
/*
#include<iostream>
#include<algorithm>
#include<stack>
#include<memory.h>
#include<vector>
#include<queue>
#define N 200

using namespace std;


struct exc{
	int des;
	double rate;
	double com;
	exc(){}
	exc(int i,double j,double k):des(i),rate(j),com(k){}
};

struct node{
	double cou;
	vector<exc> cha;
	int upt;
	bool inq;
	node &operator=(node &n){
		cou = n.cou;
		upt = n.upt;
		int siz = n.cha.size();
		for (int i = 0; i < siz; i++)
			cha.push_back(n.cha[i]);
		return *this;
	}
}cur[N];

int n, m, s;
double v;
queue<int> q;

int spf(){
	q.push(s);
	int t, i, j, k;
	double p;
	while (!q.empty()){
		t = q.front();
		q.pop();
		cur[t].inq = 0;
		int siz = cur[t].cha.size();
		for (i = 0; i < siz; i++){
			j = cur[t].cha[i].des;
			p = (cur[t].cou - cur[t].cha[i].com)*
				cur[t].cha[i].rate;
			if (!(p<0)){
				if (p>cur[j].cou){
					if (j == s)
						return 1;
					cur[j].cou = p;
					cur[j].upt++;
					if (cur[j].upt >= n)
						return 1;
					if (!cur[j].inq){
						q.push(j);
						cur[j].inq = 1;
					}
				}
			}
		}
	}
	return 0;
}



int main(){
	int i, j, k;
	double r1, c1, r2, c2;

	memset(cur, 0, sizeof(cur));
	scanf("%d%d%d%lf", &n, &m, &s, &v);

	cur[s].cou = v;

	for (i = 0; i < m; i++){
		scanf("%d%d%lf%lf%lf%lf", &j, &k, &r1, &c1, &r2, &c2);
		cur[j].cha.push_back(exc(k, r1, c1));
		cur[k].cha.push_back(exc(j, r2, c2));

	}

	if (spf())
		printf("YES\n");
	else
		printf("NO\n");
}

*/

//single point of failure in networks
/*
#include<iostream>
#include<algorithm>
#include<queue>
#include<vector>
#include<stack>
#include<memory.h>
#include<set>

#define N 1010

using namespace std;

struct node{
	int dfn;
	int low;
	int par;
	int visited;
	int follow;
	
	vector<int> con;

}nd[N];

int ti;
int ma;
set<int> spf;

void tar(int k){
	int i, j;
	nd[k].visited = 1;
	nd[k].dfn = nd[k].low = ti++;
	int siz = nd[k].con.size();
	for (i = 0; i < siz; i++){
		j = nd[k].con[i];
		if (!nd[j].visited){
			nd[j].par = k;
			tar(j);
			if (nd[j].low < nd[k].low)
				nd[k].low = nd[j].low;

			if (nd[k].dfn <= nd[j].low)
			{
				nd[k].follow++;
				if (k!=1||(k==1&&nd[1].follow>1))
				spf.insert(k);
			}
		}
		else if (nd[k].par != j){
			if (nd[j].low < nd[k].dfn){
				nd[k].low = nd[j].dfn;
			}
		}
		
	}
	
}




int main(){
	int x, y, i, j, k;
	int p, q, t;
	int flag = 0;
	k = 1;
	while (scanf("%d", &x)){
		if (!x){

			//结束
			if (flag)
				break;

			//否则是一组数据的结束
			ti = 1;
			
			tar(1);
			//proc();

			printf("Network #%d\n", k++);

			if (spf.empty())
				printf("  No SPF nodes\n");

			else{
				set<int>::iterator it;
				for (it = spf.begin(); it != spf.end(); it++){
					int sub;
					if (*it == 1)
						sub = nd[1].follow;
					else
						sub = nd[*it].follow + 1;
					printf("  SPF node %d leaves %d subnets\n", *it, sub);
				}
			}

			printf("\n");

			ma = 0;
			flag = 1;
			spf.clear();
			memset(nd, 0, sizeof(nd));
		}


		else{
			flag = 0;
			scanf("%d", &y);
			nd[x].con.push_back(y);
			nd[y].con.push_back(x);
			if (x > ma)
				ma = x;
			if (y > ma)
				ma = y;
		}
	}

	return 0;
}
*/

//greedy verson of communication systems(wrong)
/*
#include<iostream>
#include<algorithm>
#include<queue>
#include<stack>
#include<vector>
#include<memory.h>
#include<set>
#define N 10005

using namespace std;



struct node{
	double brand;
	double price;
	int deg;
	vector<int> next;
	double per;
	bool vis;
}nd[N];

struct des{
	double minb;
	double tp;
	int place;

	des(){}
	des(int i, double x, double y) :minb(x), tp(y), place(i){}
	bool operator<(const des& d)const{
		return (minb / tp) < (d.minb / d.tp);
	}
	des& operator=(const des& s){
		minb = s.minb;
		tp = s.tp;
		place = s.place;
		return *this;
	}
};

double dij(int dep){
	des p;
	priority_queue<des> q;

	q.push(des(0, 1e8, 0));
	while (!q.empty()){
		p = q.top();
		q.pop();
		int k = p.place;

		if (nd[k].deg == dep){
			return p.minb / p.tp;
		}
		if (!nd[k].vis){
			int sz = nd[k].next.size();

			for (int i = 0; i < sz; i++){
				int j = nd[k].next[i];
				double xx = p.minb < nd[j].brand ? p.minb : nd[j].brand;
				double yy = p.tp + nd[j].price;
				double zz = xx / yy;
				if (zz>nd[j].per){
					nd[j].per = zz;
					q.push(des(j, xx, yy));
				}
			}
			nd[k].vis = 1;
		}
	}
}





int main(){
	int t, n, m;
	double x, y; 
	int count;
	int old;
	int upp;
	
	scanf("%d", &t);
	for (int i = 0; i < t; i++){
		scanf("%d", &n);
		count = 1;
		old = 0;
		upp = 1;
		memset(nd, 0, sizeof(nd));
		for (int j = 0; j < n; j++){
			scanf("%d", &m);
	
			for (int k = 0; k < m; k++){
				scanf("%lf%lf", &x, &y);
				nd[count].brand = x;
				nd[count].price = y;
				nd[count].deg = j + 1;
				
				for (int p = old; p < upp; p++)
					nd[p].next.push_back(count);
				count++;
			}
			old = upp;
			upp = count;
		}

		printf("%.3lf\n", dij(n));
		
	}

}
*/

//drainage ditches
//网络流
/*
#include<iostream>
#include<algorithm>
#include<stack>
#include<queue>
#include<memory.h>
#include<deque>
#include<vector>
#define INF 1e9
#define N 205

using namespace std;

int n, m;//m汇点
int g[N][N];
int lev[N];
int visited[N];

bool bfs(){
	memset(lev, 0, sizeof(lev));
	queue<int> q;
	q.push(1);
	
	lev[1] = 1;
	while (!q.empty()){
		int t = q.front();
		q.pop();
		if (t == m)
			return 1;

		for (int i = 1; i <= m; i++){
			if (g[t][i] > 0&&!lev[i])
			{
				lev[i] = lev[t] + 1;
				q.push(i);
			}
		}
		
	}
	return 0;
}

int dinic(){
	int res, i, j, k, t, p, q;
	deque<int> s;
	res = 0;
	while (bfs()){
		s.push_back(1);
		memset(visited, 0, sizeof(visited));
		visited[1] = 1;
		while (!s.empty()){
			j = s.size();
			t = s[j - 1];
			
			if (t == m){
				q = INF;
				p = 1;
				for (i = 1; i < j; i++){
					int b = s[i - 1];
					int e = s[i];
					if (g[b][e] > 0){
						if (g[b][e] < q){
							q = g[b][e];
							p = b;
						}
					}
				}
				res += q;
				for (i = 1; i < j; i++)
				{
					g[s[i - 1]][s[i]] -= q;
					g[s[i]][s[i - 1]] += q;
				}
				while (!s.empty()&&s.back() != p)
				{
					visited[s.back()] = 0;
					s.pop_back();
				}
			}

			else{
				
				for (i = 1; i <= m; i++){
					if (g[t][i] > 0&&!visited[i])
					{
						visited[i] = 1;
						s.push_back(i);
						break;
					}
				}

				//unable to continue
				if (i > m)
					s.pop_back();
			}
		}

	}

	return res;
}




int main(){
	int i, j, k, c;


	while (scanf("%d%d",&n,&m)!=EOF){

		memset(g, 0, sizeof(g));
		for (i = 0; i < n; i++){

			scanf("%d%d%d", &j, &k, &c);
			g[j][k] += c;
			
		}
		printf("%d\n", dinic());
	}
	
	return 0;
}
*/

//acm computer factory
//网络流
/*
#include<iostream>
#include<algorithm>
#include<stack>
#include<queue>
#include<deque>
#include<memory.h>
#define INF 1e9
#define N 105
using namespace std;

int g[N][N];
int cpy[N][N];
int visited[N];
queue<int> result;



vector<int> pt[N];

void build(int x,int p){
	int i, j, k;

	//i source,j dest
	for (i = 0; i < x; i+=2)
		for (j = 1; j <= x; j+=2){
			for (k = 0; k < p; k++)
				if (!(pt[j][k] == 2 || pt[i][k] == pt[j][k]))
					break;
			if (k == p)
			{
				g[i][j] = INF;
				cpy[i][j] = INF;
			}
		}
	
}


bool bfs(int e){
	int lev[N] = { 0 };

	queue<int> q;
	q.push(0);
	lev[0] = 1;
	while (!q.empty()){
		int t = q.front();
		if (t == e)
			return 1;
		q.pop();
		for (int i = 1; i <= e;i++)
			if (g[t][i] > 0 && !lev[i])
			{
				lev[i] = lev[t] + 1;
				q.push(i);
			}

	}

	return 0;
}

int dinic(int e){
	deque<int> q;
	
	int i, j, k, res, mi, st;
	res = 0;

	while (bfs(e)){
		memset(visited, 0, sizeof(visited));
		q.push_back(0);
		visited[0] = 1;
		while (!q.empty()){
			int siz = q.size();
			int t = q.back();
			if (t == e){
				mi = INF;
				st = 1;
				for (i = 1; i < siz; i++){
					if (g[q[i - 1]][q[i]] > 0 && g[q[i - 1]][q[i]] < mi){
						mi = g[q[i - 1]][q[i]];
						st = q[i - 1];
					}
				}

				res += mi;
				for (i = 1; i < siz; i++){
					g[q[i - 1]][q[i]] -= mi;
					g[q[i]][q[i - 1]] += mi;
				}

				while (!q.empty()&&q.back() != st){
					visited[q.back()] = 0;
					q.pop_back();
				}
			}

			else{
				for (i = 1; i <= e; i++){
					if (g[t][i] > 0 && !visited[i]){
						visited[i] = 1;
						q.push_back(i);
						break;
					}
				}

				if (i > e)
					q.pop_back();
			}

		}
	}

	return res;

}


//???????????
int proc(int e){
	int i, j;
	int c = 0;
	
	for (i = 2; i < e; i+=2)
		for (j = 1; j < e; j+=2)
		{
			if (g[i][j] < cpy[i][j]){
				c++;
				result.push(i / 2);
				result.push((j + 1) / 2);
				result.push(cpy[i][j] - g[i][j]);
			}
		}
	return c;
}



int main(){
	int p, n, q, a;
	int cou;
	int i, j, k;
	while (scanf("%d%d", &p, &n) != EOF){
		memset(pt, 0, sizeof(pt));
		memset(g, 0, sizeof(g));
		memset(cpy, 0, sizeof(cpy));
		cou = 1;
		for (i = 1; i <= n; i++){
			scanf("%d", &q);
			g[cou][cou + 1] = q;
			cpy[cou][cou + 1] = q;
			for (j = 0; j < p; j++){
				scanf("%d", &a);
				pt[cou].push_back(a);
			}
			cou++;
			for (j = 0; j < p; j++){
				scanf("%d", &a);
				pt[cou].push_back(a);
			}
			cou++;
		}
		for (j = 0; j < p; j++)
		{

			pt[0].push_back(0);
			pt[cou].push_back(1);
		}

		build(cou, p);

		printf("%d ", dinic(cou)); 
		printf("%d\n", proc(cou));
		while (!result.empty()){
			printf("%d ", result.front());
			result.pop();
			printf("%d ", result.front());
			result.pop();
			printf("%d\n", result.front());
			result.pop();
		}
	}
	return 0;
}
*/

//the rabits && flowers/I love pku(the same problem)
//using dijstra algorithm
/*
#include<iostream>
#include<algorithm>
#include<vector>
#include<stack>
#include<queue>
#include<deque>
#include<string>
#include<cstring>
#include<memory.h>
#include<set>
#define INF 1e9

using namespace std;

int g[100][100];
string s[50];
vector<int> son[50];
stack<int> pro;
int reach[50];

struct way{
	int sor;
	int des;

	//distance between the source and des,not the sor && des
	int length;
	bool operator<(const way& w)const{
		return length > w.length;
	}
	way& operator=(const way& w){
		sor = w.sor;
		des = w.des;
		length = w.length;
		return *this;
	}
	way(){}
	way(int i, int j, int k) :sor(i), des(j), length(k){}
};

void dij(int source, int dest){
	priority_queue<way> q;
	int visited[50] = { 0 };
	int ini[50];
	way t;
	t.sor = -1;
	t.des = source;
	t.length = 0;
	q.push(t);
	
	while (!q.empty()){
		t = q.top();
		q.pop();
		int k = t.des;
		int ori = t.length;
		
		visited[k] = 1;
		if (k == dest){
			break;
		}
		int siz = son[k].size();
		for (int i = 0; i < siz; i++){
			if (!visited[son[k][i]])
			{
				int su = ori + g[k][son[k][i]];
				if (su < reach[son[k][i]]){
					q.push(way(k, son[k][i], su));
					reach[son[k][i]] = su;
					ini[son[k][i]] = k;
				}
			}
		}

	}
	int fin = dest;
	while (fin != source){
		pro.push(fin);
		fin = ini[fin];
	}
	pro.push(source);

}


int main(){
	string s1, s2;
	memset(g, 0x1, sizeof(g));
	int i, j, k, n, m, t;
	scanf("%d", &n);
	for (i = 0; i < n; i++)
		cin >> s[i];
	scanf("%d", &m);
	for (i = 0; i < m; i++){
		cin >> s1 >> s2 >> t;
		for (j = 0; s[j] != s1; j++){}
		for (k = 0; s[k] != s2; k++){}
		if (t < g[j][k])g[j][k] = t;
		if (t < g[k][j])g[k][j] = t;
		son[j].push_back(k);
		son[k].push_back(j);
	}


	scanf("%d", &t);
	for (i = 0; i < t; i++){
		cin >> s1 >> s2;
		for (j = 0; s[j] != s1; j++){}
		for (k = 0; s[k] != s2; k++){}
		memset(reach, 0x1, sizeof(reach));
		dij(j, k);
		while (!pro.empty()){
			m = pro.top();
			pro.pop();
			cout << s[m];
			if (pro.empty()){
				cout << endl;
				break;
			}
			n = pro.top();
			cout << "->(" << g[m][n] << ")->";
		}

	}
	return 0;
}
*/

//不知道是什么题
/*
#include<iostream>
#include<algorithm>
#include<vector>
#include<stack>
#include<queue>
#include<deque>
#include<string>
#include<cstring>
#include<memory.h>
#define N 20000

using namespace std;

int node[N];
int emp[N];

void inser(int a, int p){
	if (!emp[p])
	{
		node[p] = a;
		emp[p] = 1;
		return;
	}
	if (a == node[p])
		return;
	if (a < node[p])
		inser(a, 2 * p + 1);
	if (a>node[p])
		inser(a, 2 * p + 2);
}

void dfs(int x){
	if (!emp[x])
		return;
	printf("%d ", node[x]);
	dfs(2 * x + 1);
	dfs(2 * x + 2);
}

int main(){
	int a;
	memset(node, 0, sizeof(node));
	memset(emp, 0, sizeof(emp));
	while (scanf("%d", &a) != EOF){
		inser(a, 0);
		
	}
	dfs(0);
	return 0;
}
*/

//optimal milking(网络流)
//实力坑杀自己
/*
#include<iostream>
#include<algorithm>
#include<queue>
#include<memory.h>
#include<stack>
#include<vector>
#include<deque>
#define N 250
#define INF 0x01010101



using namespace std;


int dist[N][N];

int net[N][N];

int lev[N];


bool bfs(int t){
	memset(lev, 0, sizeof(lev));
	
	
	queue<int> q;
	q.push(0);
	lev[0] = 1;
	while (!q.empty()){
		int p = q.front();
		q.pop();
		if (p == t)
			return 1;
		for (int i = 0; i <= t; i++){
			if (net[p][i]>0&&!lev[i])
			{
				lev[i] = lev[p] + 1;
				q.push(i);
			}
		}
	}
	return 0;
}

//t 汇点
int dinic(int t){
	deque<int> q;
	int i, j, k;
	int res = 0;
	int mi, st;
	int visited[N];

	while (bfs(t)){

		memset(visited, 0, sizeof(visited));
		q.push_back(0);
		visited[0] = 1;
		while (!q.empty()){
			
			int siz = q.size();
			int nd = q.back();

			if (nd == t){
				mi = INF;
				st = 0;
				for (i = 1; i < siz; i++){
					if (net[q[i - 1]][q[i]]>0 &&
						net[q[i - 1]][q[i]] < mi)
					{
						mi = net[q[i - 1]][q[i]];
						st = q[i - 1];
					}
				}

				res += mi;
				for (i = 1; i < siz; i++){
					net[q[i - 1]][q[i]] -= mi;
					net[q[i]][q[i - 1]] += mi;
				}

				while (!q.empty() && q.back() != st){
					visited[q.back()] = 0;
					q.pop_back();
				}

				
			}

			else {
				for (i = 0; i <= t;i++)
					if (net[nd][i]>0 && (lev[i] == lev[nd] + 1) && !visited[i]){
						q.push_back(i);
						visited[i] = 1;
						break;
					}
				if (i > t){
					q.pop_back();
				}

			}
		}

	}
	//printf("%d\n", res);
	return res;
}

int enu(int cow, int mer, int con, int limit){
	int f, m, e;
	int i, j, k;
	f = 1;
	e = limit;
	while (f != e){
		
		m = (f + e) / 2;
		
		//build the net

		//s to cow
		memset(net, 0, sizeof(net));
		for (i = 1; i <= cow; i++)
			net[0][i] = 1;

		//cow to mechine
		for (i = 0; i < cow; i++)
			for (j = 0; j < mer; j++)
				if (dist[i + mer][j] <= m)
					net[i + 1][cow + j + 1] = 1;

		//mechine to t
		for (j = 1; j <= mer; j++)
			net[cow+j][cow+mer+1] = con;

		
		//get the max flow
		if (dinic(cow + mer + 1) >= cow){
			e = m;
		}
		else{
			f = m + 1;
		}
	}
	return e;
}

int main(){
	//k 机器数
	int a, i, j, k, p, q, c, m;
	memset(dist, 0x1, sizeof(dist));
	scanf("%d%d%d", &k, &c, &m);

	for (i = 0; i < k + c; i++)
		for (j = 0; j < k + c; j++)
		{
			scanf("%d", &a);
			if (a || (i == j))
				dist[i][j] = a;
		}


	//floyed
	for (p = 0; p < k + c; p++)
		for (i = 0; i < k + c; i++)
			for (j = 0; j < k + c; j++)
				if (dist[i][j]>dist[i][p] + dist[p][j])
					dist[i][j] = dist[i][p] + dist[p][j];

	int mat = 0;
	for (i = 0; i < c; i++)
		for (j = 0; j < k; j++)
			if (dist[i + k][j]<INF&&dist[i + k][j]>mat)
				mat = dist[i + k][j];

	printf("%d\n", enu(c, k, m, mat));

	return 0;

}
*/

//the perfect stall
/*
#include<iostream>
#include<algorithm>
#include<deque>
#include<queue>
#include<vector>
#include<memory.h>
#define N 500
#define INF 0x7fffffff

using namespace std;

int net[N][N];
int lev[N];

bool bfs(int sum){
	queue<int> q;
	q.push(0);
	memset(lev, 0, sizeof(lev));
	lev[0] = 1;

	while (!q.empty()){
		int t = q.front();
		q.pop();

		if (t == sum)
			return 1;

		for (int i = 0; i <= sum; i++){
			if (net[t][i] > 0 && !lev[i]){
				lev[i] = lev[t] + 1;
				q.push(i);
			}
		}
	}

	return 0;
}


int dinic(int sum){
	deque<int> q;
	int i, j, k, res, mi, st;
	int visited[N];

	res = 0;

	while (bfs(sum)){
		memset(visited, 0, sizeof(visited));
		q.push_back(0);
		visited[0] = 1;

		while (!q.empty()){
			int t = q.back();
			int siz = q.size();

			if (t == sum){
				mi = INF;
				st = 0;
				for (i = 1; i < siz; i++){
					if (net[q[i - 1]][q[i]] > 0 && net[q[i - 1]][q[i]] < mi){
						mi = net[q[i - 1]][q[i]];
						st = q[i - 1];
					}
				}

				res += mi;

				for (i = 1; i < siz; i++){
					net[q[i - 1]][q[i]] -= mi;
					net[q[i]][q[i - 1]] += mi;
				}

				while (!q.empty() && q.back() != st){
					visited[q.back()] = 0;
					q.pop_back();
				}
			}

			else{
				for (i = 0; i <= sum; i++){
					if (net[t][i] > 0 && !visited[i] && lev[i] == lev[t] + 1){
						q.push_back(i);
						visited[i] = 1;
						break;
					}
				}

				if (i > sum)
					q.pop_back();
			}

		}

	}

	return res;
}

int main(){
	int n, m, s, i, j, k;
	while (scanf("%d %d", &n, &m) != EOF){
		
		memset(net, 0, sizeof(net));

		for (i = 1; i <= n; i++){
			scanf("%d", &s);
			for (j = 0; j < s; j++){
				scanf("%d", &k);
				net[i][n + k] = 1;
			}
		}

		for (i = 1; i <= n; i++)
			net[0][i] = 1;
		for (i = n + 1; i <= n + m; i++)
			net[i][n + m + 1] = 1;

		printf("%d\n", dinic(n + m + 1));

	}


	return 0;
}


*/

//线段树求逆序数
/*
#include<iostream>
#include<algorithm>
#define N 25000

using namespace std;


struct node{
	int s;
	int e;
	int cou;

}nd[4*N];

int a[N];

//离散化使用
int cpy[N];

void build(int k,int x,int y){
	nd[k].cou = 0;
	nd[k].s = x;
	nd[k].e = y;

	if (x == y)
		return;

	int m = (x + y) / 2;
	build(2 * k + 1, x, m);
	build(2 * k + 2, m + 1, y);
}

int inser(int k, int p){
	nd[k].cou++;

	//叶节点或只存一个值的非叶节点
	if (cpy[nd[k].s] == cpy[nd[k].e]){
		if (p > cpy[nd[k].s])
			return nd[k].cou;
		return 0;
	}
		
	int s = 0;

	//位于右子树
	if (p > cpy[nd[2 * k + 1].e])
	{
		s += nd[2 * k + 1].cou;
		s += inser(2 * k + 2, p);
	}

	//
	else
		s += inser(2 * k + 1, p);
	return s;
}


int main(){
	int n, i, j, k;
	int sum;
	while (scanf("%d", &n)){
		if (!n)
			return 0;
		for (i = 0; i < n; i++){
			scanf("%d", a + i);
			cpy[i] = a[i];
		}

		sort(cpy, cpy + n);
		build(0, 0, n - 1);
		sum = 0;

		for (i = n - 1; i >= 0; i--)
			sum += inser(0, a[i]);

		printf("%d\n", sum);

	}
	return 0;
}
*/

//sequence
//问题可分步，以已得最小和为基础，每次加一组
/*
#include<iostream>
#include<algorithm>
#include<memory.h>
#include<queue>
#define MAX 0x7fffffff

using namespace std;

int a[2005], b[2005];
int p[105], np[105];
priority_queue<int> q;
//priority_queue<int, vector<int>, less<int>> cq;

int main(){
	int t, m, n;
	int i, j, k;
	scanf("%d", &t);
	while (t--){
		memset(a, 0, sizeof(a));

		scanf("%d%d", &m, &n);
		for (i = 0; i < n; i++)
			scanf("%d", a + i);

		sort(a, a + n);

		for (i = 1; i < m; i++){

			for (j = 0; j < n; j++){
				scanf("%d", b+j);
			}
			sort(b, b + n);

			for (j = 0; j < n; j++)
				q.push(a[0] + b[j]);

			for (j = 1; j < n; j++)
				for (k = 0; k < n; k++){
					if (a[j] + b[k] < q.top()){
						q.pop();
						q.push(a[j] + b[k]);
					}
				}

			for (j = n - 1; j >= 0&&!q.empty(); j--){
				a[j] = q.top();
				q.pop();
			}
		}

		for (i = 0; i < n; i++){
			printf("%d", a[i]);
			if (i == n - 1)
				printf("\n");
			else
				printf(" ");
		}

	}
	return 0;
}
*/

//树状数组求逆序数
/*
#include<iostream>
#include<algorithm>
#include<memory.h>

#define N 25000
using namespace std;

int lowbit(int x){
	return x&(-x);
}

struct pos{
	int a;
	int ori;
	bool operator<(pos &p){
		if (a == p.a)
			return ori < p.ori;
		return a < p.a;
	}
}p[N];

struct node{
	int rank;
	int old;
	bool operator<(node &n){
		return old < n.old;
	}
}nd[N];


int c[N];

int sum(int k){
	int i = k;
	int res = 0;
	while (i != 0){
		res += c[i];
		i -= lowbit(i);
	}
	return res;
}

void inser(int r,int limit){
	for (int i = r; i <= limit; i += lowbit(i))
		c[i]++;

	return;
}

int main(){
	int n, i, j, k;
	int s;
	while (scanf("%d", &n)){
		if (!n)return 0;
		s = 0;
		for (i = 0; i < n; i++)
		{
			scanf("%d", &j);
			p[i].a = j;
			p[i].ori = i;
		}
		sort(p, p + n);
		for (i = 0; i < n; i++){
			nd[i].rank = i;
			nd[i].old = p[i].ori;
		}

		sort(nd, nd + n);

		memset(c, 0, sizeof(c));

		for (i = n-1; i >= 0; i--){
			s += sum(nd[i].rank);
			inser(nd[i].rank+1, n);
		}

		printf("%d\n", s);
	}
}
*/

//后缀数组（可当作模板）
//long long message
/*
#include<iostream>
#include<algorithm>
#include<cstring>

#define N 200005


using namespace std;

//char string
char LL[N];
//后缀
int s[N];
//rank
int r[N];

int height[N];

//build
int wa[N], wb[N], wv[N], Ws[N]; 
//辅助数组 int sa[MAXN]; 
//sa[i]是名次为i的后缀的位置 
void BuildSA(const char * s, int sa[], int n, int m) {
	int i, j, p, *pm = wa, *k2sa = wb, *t;
	for (i = 0; i<m; i++) Ws[i] = 0;
	for (i = 0; i<n; i++) Ws[pm[i] = s[i]]++;
	for (i = 1; i<m; i++) Ws[i] += Ws[i - 1];
	for (i = n - 1; i >= 0; i--) sa[--Ws[pm[i]]] = i;
	for (j = p = 1; p < n; j <<= 1, m = p) {
		for (p = 0, i = n - j; i < n; i++) k2sa[p++] = i;
		for (i = 0; i < n; i++)
			if (sa[i] >= j)
				k2sa[p++] = sa[i] - j;
		for (i = 0; i < m; i++)
			Ws[i] = 0;
		for (i = 0; i < n; i++)
			Ws[wv[i] = pm[k2sa[i]]]++;
		for (i = 1; i < m; i++)
			Ws[i] += Ws[i - 1];
		for (i = n - 1; i >= 0; i--)
			sa[--Ws[wv[i]]] = k2sa[i];
		for (t = pm, pm = k2sa, k2sa = t, pm[sa[0]] = 0, p = i = 1; i<n; i++){
			int a = sa[i - 1], b = sa[i];
			if (k2sa[a] == k2sa[b] && k2sa[a + j] == k2sa[b + j])
				pm[sa[i]] = p - 1;
			else pm[sa[i]] = p++; 
		}
	}
	return;
}

void BuildHeight(char * str, int n, int * sa, int * Rank) {
	int i, j, k; for (int i = 0; i < n + 1; ++i) //i 是名次,n是字符串长度
		Rank[sa[i]] =  i; 
	for (i = k = 0; i < n; height[Rank[i++]] = k)//i是位置
		for (k ? k-- : 0, j = sa[Rank[i] - 1]; //Rank[0]>0才不越界
			str[i+k]==str[j+k]; k++); //k相当于是 H[i]; height[Rank[i]] =  H[i] ;
}



int main(){
	int i, j, k;
	char c;
	int pos = 0;
	int cnt = 0;
	int mid;

	//get it
	while (cnt != 2){
		scanf("%c", &c);
		if (c == '\n'){
			cnt++;
			if (cnt == 1){
				mid = pos;
				LL[pos++] = '-';
			}
		}
		else
			LL[pos++] = c;
	}
	LL[pos] = '\0';

	BuildSA(LL, s, pos+1, 127);

	BuildHeight(LL, pos, s, r);

	int mi = 0;

	for (i = 0; i < pos; i++){
		if (height[i]>mi)
			if (s[i - 1]<mid&&s[i]>mid || s[i - 1] > mid&&s[i] < mid)
				mi = height[i];
	}

	printf("%d\n", mi);
	return 0;
}
*/

//corporate identity
//后缀数组
//字符串视为数组，拼合，用不在串中的数分隔
//二分找可能的最大公共字串
/*
#include<iostream>
#include<cstring>
#include<algorithm>
#include<memory.h>
#define N 4005
#define L 205

using namespace std;

int str[N*L];
int s[N*L];
int ran[N*L];
int height[N*L];
int symbol[N];
int mrk[N*L];

int *theptr;
int commonlen;

//build
int wa[N*L], wb[N*L], wv[N*L], Ws[N*L];
//辅助数组 int sa[MAXN]; 
//sa[i]是名次为i的后缀的位置 
void BuildSA(const int * s, int sa[], int n, int m) {
	int i, j, p, *pm = wa, *k2sa = wb, *t;
	for (i = 0; i<m; i++) Ws[i] = 0;
	for (i = 0; i<n; i++) Ws[pm[i] = s[i]]++;
	for (i = 1; i<m; i++) Ws[i] += Ws[i - 1];
	for (i = n - 1; i >= 0; i--) sa[--Ws[pm[i]]] = i;
	for (j = p = 1; p < n; j <<= 1, m = p) {
		for (p = 0, i = n - j; i < n; i++) k2sa[p++] = i;
		for (i = 0; i < n; i++)
			if (sa[i] >= j)
				k2sa[p++] = sa[i] - j;
		for (i = 0; i < m; i++)
			Ws[i] = 0;
		for (i = 0; i < n; i++)
			Ws[wv[i] = pm[k2sa[i]]]++;
		for (i = 1; i < m; i++)
			Ws[i] += Ws[i - 1];
		for (i = n - 1; i >= 0; i--)
			sa[--Ws[wv[i]]] = k2sa[i];
		for (t = pm, pm = k2sa, k2sa = t, pm[sa[0]] = 0, p = i = 1; i<n; i++){
			int a = sa[i - 1], b = sa[i];
			if (k2sa[a] == k2sa[b] && k2sa[a + j] == k2sa[b + j])
				pm[sa[i]] = p - 1;
			else pm[sa[i]] = p++;
		}
	}
	return;
}

void BuildHeight(int * str, int n, int * sa, int * Rank) {
	int i, j, k; for (int i = 0; i < n + 1; ++i) //i 是名次,n是字符串长度
		Rank[sa[i]] = i;
	for (i = k = 0; i < n; height[Rank[i++]] = k)//i是位置
		for (k ? k-- : 0, j = sa[Rank[i] - 1]; //Rank[0]>0才不越界
			str[i + k] == str[j + k]; k++); //k相当于是 H[i]; height[Rank[i]] =  H[i] ;
}

int* check(int test,int totlen,int strcou){
	int sum;
	int beg;
	int *p = NULL;
	for (int i = 0; i <= totlen - strcou; i++){
		memset(symbol, 0, sizeof(symbol));
		sum = 0;
		beg = i;
		while (i < totlen&&height[i + 1] >= test){
			if (!symbol[mrk[s[i]]]){
				symbol[mrk[s[i]]] = 1;
				sum++;
			}
			i++;
		}

		//last one
		if (!symbol[mrk[s[i]]]){
			symbol[mrk[s[i]]] = 1;
			sum++;
		}

		if (sum == strcou){
			p = str + s[beg];
			return p;
		}
	}

	return NULL;
}

bool binser(int totlen,int strcou){
	int i, j, k;
	int left = 1;
	int right = 200;
	int mid;
	int *p = NULL;
	while (left != right){
		mid = (left + right) / 2 + 1;
		if ((p = check(mid, totlen, strcou)) != NULL){
			left = mid;
		}
		else
			right = mid - 1;
	}
	theptr = check(left, totlen, strcou);
	commonlen = left;
	return theptr != NULL;

}

int main(){
	int n, i, j, k;
	char c;
	int sc;
	int oth;
	int pos;
	while (scanf("%d", &n)){
		getchar();
		if (!n)break;
		//string count 
		sc = 0;
		//more chars
		oth = 130;
		//length
		pos = 0;
		//mark
		memset(mrk, 0, sizeof(mrk));
		while (scanf("%c", &c)){
			mrk[pos] = sc;
			if (c == '\n'){
				sc++;
				if (sc == n)break;
				str[pos++] = oth++;
			}
			else
				str[pos++] = c;
		}
		str[pos] = 0;

		BuildSA(str, s, pos + 1, oth);
		BuildHeight(str, pos, s, ran);

		if (binser(pos, n)){
			for (i = 0; i < commonlen; i++){
				printf("%c", *(theptr + i));
			}
			printf("\n");
		}

		else
			printf("IDENTITY LOST\n");
	}

	return 0;
}
*/

//take care: the result might be a long long int!!!!
/*
#include<iostream>
#include<algorithm>
#include<memory.h>

#define N 100005
using namespace std;

int lowbit(int x){
	return x&(-x);
}

struct pos{
	int a;
	int ori;
	bool operator<(pos &p){
		if (a == p.a)
			return ori > p.ori;
		return a < p.a;
	}
}p[N];

struct node{
	int rank;
	int old;
	bool operator<(node &n){
		return old < n.old;
	}
}nd[N];


int c[N];

int sum(int k){
	int i = k;
	int res = 0;
	while (i != 0){
		res += c[i];
		i -= lowbit(i);
	}
	return res;
}

void inser(int r, int limit){
	for (int i = r; i <= limit; i += lowbit(i))
		c[i]++;

	return;
}

int main(){
	int n, i, k;
	int j;
	long long s;
	scanf("%d", &n);
	s = 0;
	for (i = 0; i < n; i++)
	{
		scanf("%d", &j);
		p[i].a = j;
		p[i].ori = i;
	}
	sort(p, p + n);
	for (i = 0; i < n; i++){
		nd[i].rank = i;
		nd[i].old = p[i].ori;
	}

	sort(nd, nd + n);

	memset(c, 0, sizeof(c));

	for (i = 0; i < n; i++){
		s += sum(nd[i].rank);
		inser(nd[i].rank + 1, n);
	}

	printf("%lld\n", s);

}
*/

//???
/*
#include<iostream>
#include<algorithm>
#include<memory.h>
#include<vector>
#define N 1000005
#define MA 0x7fffffff
#define MI 0x80000000
using namespace std;

struct node{
	int bg;
	int ed;
	int mi;
	int ma;
}nd[4*N];

int tmin, tmax;

vector<int> minmun;
vector<int> maxmun;

void build(int b, int e, int k){
	nd[k].bg = b;
	nd[k].ed = e;
	nd[k].ma = MI;
	nd[k].mi = MA;
	if (b == e)return;
	int mid = (b + e) / 2;
	build(b, mid, 2 * k + 1);
	build(mid + 1, e, 2 * k + 2);
}

void inser(int k, int p, int num){
	if (num > nd[k].ma)
		nd[k].ma = num;
	if (num < nd[k].mi)
		nd[k].mi = num;

	if (nd[k].bg == nd[k].ed)
		return;
	if (p <= nd[2 * k + 1].ed)
		inser(2 * k + 1, p, num);
	else
		inser(2 * k + 2, p, num);
}

void query(int b, int e,int k){

	if (b <= nd[k].bg&&e >= nd[k].ed){
		if (nd[k].ma>tmax)
			tmax = nd[k].ma;
		if (nd[k].mi < tmin)
			tmin = nd[k].mi;
		return;
	}
	if (b <= nd[2 * k + 1].ed)
		query(b, e, 2 * k + 1);
	if (e >= nd[2 * k + 2].bg)
		query(b, e, 2 * k + 2);
}

int main(){
	int n, i, j, k, a;
	scanf("%d%d", &n, &k);
	build(0, n - 1, 0);
	for (i = 0; i < n; i++){
		scanf("%d", &a);
		inser(0, i, a);
	}

	for (i = 0; i < n - k + 1; i++){
		tmax = MI;
		tmin = MA;
		query(i, i + k - 1, 0);
		maxmun.push_back(tmax);
		minmun.push_back(tmin);
	}

	for (i = 0; i < n - k + 1; i++){
		printf("%d", minmun[i]);
		if (i == n - k)
			printf("\n");
		else
			printf(" ");
	}
	for (i = 0; i < n - k + 1; i++){
		printf("%d", maxmun[i]);
		if (i == n - k)
			printf("\n");
		else
			printf(" ");
	}
	return 0;
}
*/

//???
/*
#include<iostream>
#include<algorithm>
#include<set>

using namespace std;


struct guest{
	int priority;
	int num;
	bool operator<(const guest &g)const{
		return priority < g.priority;
	}
}g;

multiset<guest> s;
int main(){
	int n, i, j, k, p;
	while (scanf("%d", &j)){
		if (!j)
			return 0;
		if (j == 1){
			scanf("%d%d", &k, &p);
			g.num = k;
			g.priority = p;
			s.insert(g);
		}
		else {
			if (s.empty())
				printf("0\n");
			else{
				set<guest>::iterator it;
				if (j == 2){
					it = s.end();
					it--;
				}
				else
					it = s.begin();
				printf("%d\n", it->num);
				s.erase(it);
			}
		}
		
	}

	return 0;
}
*/

//???
/*
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;

const int MAXN = 100010;
const int MAXM = 200010;
struct Trie {
	Trie *children[2];
	int value;
} root, pool[MAXN << 5], *ptr, emptyNode;
// Trie树是一棵二叉树，左孩子表示当前二进制位为0，右孩子表示当前二进制位为1
// 叶子节点额外保存value，表示当前路径对应的数，如根-左-右-右这个节点的value为(011)2=3

int n;
// 节点个数

int head[MAXN], succeed[MAXM], vertex[MAXM], weight[MAXM], now;
// head,succeed使用链表保存边表；vertex,weight存每条边的终点以及权值

int rootXor[MAXN];
// 每个点到根的路径的xor值

void addEdge(int u, int v, int w) {
	// 往链表中添加一条u->v，权值为w的边
	succeed[now] = head[u];
		// 新的节点的后继为u节点原来的链表头
		vertex[now] = v;
	// 新的边的终点为v
	weight[now] = w;
	// 新的边的权值为w
	head[u] = now++;
	// 指定u节点新的链表头为now
}

void readTree() {
	// 读入树
	memset(head, -1, sizeof head);
	now = 0;
	for (int i = 0; i < n - 1; i++) {
		int u, v, w;
		scanf("%d%d%d", &u, &v, &w);
		addEdge(u, v, w);
		addEdge(v, u, w);
	}
}

void DFS(int x = 0, int father = -1) {
	// DFS求出每个点到0号点的路径的xor值，并填到rootXor数组中
	// 参数x表示当前递归到的节点，father表示当前点的父亲
	for (int now = head[x]; now != -1; now = succeed[now]) {
		// 遍历x节点的链表，找到所有x的出边
		int y = vertex[now];
		int w = weight[now];
		// 从x到y有一条权值为w的边
		if (y==father)
			continue;
		// 忽略返回父亲的边
		rootXor[y] = w^rootXor[x];
		// 计算rootXor[y]的值，为y到0号点的路径的xor值
		DFS(y, x);
	}
}

int getNthBit(int value, int nBit) {
	return value >> nBit & 1;
}

void insertTrie(Trie *node, int value, int nBit) {
	// 往Trie中插入value
	// 当前在node这个节点，处理到value的第nBit个bit
	// 如果nBit为-1，说明已经处理完毕，node为叶子节点
	if (nBit == -1)
		node->value = value;
	else {
		bool bit = getNthBit(value,nBit);
		// bit为0说明当前的节点应为node的左孩子
		// bit为1说明当前的节点应为node的右孩子
		if (!node->children[bit]) {
			node->children[bit] = ptr++;
			*node->children[bit] = emptyNode;
		}
		insertTrie(node->children[bit], value, nBit - 1);
	}
}

void buildTrie() {
	// 将所有rootXor[i]插入Trie中
	root = emptyNode;
	ptr = pool;
	for (int i = 0; i < n; i++)
		insertTrie(&root, rootXor[i], 30);
}

int queryTrie(Trie *node, int value, int nBit) {
	// 在Trie中查询与value的xor值最大的答案
	// 当前在node这个节点，处理到value的第nBit个bit
	// 如果nBit为-1，说明已经处理完毕，node为叶子节点
	if (nBit == -1)
		return node->value ^ value;
	else {
		bool bit = getNthBit(value, nBit);
		// 同上
		if (node->children[!bit])
			return queryTrie(node->children[!bit], value, nBit - 1);
		else
			return queryTrie(node->children[bit], value, nBit - 1);
	}
}

int getAns() {
	// 将所有rootXor[i]在Trie中查询最大的xor答案
	int ans = 0;
	for (int i = 0; i < n; i++)
		ans = max(ans, queryTrie(&root, rootXor[i], 30));
	return ans;
}

int main() {
	while (scanf("%d", &n) == 1) {
		readTree();
		DFS();
		buildTrie();
		printf("%d\n", getAns());
	}
	return 0;
}
*/

//kmp
/*
#include<iostream>
#include<algorithm>
#include<stack>
#include<cstring>
#include<vector>
#include<memory.h>
#define N 400005
#define MIN(a,b) (a<b?a:b)

using namespace std;
stack<int> s;
char str[N];
int nxt[N];

void getfix(int m){
	int i, k;
	memset(nxt, 0, sizeof(nxt));
	nxt[0] = -1;
	i = 0; k = -1;
	while (i < m){
		while (k >= 0 && str[k] != str[i]){
			k = nxt[k];
		}
		i++; k++;
		nxt[i] = k;
	}

}

int main(){
	int len;
	int j;
	int mrk;
	while (scanf("%s", str) != EOF){

		len = strlen(str);

		getfix(len);

		s.push(len);
		for (j = nxt[len]; j > 0; j = nxt[j]){
			s.push(j);
		}
		mrk = 0;

		while (!s.empty()){
			if (!mrk)
			{
				printf("%d", s.top());
				mrk = 1;
			}
			else
				printf(" %d", s.top());
			s.pop();
		}
		printf("\n");
	}

}
*/

//square
//enum with brute force 
/*
#include<iostream>
#include<algorithm>
#include<set>
#define N 1010

using namespace std;


struct point{
	int x, y;
	bool operator<(const point &p)const
	{
		if (x == p.x)
			return y < p.y;
		return x < p.x;
	}
	bool operator==(const point &p)const{
		return (x == p.x) && (y == p.y);
	}
}tp,p[N];


int n;

int binser(int x,int y){
	point tp;
	int l, r, m;
	tp.x = x;
	tp.y = y;
	l = 0; r = n - 1;
	while (l <= r){
		m = (l + r) / 2;
		if (tp == p[m])
			return m;
		if (tp < p[m])
			r = m - 1;
		else
			l = m + 1;
	}
	
	return -1;

}

int ser(int a, int b){
	int x1, y1, x2, y2, x3, y3, x4, y4;
	x1 = p[a].x + p[a].y - p[b].y;
	y1 = p[a].y + p[b].x - p[a].x;
	x2 = p[b].x + p[a].y - p[b].y;
	y2 = p[b].y + p[b].x - p[a].x;
	x3 = p[a].x + p[b].y - p[a].y;
	y3 = p[a].y + p[a].x - p[b].x;
	x4 = p[b].x + p[b].y - p[a].y;
	y4 = p[b].y + p[a].x - p[b].x;

	int p1 = -1, p2 = -1, p3 = -1, p4 = -1;
	int res = 0;
	p1 = binser(x1, y1);
	if (p1 >= 0)
		p2 = binser(x2, y2);
	if (p2 >= 0)
		res++;
	p3 = binser(x3, y3);
	if (p3 >= 0)
		p4 = binser(x4, y4);
	if (p4 >= 0)
		res++;

	return res;

}


int main(){
	int i, j,  x, y, c;
	long long k;

	while (scanf("%d", &n)){
		if (!n)break;
		for (i = 0; i < n; i++){
			scanf("%d%d", &x, &y);
			p[i].x = x, p[i].y = y;
		}
		k = 0;
		sort(p, p + n);
	

		for (i = 0; i < n - 1; i++){
			for (j = i + 1; j < n; j++){
				k += ser(i, j);
			}
		}

		printf("%lld\n", k/4);
	}

}
*/

//intersecting lines
//计算几何
//运用类/结构体 谨防编译错误
/*
#include<iostream>
#include<algorithm>
#include<stack>
#include<vector>
#include<memory.h>
#include<string>
#include<math.h>
#define INF 0x7fffffff


using namespace std;
double EPS = 1e-6;

struct point{
	double x, y;
	point(){}
	point(double i,double j){
		x = i; y = j;
	}
	point operator+(point p){
		return point(x + p.x, y + p.y);
	}
	point operator-(point p){
		return point(x - p.x, y - p.y);
	}
	point operator*(double a){
		return point(x*a, y*a);
	}
	point &operator=(point p){
		x = p.x; y = p.y;
		return *this;
	}
};

double operator*(point p1, point p2){
	return p1.x*p2.x + p1.y*p2.y;
}

double operator^(point p1,point p2){
	return p1.x*p2.y - p2.x*p1.y;
}

double area(point p, point q){
	return (p^q) / 2;
}

struct line{
	point a, b;
	line(){}
	line(point c, point d){
		a = c;
		b = d;
	}
};

bool isz(double t){
	return t > -EPS&&t < EPS;
}

double len(point l){
	return sqrt(l*l);
}

double dist(point p, line l){
	return fabs((p - l.a) ^ (l.b - l.a)) / len(l.b - l.a);
}

point intersect(line l,line m, string &g){
	double x = area(m.a - l.a, l.b - l.a);
	double y = area(l.b - l.a, m.b - l.a);
	
	if (isz(x + y)){
		if (isz(dist(l.a, m)))
			g = "LINE";
		else
			g = "NONE";
		return point(INF,INF);
	}

	g = "POINT";
	return m.a + (m.b - m.a)*(x / (x + y));

}

int main(){
	int n, i, j, k;
	string msg;
	point res = point(INF, INF);
	double st[8];
	scanf("%d", &n);
	printf("INTERSECTING LINES OUTPUT\n");
	for (i = 0; i < n; i++){
		for (j = 0; j < 8; j++)
			scanf("%lf", st + j);
		res = intersect(line(point(st[0],st[1]), point(st[2],st[3])),
			line(point(st[4],st[5]), point(st[6],st[7])), msg);
		
		printf("%s", msg.c_str());

		if (msg=="POINT")
			printf(" %.2lf %.2lf", res.x, res.y);
		printf("\n");
	}
	printf("END OF OUTPUT\n");
}
*/

//Myacm Triangles
/*
#include<iostream>
#include<algorithm>
#include<math.h>
#include<memory.h>
#define INF 0x7fffffff

using namespace std;
double EPS = 1e-6;

double maxs = 0;
char tri[4];

struct point{
	char syb;
	double x, y;
	point(){}
	point(double i, double j){
		x = i; y = j;
	}
	point operator+(point p){
		return point(x + p.x, y + p.y);
	}
	point operator-(point p){
		return point(x - p.x, y - p.y);
	}
	point operator*(double a){
		return point(x*a, y*a);
	}
	point &operator=(point p){
		x = p.x; y = p.y;
		return *this;
	}
}p[20];

double operator*(point p1, point p2){
	return p1.x*p2.x + p1.y*p2.y;
}

double operator^(point p1, point p2){
	return p1.x*p2.y - p2.x*p1.y;
}

//square (between two vectors)
double area(point p, point q){
	return fabs((p^q) / 2);
}

struct line{
	point a, b;
	line(){}
	line(point c, point d){
		a = c;
		b = d;
	}
};

bool isz(double t){
	return t > -EPS&&t < EPS;
}

double len(point l){
	return sqrt(l*l);
}

bool order(int a,int b,int c){
	return ((p[b] - p[a]) ^ (p[c] - p[b])) > EPS;
}

double dist(point p, line l){
	return fabs((p - l.a) ^ (l.b - l.a)) / len(l.b - l.a);
}

//if a point were on the leftside of a line or on the line 
bool leftoron(int a,int b,int k){
	return !(((p[b] - p[a]) ^ (p[k] - p[a])) < -EPS);
}

double init(int a,int b,int c,int k){
	return leftoron(a, b, k) && leftoron(b, c, k)
		&& leftoron(c, a, k);
}

void test(int a,int b,int c,int lim){
	int i, j, k, m;
	for (m = 0; m < lim; m++){
		if (!(m == a || m == b || m == c)){
			if (init(a, b, c, m))
				return;
		}
	}
	double s = area(p[a] - p[b], p[c] - p[b]);
	if (s > maxs){
		maxs = s;
		tri[0] = p[a].syb;
		tri[1] = p[b].syb;
		tri[2] = p[c].syb;
		tri[3] = '\0';
	}
}


int main(){
	int i, j, k, n;
	double x, y;
	char c;
	while (scanf("%d", &n)){
		if (!n)
			break;
		for (i = 0; i < n; i++){
			getchar();
			scanf("%c%lf%lf", &c, &x, &y);
			p[i].syb = c; p[i].x = x; p[i].y = y;
		}
		maxs = 0;
		for (i = 0; i < n - 2; i++)
			for (j = i + 1; j < n - 1; j++)
				for (k = j + 1; k < n; k++){
					if (order(i, j, k))
						test(i, j, k, n);
					else
						test(i, k, j, n);
				}
		sort(tri, tri + 3);
		printf("%s\n", tri);
	}
}

*/

//field reduction
//有机会再交一次（哇，大哭）
//这个是wa
//可以写个树状数组版本的（那个好写些啊啊）
/*
#include<iostream>
#include<algorithm>
#include<stack>
#include<vector>
#include<string>
#include<set>
#include<map>
#include<queue>
#include<memory.h>

#define N 50005
#define MIN 0x80000000
#define MAX 0x7fffffff

#define big(x,y) (x>y?x:y)
#define smal(x,y) (x<y?x:y)

using namespace std;


struct point {
	int x, y;
	bool operator<(const point &p)const {
		if (x == p.x)
			return y < p.y;
		return x < p.x;
	}
}p[N];


int range;
int n, k;
int ma = MIN;
int mi = MAX;
int mrk[N];

struct node {
	int b, e;

	int xma;
	int yma;
	int xmi;
	int ymi;

}nd[4 * N];

void build(int b, int e, int k) {
	nd[k].b = b;
	nd[k].e = e;

	nd[k].xma = nd[k].yma = MIN;
	nd[k].xmi = nd[k].ymi = MAX;

	if (b == e)
		return;
	int m = (b + e) / 2;
	build(b, m, 2 * k + 1);
	build(m + 1, e, 2 * k + 2);
}


void dele(int i, int k) {
	if (nd[k].b == nd[k].e) {
		nd[k].xma = nd[k].yma = MIN;
		nd[k].xmi = nd[k].ymi = MAX;
		return;
	}

	if (i < nd[2 * k + 2].b)
		dele(i, 2 * k + 1);
	else
	{
		dele(i, 2 * k + 2);
	}

	nd[k].xma = big(nd[2 * k + 1].xma, nd[2 * k + 2].xma);
	nd[k].yma = big(nd[2 * k + 1].yma, nd[2 * k + 2].yma);
	nd[k].xmi = smal(nd[2 * k + 1].xmi, nd[2 * k + 2].xmi);
	nd[k].ymi = smal(nd[2 * k + 1].ymi, nd[2 * k + 2].ymi);


}

void inser(int i, int k) {
	if (nd[k].xma < p[i].x)
		nd[k].xma = p[i].x;
	if (nd[k].xmi > p[i].x)
		nd[k].xmi = p[i].x;
	if (nd[k].ymi > p[i].y)
		nd[k].ymi = p[i].y;
	if (nd[k].yma < p[i].y)
		nd[k].yma = p[i].y;

	if (nd[k].b == nd[k].e)
		return;

	if (i < nd[2 * k + 2].b)
		inser(i, 2 * k + 1);
	else
	{
		inser(i, 2 * k + 2);
	}

}


int main() {
	int x, y;
	stack<int> rm;
	int sq;
	int mk;
	scanf("%d", &n);
	if (n <= 3) {
		printf("0\n");
		return 0;
	}
	for (int i = 0; i < n; i++) {
		scanf("%d%d", &x, &y);
		p[i].x = x;
		p[i].y = y;
	}
	build(0, n - 1, 0);

	sort(p, p + n);
	for (int i = 0; i < n; i++)
		inser(i, 0);

	memset(mrk, 0, sizeof(mrk));
	//printf("%d %d %d %d\n", nd[0].xma, nd[0].xmi, nd[0].yma, nd[0].ymi);
	for (int o = 0; o < 3; o++) {
		sq = MAX;
		for (int i = 0; i < n; i++) {
			if (mrk[i])
				continue;
			dele(i, 0);
			if (nd[0].xma <= nd[0].xmi || nd[0].yma <= nd[0].ymi) {
				printf("0\n");
				return 0;
			}

			int tp = (nd[0].xma - nd[0].xmi)*(nd[0].yma - nd[0].ymi);
			if (tp < sq) {
				sq = tp;
				mk = i;
			}
			inser(i, 0);

		}
		dele(mk, 0);
		mrk[mk] = 1;

		//printf("%d %d %d %d\n", nd[0].xma, nd[0].xmi, nd[0].yma, nd[0].ymi);
	}
	int tp = (nd[0].xma - nd[0].xmi)*(nd[0].yma - nd[0].ymi);
	printf("%d\n", tp);

	return 0;
}
*/

//granpa'estate
/*
#include<iostream>
#include<algorithm>
#include<memory.h>
#include<vector>
#include<string>
#include<set>
#include<math.h>
#define EPS 1e-6

using namespace std;


struct point{
	double x, y;
	bool v;
	point():v(-1){}
	point(double i, double j){
		x = i; y = j;
		v = -1;
	}
	bool operator<(const point p)const{
		if (x < p.x)
			return 1;
		else if (x > p.x)
			return 0;
		else
			return y < p.y;
	}
	point operator+(point p){
		return point(x + p.x, y + p.y);
	}
	point operator-(point p){
		return point(x - p.x, y - p.y);
	}
	point operator*(double a){
		return point(x*a, y*a);
	}
	point &operator=(point p){
		x = p.x; y = p.y;
		v = p.v;
		return *this;
	}
};

double operator*(point p1, point p2){
	return p1.x*p2.x + p1.y*p2.y;
}

double operator^(point p1, point p2){
	return p1.x*p2.y - p2.x*p1.y;
}

double distance(const point& p1, const point& p2){
	return sqrt((p1.x - p2.x)*(p1.x - p2.x)
		+ (p1.y - p2.y)*(p1.y - p2.y));
}

int sign(double x){
	return fabs(x) < EPS ? 0 : x > 0 ? 1 : -1;
}

void graham(vector<point> &points, vector<point> &convex){
	convex.clear();
	if (points.size() < 6)
		return;
	int n = points.size();
	int si;
	
	sort(points.begin(), points.end());
	points[0].v = 1;
	points[n - 1].v = 1;
	convex.push_back(points[0]);
	convex.push_back(points[1]);
	

	for (int i = 2; i < n; i++){
		while (convex.size()>1){
			point p2 = *(convex.end() - 1);
			point p1 = *(convex.end() - 2);

			si = sign((p2 - p1) ^ (points[i] - p2));
			if (si < 0)
				convex.pop_back();

			else{
				if (si > 0)
					(*(convex.end() - 1)).v = 1;
				else
					(*(convex.end() - 1)).v = 0;

				break;
			}
		}
		convex.push_back(points[i]);
	}

	int siz = convex.size();
	convex.push_back(points[n - 2]);
	for (int i = n - 3; i >= 0; i--){
		while (convex.size() > siz){
			point p2 = *(convex.end() - 1);
			point p1 = *(convex.end() - 2);

			si = sign((p2 - p1) ^ (points[i] - p2));
			if (si < 0)
				convex.pop_back();
			else{
				if (si > 0)
					(*(convex.end() - 1)).v = 1;
				else
					(*(convex.end() - 1)).v = 0;

				break;
			}
		}
		convex.push_back(points[i]);
	}
	convex.pop_back();
}

bool deter(vector<point>& p){
	int l, v;
	l = v = 0;
	int siz = p.size();
	if (!siz)
		return 0;
	
	for (int i = 0; i < siz; i++){

		if (p[i].v)
			v++;
		else if (i>0 && p[i - 1].v)
			l++;
	}

	return l == v;
}
int main(){
	int i, j, k;
	int n, t;
	vector<point> ps;
	vector<point> cv;
	double x, y;

	scanf("%d", &t);
	for (i = 0; i < t; i++){
		scanf("%d", &n);
		ps.clear();
		for (j = 0; j < n; j++){
			scanf("%lf%lf", &x, &y);
			ps.push_back(point(x, y));
		}

		graham(ps, cv);
		if (deter(cv))
			printf("YES\n");
		else
			printf("NO\n");
	}
}
*/

//the most distant point from the sea
//冗长
/*
#include<iostream>
#include<algorithm>
#include<memory.h>
#include<vector>
#include<string>
#include<set>
#include<math.h>
#define EPS 1e-6

using namespace std;


struct point{
	double x, y;
	point(){}
	point(double i, double j){
		x = i; y = j;
	}
	bool operator<(const point p)const{
		if (x < p.x)
			return 1;
		else if (x > p.x)
			return 0;
		else
			return y < p.y;
	}
	bool operator==(const point p)const{
		return fabs(x - p.x) < EPS&&fabs(y - p.y) < EPS;
	}
	point operator+(point p){
		return point(x + p.x, y + p.y);
	}
	point operator-(point p){
		return point(x - p.x, y - p.y);
	}
	point operator*(double a){
		return point(x*a, y*a);
	}
	point &operator=(point p){
		x = p.x; y = p.y;
		return *this;
	}
};

struct seg{
	point a, b;
	seg(const point aa, const point bb) :a(aa), b(bb){}
};

typedef vector<point> polygon;

double operator*(point p1, point p2){
	return p1.x*p2.x + p1.y*p2.y;
}

double operator^(point p1, point p2){
	return p1.x*p2.y - p2.x*p1.y;
}

double length(point l){
	return sqrt(l*l);
}

double distance(point p, seg s){
	return fabs((p - s.a) ^ (s.b - s.a)) / length(s.b - s.a);
}

double area(point p, point q){
	return (p^q);
}

int sign(double x){
	return fabs(x) < EPS ? 0 : x > 0 ? 1 : -1;
}

bool isz(double t){
	return t > -EPS&&t < EPS;
}
bool lessoreq(double a, double b){
	return b - a>EPS || isz(b - a);
}

point f(point p, double d){
	return point(-p.y, p.x)*(d / length(p));
}

bool pointinseg(point p, seg s){
	double tmp = (s.a - p) ^ (s.a - s.b);
	if (!isz(tmp))
		return 0;
	if (lessoreq(min(s.a.x, s.b.x), p.x) &&
		lessoreq(p.x, max(s.a.x, s.b.x)) &&
		lessoreq(min(s.a.y, s.b.y), p.y) &&
		lessoreq(p.y, max(s.a.y, s.b.y)))
		return 1;
	return 0;
}

pair<int, point> crosspoint(seg s1, seg s2){
	point p1 = s1.a;
	point p2 = s1.b;
	point p3 = s2.a;
	point p4 = s2.b;

	double a1 = area(p3 - p1, p4 - p1);
	double a2 = area(p4 - p2, p3 - p2);

	if (((p2 - p1) ^ (p3 - p1))*((p2 - p1) ^ (p4 - p1)) < -EPS &&
		((p4 - p3) ^ (p1 - p3))*((p4 - p3) ^ (p2 - p3)) < -EPS)
		return make_pair(0, p1 + ((p2 - p1)*(a1 / (a1 + a2))));

	if (!isz((p2 - p1) ^ (p3 - p4))){
		if (p1 == p3 || p1 == p4)
			return make_pair(1, p1);
		if (p2 == p3 || p2 == p4)
			return make_pair(1, p2);
		if (pointinseg(p1, s2))
			return make_pair(2, p1);
		if (pointinseg(p2, s2))
			return make_pair(2, p2);
		if (pointinseg(p3, s1))
			return make_pair(3, p3);
		if (pointinseg(p4, s1))
			return make_pair(3, p4);

		point crossp = p1 + ((p2 - p1)*(a1 / (a1 + a2)));

		if (pointinseg(crossp, s1))
			return make_pair(8, crossp);
		if (pointinseg(crossp, s2))
			return make_pair(9, crossp);
		return make_pair(4, crossp);
	}

	if (!isz(distance(p1, s2)))
		return make_pair(5, point(0, 0));
	if (pointinseg(p1, s2))
		return make_pair(6, p1);
	if (pointinseg(p2, s2))
		return make_pair(6, p2);
	if (pointinseg(p3, s1))
		return make_pair(6, p3);
	if (pointinseg(p4, s1))
		return make_pair(6, p4);
	return make_pair(7, point(0, 0));

}

//用直线a->b 切割src,返回其左侧的多边形，放到result
//src里的点什么顺序无所谓 
int cutpolygon(const polygon &src, point a, point b, polygon& result){
	int n = src.size();
	result.clear();
	for (int i = 0; i < n; i++){
		point c = src[i];
		point d = src[(i + 1) % n];
		if (sign((b - a) ^ (c - a)) >= 0){
			result.push_back(c);
		}
		pair<int, point> r = crosspoint(seg(c, d), seg(a, b));
		if (r.first == 0 || r.first == 8 || r.first == 3)
			result.push_back(r.second);
	}
	return result.size();
}

//切割，返回点数
int ope(polygon& p, double dis){
	int i, j, k;
	int siz = p.size();
	point cp;
	polygon res, tp;
	polygon::iterator it;
	for (it = p.begin(); it != p.end(); it++)
		tp.push_back(*it);

	for (i = 0; i < siz; i++){
		cp = f(p[(i + 1)%siz] - p[i%siz], dis);
		cutpolygon(tp, p[i%siz] + cp, p[(i + 1)%siz] + cp, res);
		if (!res.size())
			return 0;
		tp.clear();
		for (it = res.begin(); it != res.end(); it++)
			tp.push_back(*it);
	}

	return res.size();
}

double mostdis(polygon& p){
	double l = 0;
	double r = 10000;
	double m ;
	int k;

	while (!(r-l<EPS)){
		m = (l + r) / 2;
		k = ope(p, m);
		if (k == 1)
			return m;
		if (k > 1)
			l = m;
		else
			r = m;
	}
	m = (l + r) / 2;
	return m;
}

int main(){
	int n, i, j, k;
	double x, y;
	polygon isl;
	while (scanf("%d", &n)){
		if (!n)
			return 0;
		isl.clear();
		for (i = 0; i < n; i++){
			scanf("%lf%lf", &x, &y);
			isl.push_back(point(x, y));
		}
		printf("%lf\n", mostdis(isl));
	}
	return 0;
}
*/

//pipe
//****** 注意一般不能直接用'<', '>'比较浮点类数大小 ******
/*
#include<iostream>
#include<algorithm>
#include<memory.h>
#include<vector>
#include<string>
#include<set>
#include<math.h>
#define EPS 1e-6
#define MI -1e9

using namespace std;


struct point{
	double x, y;
	point(){}
	point(double i, double j){
		x = i; y = j;
	}
	
	point operator+(point p){
		return point(x + p.x, y + p.y);
	}
	point operator-(point p){
		return point(x - p.x, y - p.y);
	}
	point operator*(double a){
		return point(x*a, y*a);
	}
	point &operator=(point p){
		x = p.x; y = p.y;
		return *this;
	}
};

struct seg{
	point a, b;
	seg(point aa, point bb) :a(aa), b(bb){}
};

double operator^(point p1, point p2){
	return p1.x*p2.y - p2.x*p1.y;
}

double area(point p, point q){
	return (p^q);
}

bool isz(double t){
	return t > -EPS&&t < EPS;
}

int sign(double x){
	return fabs(x) < EPS ? 0 : x > 0 ? 1 : -1;
}

point cp(seg l, seg m){
	double x = area(m.a - l.a, l.b - l.a);
	double y = area(l.b - l.a, m.b - l.a);
	if (isz(x + y))
		return point(MI, MI);
	return m.a + (m.b - m.a)*(x/(x+y));
}




pair<int, double> f(vector<point> &v1, vector<seg> &v2){
	int i, j, k;
	double x_max=v2[0].a.x;
	int siz = v1.size();
	int ssiz = v2.size();
	for (i = 0; i < siz - 1; i++){
		for (j = i + 1; j < siz; j++){
			if ((i>>1)==(j>>1))
				continue;
			seg s(v1[i], v1[j]);
			for (k = 0; k < ssiz; k++){
				point pa = cp(s, v2[k]);
				
				if (sign(pa.y - v2[k].a.y) < 0 || sign(v2[k].b.y - pa.y) < 0){
					if (k){
						seg s1(v2[k - 1].a, v2[k].a);
						seg s2(v2[k - 1].b, v2[k].b);
						point p1 = cp(s, s1);
						point p2 = cp(s, s2);

						x_max = max(max(p1.x, p2.x), x_max);
					}
					break;
				}
			}
			if (k == ssiz)
				return make_pair(1, 0);
		}
	}
	
	return make_pair(0, x_max);
}

int main(){
	int n, i, j;
	double x, y;
	vector<point> pi;
	vector<seg> ss;
	while (scanf("%d", &n)){
		if (!n)
			break;
		pi.clear();
		ss.clear();
		for (i = 0; i < n; i++){
			scanf("%lf%lf", &x, &y);
			point a(x, y-1);
			point b(x, y);
			pi.push_back(a);
			pi.push_back(b);
			ss.push_back(seg(a,b));
		}
		pair<int, double> pp = f(pi, ss);
		if (pp.first)
			printf("Through all the pipe.\n");
		else
			printf("%.2lf\n", pp.second);

	}
	return 0;
}
*/

/*
#include<iostream>
#include<algorithm>
#include<stack>
#define N 50005
using namespace std;

int par[N];

int getroot(int i){
	if (par[i] == i)
		return i;

	int k = getroot(par[i]);
	par[i] = k;
	return k;
}

int merge(int i, int j){
	int p = getroot(i), q = getroot(j);
	if (p == q)
		return 0;

	par[q] = p;
	return 1;
}


int main(){
	int n, m, i, j, k;
	int t = 1;
	while (scanf("%d%d", &n, &m)){
		if (n == 0 && m == 0)
			return 0;

		for (i = 1; i <= n; i++)
			par[i] = i;
		for (i = 0; i < m; i++){
			scanf("%d%d", &j, &k);
			n -= merge(k, j);
		}
		printf("Case %d: %d\n", t++, n);
	}
}


*/

/*
#include<iostream>
#include<algorithm>
#include<memory.h>
#define N 500005
using namespace std;

struct num{
	int v;
	int ori;
	bool operator<(num &n){
		if (v == n.v)
			return ori < n.ori;
		return v < n.v;
	}
}nu[N];

struct pp{
	int ori;
	int rank;
	bool operator<(pp &p){
		return ori < p.ori;
	}
}ra[N];


int lowbit(int i){
	return i&(-i);
}

int c[N];

int query(int k){
	int res = 0;
	for (int p = k; p > 0; p -= lowbit(p)){
		res += c[p];
	}
	return res;
}

void inser(int k, int lim){
	for (int p = k; p <= lim; p+=lowbit(p)){
		c[p]++;
	}
}


int main(){
	int n, i, j, k;
	while (scanf("%d", &n)){
		if (n == 0)
			return 0;


		memset(c, 0, sizeof(c));
		for (i = 0; i < n; i++){
			scanf("%d", &j);
			nu[i].ori = i;
			nu[i].v = j;
		}
		sort(nu, nu + n);

		for (i = 0; i < n; i++){
			ra[i].rank = i;
			ra[i].ori = nu[i].ori;
		}

		sort(ra, ra + n);
		long long res = 0;

		for (i = n - 1; i >= 0; i--){
			res += query(ra[i].rank);
			inser(ra[i].rank + 1, n);
		}
		printf("%lld\n", res);
	}
}
*/

/*
#include<iostream>
#include<algorithm>
#include<memory.h>
#define N 100005

using namespace std;

long long res;
struct node{
	int b, e;
	long long val;
	int inc;
}nd[4*N];

void build(int k, int b, int e){
	nd[k].b = b;
	nd[k].e = e;
	nd[k].val = 0;
	nd[k].inc = 0;
	if (b == e)
		return;
	int m = (b + e) / 2;
	build(2 * k + 1, b, m);
	build(2 * k + 2, m + 1, e);
}

void adder(int b, int e, int val, int k){
	if (b <= nd[k].b&&e >= nd[k].e){
		nd[k].inc += val;
		return;
	}
	int l = max(b, nd[k].b);
	int r = min(e, nd[k].e);

	if (l <= r)
		nd[k].val += (long long)(r - l + 1)*(long long)(val);

	if (b < nd[2 * k + 2].b)
		adder(b, e, val, 2 * k + 1);
	if (e>nd[2 * k + 1].e)
		adder(b, e, val, 2 * k + 2);
}

void sum(int b, int e, int k){
	if (b <= nd[k].b&&e >= nd[k].e){
		res += (long long)(nd[k].e - nd[k].b + 1)*
			(long long)(nd[k].inc) + nd[k].val;
		return;
	}
	nd[2 * k + 1].inc += nd[k].inc;
	nd[2 * k + 2].inc += nd[k].inc;
	nd[k].val += (long long)(nd[k].e - nd[k].b + 1)*
		(long long)(nd[k].inc);
	nd[k].inc = 0;

	if (b < nd[2 * k + 2].b)
		sum(b, e, 2 * k + 1);
	if (e>nd[2 * k + 1].e)
		sum(b, e, 2 * k + 2);
}

int main(){
	int n, i, j, k, q,l;
	char c;
	scanf("%d%d", &n, &q);
	build(0, 1, n);
	for (int i = 0; i < n; i++){
		scanf("%d", &j);
		adder(i + 1, i + 1, j, 0);
	}

	for (int i = 0; i < q; i++){
		getchar();
		scanf("%c", &c);
		if (c == 'Q'){
			res = 0;
			scanf("%d%d", &j, &k);
			sum(j, k, 0);
			printf("%lld\n", res);
		}
		else{
			scanf("%d%d%d", &j, &k, &l);
			adder(j, k, l, 0);
		}
	}
}
*/

/*
#include<iostream>
#include<cstring>
#include<queue>
#define N 120005

using namespace std;

int nodecnt;
struct node{
	int prev;
	int child[26];
	bool bad;

	node(){
		prev = -1;
		bad = 0;
		memset(child, 0xff, sizeof(child));
	}
}nd[N];


void inse(char* s);
void builddf();
bool cp(char *s);


int main(){
	int n, i, j, k, m;
	char s[125], ss[1005];
	scanf("%d", &n);
	nodecnt = 1;
	for (i = 0; i < n; i++){
		cin >> s;
		inse(s);
	}
	builddf();
	scanf("%d", &m);
	for (i = 0; i < m; i++){
		cin >> ss;
		if (cp(ss))
			printf("YES\n");
		else
			printf("NO\n");
	}
}


void inse(char *s){
	int k = 1;
	for (int i = 0; s[i]; i++){
		if (nd[k].child[s[i] - 'a'] < 0)
			nd[k].child[s[i] - 'a'] = ++nodecnt;
		k = nd[k].child[s[i] - 'a'];
	}
	nd[k].bad = 1;
}

void builddf(){
	for (int i = 0; i < 26; i++){
		nd[0].child[i] = 1;
	}
	nd[0].prev = -1;
	nd[1].prev = 0;
	queue<int> q;
	q.push(1);
	while (!q.empty()){
		int pr = q.front();
		q.pop();

		for (int i = 0; i < 26; i++){
			int tp = nd[pr].child[i];
			if (tp > 0){
				int pp = nd[pr].prev;
				while (nd[pp].child[i] < 0)
					pp = nd[pp].prev;

				nd[tp].prev = nd[pp].child[i];
				if (nd[nd[tp].prev].bad)
					nd[tp].bad = 1;

				q.push(tp);
			}
		}
	}
}

bool cp(char *s){
	int k = 1;
	for (int i = 0; s[i]; i++){
		while (nd[k].child[s[i] - 'a'] < 0)
			k = nd[k].prev;

		k = nd[k].child[s[i] - 'a'];
		if (nd[k].bad)
			return 1;
	}
	return 0;
}
*/

//LCS(O(nlogn)dynamic programming)
/*
#include<stdio.h>
#include<algorithm>
#include<memory.h>

int a[1005];

int binser(int b,int e,int k){
	int l = b, r = e;
	int m;
	while (l != r){
		m = (l + r) / 2;
		if (k > a[m])
			l = m + 1;
		else
			r = m;
	}
	return l;
}

int main(){
	int n, m, p;
	int len = 0;
	memset(a, 0x7f, sizeof(a));
	a[0] = -1;
	scanf("%d", &n);
	for (int i = 0; i < n; i++){
		scanf("%d", &m);
		if (m>a[len]){
			len++;
			a[len] = m;
		}
		else{
			p = binser(0, len, m);
			if (m < a[p])
				a[p] = m;
		}
	}
	
	printf("%d\n", len);
}
*/

/* the coin */
/*
#include<stdio.h>
#include<algorithm>
using namespace std;
int f[10005] = { 1 };
int a[205];
int m[205] = { 0 };

//the number of ways to consist sum r without element p
int ser(int p, int r){
	if (r < p){
		return f[r];
	}
	return f[r] - ser(p, r - p);
}

int main(){
	int i, j, x, n, c;
	scanf("%d%d", &n, &x);
	for (i = 0; i < n; i++){
		scanf("%d", a + i);
		for (j = x; j >= a[i]; j--){
			f[j] += f[j - a[i]];
		}
	}
	sort(a, a + n);
	c = 0;
	for (i = 0; i < n; i++){
		if (!(m[i] = ser(a[i], x)))
			c++;
	}
	printf("%d\n", c);
	for (i = 0; i < n; i++){
		if (!m[i])
			printf("%d ", a[i]);
	}
	printf("\n");
	return 0;
}
*/

/*
#include<stdio.h>
#include<math.h>

using namespace std;

double x[10], y[10];
double avex=0;
double avey = 0;

double sxy;
double sxx;

double sx;
double sy;

int main(){
	for (int i = 0; i < 10; i++){
		scanf("%lf", x + i);
		avex += x[i];
	}
	for (int i = 0; i < 10; i++){
		scanf("%lf", y + i);
		avey += y[i];
	}
	avex /= 10;
	avey /= 10;

	sxy = sxx = 0;

	sy = 0;

	for (int i = 0; i < 10; i++){
		sxx += (x[i] - avex)*(x[i] - avex);
		sxy += (x[i] - avex)*(y[i] - avey);
		sy += (y[i] - avey)*(y[i] - avey);
	}

	printf("%lf\n", sqrt(((2 * sxy*sxy) / sxx - (sxy / sxx)*(sxy / sxx)*sxx) / sy));
}
*/

/*
#include<stdio.h>

int t[4] = { 53,8,68,24 };

void pro(int q) {
	int temp[4], i, j, k;
	int round = 0;
	int sum = 0;
	int waitime[4] = { 0 };
	int decinc;

	for (i = 0; i < 4; i++) {
		temp[i] = t[i];
		sum += temp[i];
	}

	//模拟轮转
	while (sum) {
		decinc = (temp[round % 4] < q ? temp[round % 4] : q);
		temp[round % 4] -= decinc;
		for (i = 0; i < 4; i++)
			if (i != round % 4 && temp[i])
				waitime[i] += decinc;

		sum -= decinc;
		round++;

	}

	double ave = 0;
	for (i = 0; i < 4; i++) {
		printf("%d ", waitime[i]);
		ave += waitime[i];
	}
	printf("%lf\n", ave / 4.0);
	return;
}

int main() {
	int q, i, j;
	while (scanf("%d", &q)) {
		if (!q)return 0;
		pro(q);
	}

	return 0;
}
*/

/*
#include<stdio.h>
#include<string.h>


#define EXP 1e-6









bool eq(double x, double y) {
	return (-EXP < x - y) && (x - y < EXP);
}

bool chk(int t, double* p) {
	if (t == 1)return eq(p[0], 24);

	bool res = 0;
	
	int mrk[4] = { 0 };
	for (int i = 0; i < t; i++) {
		mrk[i] = 1;
		for (int j = 0; j < t; j++) {
			if (mrk[j])continue;
			mrk[j] = 1;
			double np[4];
			int f = 0;
			for (int k = 0; k < t; k++)
				if (!mrk[k])np[f++] = p[k];

			
			np[f] = p[i] + p[j];
			res = res || chk(t - 1, np);
			np[f] = p[i] - p[j];
			res = res || chk(t - 1, np);
			np[f] = p[i] * p[j];
			res = res || chk(t - 1, np);
			np[f] = p[i] / p[j];
			res = res || chk(t - 1, np);

			mrk[j] = 0;
		}
		mrk[i] = 0;
	}

	return res;
}

int main() {
	double d[4];
	while (scanf("%lf%lf%lf%lf", d, d + 1, d + 2, d + 3)) {
		if (eq(d[0], 0))return 0;
		if (chk(4, d))printf("YES\n");
		else printf("NO\n");
	}
}
*/

/*
#include<stdio.h>

#include<string.h>


char ch[205];
int lv[205];
int left[205];
int right[205];

void fr(int nd) {
	printf("%c", ch[nd]);
	if (left[nd]>0)fr(left[nd]);
	if (right[nd]>0)fr(right[nd]);
}

void mi(int nd) {
	if (left[nd]>0)mi(left[nd]);
	printf("%c", ch[nd]);
	if (right[nd]>0)mi(right[nd]);
}

void bh(int nd) {
	if (left[nd]>0)bh(left[nd]);
	if (right[nd]>0)bh(right[nd]);
	printf("%c", ch[nd]);
}



int main() {
	char s[105];
	int las[105];
	int cnt, i, j;
	int ni;
	scanf("%d", &cnt);
	for (i = 0; i < cnt; i++) {

		memset(las, 0xff, sizeof(las));
		memset(left, 0xff, sizeof(left));
		memset(right, 0xff, sizeof(right));
		ni = 0;

		while (scanf("%s", s)) {
			if (s[0] == '0')break;
			int l = 0;
			while (s[l] == '-')l++;
			ch[ni] = s[l];
			lv[ni] = l;


			if (s[l] != '*') {
				las[l] = ni;
				if (l) {
					int ft = las[l - 1];
					if (ft + 1 != ni && lv[ft + 1] == l)right[ft] = ni;
					else left[ft] = ni;
				}
			}

			ni++;
		}

		fr(0); printf("\n");
		bh(0); printf("\n");
		mi(0); printf("\n");
		printf("\n");


	}
}
*/

/*
#include<stdio.h>

#define M 100000000

bool isp(int k) {
	int r = 0;
	for (int i = 2; (i * i) <= k; i++) {
		if (k % i==0) {
			r = 1; break;
		}
	}
	return (r == 0);
}

int main() {
	int sum = 0;
	for (int i = 2; i <= M; i++) {
		if (isp(i)) {
			//printf("%d\n", i);
			sum++;
		}
	}
	printf("%d\n", sum);
}
*/

/*
#include<stdio.h>


int main() {
	int n, m, i, j, k;
	int f[105][105] = { 0 };
	int p[105] = { 0 };
	int tmp;
	scanf("%d%d", &n, &m);
	for (i = 1; i <= n; i++) {
		scanf("%d", p + i);
		f[i][1] = p[i];
	}
	for (i = 1; i <= n; i++) {
		for (j = 2; j <= i; j++) {
			tmp = 0;
			for (k = j - 1; k < i; k++)
				if (tmp < f[k][j - 1] + p[i - k])
					tmp = f[k][j - 1] + p[i - k];

			f[i][j] = tmp;
		}
	}


	printf("%d\n", f[n][m]);
}
*/


//Floyd算法
// f[i][j][k] = min{f[i][j][k-1],f[i][k][k-1]+f[k][j][k-1]}
// f[i][j][k]为从i到j，途中经过的点的index<=k的最短路径长
/*
#include<stdio.h>
#include<string.h>

int f[12][12][12] = {0};
int r[12][12] = { 0 };

void path(int i,int j) {
	int k;
	if (k=r[i][j]) {
		path(i, k);
		printf("%d ", k);
		path(k, j);
	}
}

int main() {
	int n;

	memset(f, 0x7f, sizeof(f));
	scanf("%d", &n);
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++) {
			scanf("%d", &f[i][j][0]);
		}

	for(int k=1;k<=n;k++)
		for(int i=1;i<=n;i++)
			for (int j = 1; j <= n; j++)
			{
				int a = f[i][k][k - 1] + f[k][j][k - 1],
					b = f[i][j][k - 1];
				
				f[i][j][k] = a < b ? a : b;
				r[i][j] = a < b ? k : r[i][j];
					
			}
			

	printf("%d\n", f[1][n][n-1]);
	printf("%d ", 1);
	path(1,n);
	printf("%d\n", n);
}
*/

/*
#include<stdio.h>
#include<string.h>

int f[12][12][12] = { 0 };

int main() {
	int n, K, m;

	memset(f, 0x7f, sizeof(f));
	scanf("%d%d", &n, &K);
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++) {
			scanf("%d", &f[i][j][0]);
		}

	for (int k = 1; k <= n; k++)
		for (int i = 1; i <= n; i++)
			for (int j = 1; j <= n; j++)
			{
				int a = f[i][k][k - 1] + f[k][j][k - 1],
					b = f[i][j][k - 1];

				f[i][j][k] = a < b ? a : b;

			}


	for (m = n; m > 0;m--) {
		if (f[1][m][n] <= K)
			break;
	}
	
	printf("%d\n", m);
}
*/

/*
#include<stdio.h>
#define N 1005


int f[N][N];

int main() {
	int n, K, m;
	
	scanf("%d%d", &n, &K);
	m = K;

	for (int i = 1; i <= n; i++) {
		scanf("%d", &f[i][0]);
	}
	
	for(int i=1;i<=n;i++)
		for (int j = 1; j <= n; j++) {
			f[i][j] = K - f[i][0] - f[j][0];
			if (i!=j && f[i][j] >= 0 && f[i][j] < m)
				m = f[i][j];
		}

	printf("%d\n", K - m);
}

*/

//多重背包问题？
/*
#include<stdio.h>
#include<string.h>
#define N 100005

int main() {
	int n, c, w[105], v[105], m = 0;

	int f[1005];
	memset(f, 0xff, sizeof(f));
	f[0] = 0;

	scanf("%d%d", &n, &c);
	for (int i = 0; i < n; i++)
		scanf("%d", w + i);
	for (int i = 0; i < n; i++)
		scanf("%d", v + i);

	for (int i = 0; i < n; i++)
		for (int j = w[i]; j <= c; j++) {
			if (f[j - w[i]] >= 0) {
				if (f[j - w[i]] + v[i] > f[j])
					f[j] = f[j - w[i]] + v[i];
			}
			if (f[j] > m)
				m = f[j];
		}

	printf("%d\n", m);

}
*/

/*
#include<stdio.h>

int main() {
	int m = 0;
	int a, b, c, d, e;
	int s[2] = {0,0};

	scanf("%d%d%d%d%d", &a, &b, &c, &d, &e);

	for(int i=0;i<=100;i++)
		for (int j = 0; j <= 100; j++) {
			if (2 * i <= b && 75 * j <= e &&
				250 * i + 200 * j <= a && 75 * i + 150 * j <= c &&
				100 * i + 150 * j <= d) {
				if (400 * i + 450 * j > m) {
					m = 400 * i + 450 * j;
					s[0] = i; s[1] = j;
				}
			}
		}

	printf("%d\n%d\n%d\n", m,s[0],s[1]);
}
*/

/*
#include<stdio.h>
unsigned f(unsigned n) {
	unsigned temp = n, r = 0, p = 0, q = 0;
	while (!(temp & 1)) {
		temp = temp >> 1;
		p++;
	}
	while (temp & 1) {
		temp = temp >> 1;
		q++;
	}

	temp += 1;
	temp = temp << (p + q);
	if (q > 1) {
		r = n << (33 - p - q);
		r = r >> (33 - q);
	}

	return temp | r;
}

int main() {
	unsigned n;
	
	while (scanf("%d", &n)) {
		if (!n)break;
		printf("%d\n", f(n));
	}
}
*/

/*
#include<stdio.h>
#include<math.h>
#include<queue>
using namespace std;
int main() {
	queue<int> q1, q2;
	int m, n, k;
	scanf("%d%d", &n, &m);
	for(int i=0;i<n;i++){
		scanf("%d", &k);
		q1.push(k);
	}

	for (int i = 0; i < n - m; i++) {
		k = q1.front();
		q1.pop();
		q2.push(k);
	}

	while (!q1.empty()) {
		k = q1.front();
		q1.pop();
		printf("%d ", k);
	}

	while (!q2.empty()) {
		k = q2.front();
		q2.pop();
		printf("%d ", k);
	}
}
*/

/*
#include"stdio.h"
#include"string.h"

int g[20005] = { 0 };

long fil(int l,int r) {
	int i;
	long res = 0;
	while (l < r) {
		if (g[l] <= g[r]) {
			for (i = l; i==l||g[i] < g[l]; i++) {
				res += g[l] - g[i];
			}
			l = i;
		}
		else {
			for (i = r; i==r||g[i] < g[r]; i--) {
				res += g[r] - g[i];
			}
			r = i;
		}
	}
	
	return res;
}

int main() {
	int m, n, i, j, l, r;
	scanf("%d", &m);

	while (m--) {
		scanf("%d", &n);

		for (i = 0; i < n; i++) {
			scanf("%d", g + i);
		}
	
		printf("%ld\n", fil(0, n - 1));
		
	}
}*/

/*
#include<stdio.h>
#include<queue>
#include<string.h>
using namespace std;

struct way {
	
	int board[4][4];
	queue<int> q;
	way& operator = (way w) {
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				board[i][j] = w.board[i][j];

		while (!q.empty())q.pop();
		while (!w.q.empty()) {
			int a = w.q.front();
			w.q.pop();
			q.push(a);
		}
		return *this;
	}
};

int to[4][2] = { -1,0,0,-1,1,0,0,1 };
int mrk[4][4];

int use[65536] = { 0 };

int trans(way& w) {
	int res=0;
	for(int i=0;i<4;i++)
		for (int j = 0; j < 4; j++) {
			res = res << 1;
			res = res | w.board[i][j];
		}
	return res;
}

int main() {
	int d[2][4][4];
	int i, j, a, k;
	way w, temp;
	queue<way> qq;
	
	for (k = 0; k < 2; k++) {
		for (i = 0; i < 4; i++) {
			scanf("%d", &a);
			for (j = 0; j < 4; j++) {
				d[k][i][3 - j] = a % 10;
				a = a / 10;
			}
		}
	}

	

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			w.board[i][j] = d[0][i][j];
	while (!w.q.empty())w.q.pop();
	use[trans(w)] = 1;
	qq.push(w);
	while (!qq.empty()) {
		
		w = qq.front();
		qq.pop();


		for (i = 0; i < 4; i++) {
			for (j = 0; j < 4; j++)
				if (w.board[i][j] != d[1][i][j])
					break;
			if (j < 4)break;
		}


		if (i == 4 && j == 4) break;

		memset(mrk,0, sizeof(mrk));

		//move
		for(i=0;i<4;i++)
			for (j = 0; j < 4; j++) {
				for (k = 0; k < 4; k++) {
					int ii = i + to[k][0], jj = j + to[k][1];
					if (ii >= 0 && ii < 4 && jj >= 0 && jj < 4) {
						if (!mrk[ii][jj] && w.board[i][j] != w.board[ii][jj]) {
							temp = w;
							temp.board[i][j] = w.board[ii][jj];
							temp.board[ii][jj] = w.board[i][j];

							
							temp.q.push(i);
							temp.q.push(j);
							temp.q.push(ii);
							temp.q.push(jj);

							int t = trans(temp);
							if (!use[t]) {
								use[t] = 1;
								qq.push(temp);
							}


						}
					}
				}
				mrk[i][j] = 1;
			}
		//printf("%d\n", qq.size());
	}
	

	printf("%d\n", w.q.size() / 4);

	for (int i = 0; !w.q.empty(); i++) {
		printf("%d", w.q.front()+1);
		w.q.pop();
		if (i % 4 == 3)printf("\n");
		else printf(" ");
	}

	
}
*/

//上poj(openjudge不行)，广搜+数据压缩技巧
/*
#include"stdio.h"
#include"stack"
#include"queue"
#define N 1<<16

using namespace std;

int use[N] = { 0 };
int pre[N] = { 0 };
int to[4][2] = { 0,-1,-1,0,0,1,1,0 };

stack<int> st;

int tran(int a, int i, int j, int bia) {
	int ii = i + to[bia][0], jj = j + to[bia][1];
	if (ii < 0 || ii >= 4 || jj < 0 || jj >= 4)return  -1;
	int s = 15 - 4 * i - j, d = 15 - 4 * ii - jj;

	int p1 = ((1 << s) & a) >> s, p2 = ((1 << d) & a) >> d;
	if (p1==p2)return -1;

	int c = ~((1 << s) | (1 << d));

	return ((a & c) | (p1 << d) | (p2 << s));
}

void pick(int a,int b) {
	int u = a ^ b;
	int p1 = 0, p2;
	while (!((u >> p1) & 1))p1++;
	p2 = p1 + 1;
	while (!((u >> p2) & 1))p2++;

	st.push((15 - p1) % 4);
	st.push((15 - p1) / 4);
	st.push((15 - p2) % 4);
	st.push((15 - p2) / 4);
}

int main() {
	int d[2] = { 0 };
	int i, j, k, a, t;
	queue<int> q;

	for (k = 0; k < 2; k++) {
		for (i = 0; i < 4; i++) {
			d[k] = d[k] << 4;
			scanf("%d", &a);
			int t = 0;
			for (j = 0; j < 4; j++) {
				t = t |((a % 10)<<j);
				a = a / 10;
			}
			d[k] = d[k] | t;
		}
	}

	//printf("%x %x\n", d[0], d[1]);

	q.push(d[0]);

	while (!q.empty()) {
		t = q.front();
		q.pop();
		if (t == d[1])break;

		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				for (k = 0; k < 4; k++) {
					int temp = tran(t, i, j, k);
					if (temp >= 0) {
						if (!use[temp]) {
							use[temp] = 1;
							pre[temp] = t;
							q.push(temp);
						}
					}

				}
	}

	while (t != d[0]) {
		pick(pre[t], t);
		t = pre[t];
	}

	printf("%d\n", st.size() >> 2);
	for (i = 0; !st.empty(); i++) {
		printf("%d ", st.top() + 1);
		st.pop();
		if (i % 4 == 3)printf("\n");
	}
}
*/

/*
#include<stdio.h>
#define N 1<<28

int mrk[N] = { 0 };

int main() {
	int k = 0, i = 1;
	while (k < 1500) {
		int temp = i;
		if (i <= N) {
			if ((i % 2 == 0 && mrk[i / 2]) || (i % 5 == 0 && mrk[i / 5]) || (i % 3 == 0 && mrk[i / 3])||i==1) {
				mrk[i] = 1;
			}
		}

		else {
			while (temp > N) {
				if (temp % 2 == 0)temp = temp / 2;
				else if (temp % 3 == 0)temp = temp / 3;
				else if (temp % 5 == 0)temp = temp / 5;
				else break;
			}
			if (temp < N)mrk[i] = mrk[temp];
		}

		if (mrk[i]) {
			printf("%d %d\n", ++k, i);
		}
		i++;

	}
}
*/

//ugly numbers
/*
#include"stdio.h"
#define M 0x7fffffff

int u[1505];
int p[3] = { 2,3,5 };

int main() {
	int i, j, ind;
	int k = 1;
	u[1] = 1;
	

	while (k < 1500) {
		int t = M;
		for (j = 0; j < 3; j++) {
			for (i = k; i > 0 && (p[j] * u[i]<0 || p[j] * u[i] > u[k]); i--);
			i++;
			if (p[j] * u[i] < t)t = p[j] * u[i];
		}
		u[++k] = t;
	}

	while (scanf("%d", &ind)) {
		if (!ind)break;
		printf("%d\n", u[ind]);
	}
}
*/

//线段树法（忘的差不多了）
/*
#include<stdio.h>
#include<algorithm>
#include<queue>
#define N 50005

using namespace std;

int a[N], b[N];

struct node {
	int left, right;
	int cnt, inc;
}nd[8*N];

struct T {
	int tim;
	int act;//0 for start. 1 for finish
	int cow;

	bool operator < (T& _t) {
		if (tim == _t.tim) {
			if (act == _t.act) {
				return cow < _t.cow;
			}
			return act < _t.act;
		}
		return tim < _t.tim;
	}

	void set(int t, int a, int c) {
		tim = t;
		act = a;
		cow = c;
	}
}tp[2*N];

void build(int ind, int s, int e) {
	nd[ind].cnt = nd[ind].inc = 0;

	if (s == e) {
		nd[ind].left = nd[ind].right = s;
		return;
	}
	nd[ind].left = s;
	nd[ind].right = e;
	int m = (s + e) / 2;
	build(2 * ind + 1, s, m);
	build(2 * ind + 2, m + 1, e);
}

void ser(int ind, int cowind, int* rd) {
	T t1, t2;
	t1.set(a[cowind], 0, cowind);
	t2.set(b[cowind], 1, cowind);

	//completely
	int ll = nd[ind].left, rr = nd[ind].right;
	if (!(tp[ll] < t1) && !(t2 < tp[rr])) {
		nd[ind].inc++;
		if (nd[ind].cnt+nd[ind].inc > * rd)
			*rd = nd[ind].cnt+nd[ind].inc;
	}

	else {
		nd[ind].cnt += nd[ind].inc;
		nd[2 * ind + 1].inc += nd[ind].inc;
		nd[2 * ind + 2].inc += nd[ind].inc;
		nd[ind].inc = 0;

		int lr = nd[2 * ind + 1].right;
		int rl = nd[2 * ind + 2].left;

		if (!(tp[lr] < t1)) {
			ser(2 * ind + 1, cowind, rd);
		}

		if (!(t2 < tp[rl])) {
			ser(2 * ind + 2, cowind, rd);
		}
		
		int tl = nd[2 * ind + 1].inc + nd[2 * ind + 1].cnt,
			tr = nd[2 * ind + 2].inc + nd[2 * ind + 2].cnt;

		if (nd[ind].cnt < tl)nd[ind].cnt = tl;
		if (nd[ind].cnt < tr)nd[ind].cnt = tr;

		if (*rd < nd[ind].cnt)*rd = nd[ind].cnt;
	}
}

void PR(int i) {
	printf("range:%d %d, sum:%d\n", nd[i].left, nd[i].right, nd[i].cnt + nd[i].inc);
	if (nd[i].left != nd[i].right) {
		PR(2 * i + 1);
		PR(2 * i + 2);
	}

}


int main() {
	int n, ma = 0;
	queue<int> q;

	scanf("%d", &n);

	for (int i = 0; i < n; i++) {
		scanf("%d%d", a+i, b+i);
		tp[2 * i].set(a[i], 0, i);
		tp[2 * i + 1].set(b[i], 1, i);
	}

	sort(tp, tp + 2 * n);

	build(0, 0, 2 * n - 1);

	for (int i = 0; i < n; i++) {
		ser(0, i, &ma);
		//printf("%d\n", ma);
	}

	//PR(0);

	//printf("%d\n", ma); return 0;

	for (int i = 1; i <= ma; i++)
		q.push(i);

	for (int i = 0; i < 2 * n; i++) {
		if (tp[i].act) {
			q.push(a[tp[i].cow]);
		}
		else {
			a[tp[i].cow] = q.front();
			q.pop();
		}
	}

	printf("%d\n", ma);

	for (int i = 0; i < n; i++) {
		printf("%d\n", a[i]);
	}
}

*/

/*
#include<stdio.h>
#include<algorithm>
#include<queue>
#define N 50005

using namespace std;

int a[N];



struct T {
	int tim;
	int act;//-1 for start. 1 for finish
	int cow;

	bool operator < (T& _t) {
		if (tim == _t.tim) {
			if (act == _t.act) {
				return cow < _t.cow;
			}
			return act < _t.act;
		}
		return tim < _t.tim;
	}

	void set(int t, int a, int c) {
		tim = t;
		act = a;
		cow = c;
	}
}tp[2 * N];






int main() {
	int n, ma = 0;
	int x, y, k;
	queue<int> q;

	scanf("%d", &n);

	for (int i = 0; i < n; i++) {
		scanf("%d%d", &x, &y);
		tp[2 * i].set(x, -1, i);
		tp[2 * i + 1].set(y, 1, i);
	}

	sort(tp, tp + 2 * n);

	for (int i = 0; i < 2 * n; i++) {
		int t = tp[i].cow;
		if (tp[i].act < 0) {
			if (q.empty()) {
				a[t] = ++ma;
			}
			else {
				a[t] = q.front();
				q.pop();
			}
		}
		else {
			q.push(a[t]);
		}

	}

	printf("%d\n", ma);

	for (int i = 0; i < n; i++) {
		printf("%d\n", a[i]);
	}
}
*/


/*
#include<stdio.h>
#include<string.h>
#include<algorithm>
#include<malloc.h>
#define N 5000
#define eps 0.0001
#define stdic 0.065

using namespace std;

char str[N];
int preserve[N][N];
int biaK[N] = {0};

bool testk(char* str, int step) {

	int length = strlen(str);
	int num[26] = { 0 };
	int pl = 0;	//part of ...
	double total = 0, single;

	memset(preserve, 0, sizeof(preserve));

	for (int i = 0; i < step; i++) {
		memset(num, 0, sizeof(num));
		pl = 0;
		single = 0;
		for (int j = i; j < length; j += step) {
			num[str[j] - 'A']++;
			preserve[i][pl] = str[j] - 'A';
			pl++;
		}
		preserve[i][pl] = -1;
		printf("%d\n", pl);

		for (int j = 0; j < 26; j++) {

			if (num[j]) {
				single += double(num[j]) * double(num[j]);
			}//0 is special
		}

		single = single / (double(pl) * double(pl));
		printf("%lf\n", single);

		total += (single - stdic) * (single - stdic);

	}
	total = total / step;
	printf("%lf\n", total);
	return total < eps;

}

double testt(int p, int q, int bia) {
	int m1[26] = { 0 }, m2[26] = { 0 };
	int i, j, k;
	double res = 0;

	for (i = 0; preserve[p][i] != -1; i++) {
		m1[(preserve[p][i]+bia)%26]++;
	}

	for (j = 0; preserve[q][j] != -1; j++) {
		m2[preserve[q][j]]++;
	}

	for (k = 0; k < 26; k++) {
		res += double(m1[k]) * m2[k];
	}
	res = res / (double(i) * j);
	return res;
}


void generate(int step) {
	int stp = step;

	double dest = 1;

	for (int i = 1; i < step; i++) {
		int temp = 0;
		dest = 1;
		for (int j = 0; j < 26; j++) {
			double res = testt(0, i, j) > stdic ? testt(0, i, j) - stdic : stdic - testt(0, i, j);
			printf("%d %lf\n", j, res);
			if (res < dest) {
				dest = res;
				temp = j;
			}
		}
		biaK[i] = temp;
	}
}

int main() {
	int* bias, step = 1;

	scanf("%s", str);

	int len = strlen(str);
	while (step < len && !testk(str, step))
		step++;

	//printf("%d\n", step);
	generate(step);

	for (int i = 0; i < step; i++)
		printf("%d ", biaK[i]);
	printf("\n");

	for (int i = 0; i < 26; i++) {
		for (int pos = 0; pos < len; pos++) {
			char c = 'A' + ((str[pos] - 'A' - (i + biaK[pos%step])+52) % 26);

			printf("%c", c);
		}
		printf("\n\n");
	}

}

*/

/*
#include"iostream"
#include"algorithm"

using namespace std;

struct player {
	int score;
	int num;
	int rank;

	bool operator<(player& P) {
		if (score == P.score)
			return num < P.num;
		return score > P.score;
	}
}PL[1<<8];


int main() {
	int n, p, i, j, k;
	
	while (scanf("%d", &n)) {
		if (!n)return 0;
		for (i = 0; i < n; i++) {
			scanf("%d", PL + i);
			PL[i].num = i + 1;
		}

		sort(PL, PL + n);

		PL[0].rank = 1;
		for (i = 1; i < n; i++) {
			if (PL[i].score == PL[i - 1].score)PL[i].rank = PL[i - 1].rank;
			else PL[i].rank = i + 1;
		}

		scanf("%d", &p);
		
		if (PL[p - 1].rank != p)printf("0\n");
		else {
			for (i = p - 1; i < n && PL[i].rank == p; i++)
				printf("%d ", PL[i].num);

			printf("\n");
		}
	}
}
*/


/*
#include"stdio.h"
#include"math.h"

bool check(int a, int b) {
	int mi = a < b ? a : b;
	int ma = a < b ? b : a;
	int k = ma-mi;
	double fs = (sqrt(5) + 1) / 2;
	return mi == int(fs * k);
}


int main() {
	int a, b;
	while (scanf("%d%d", &a, &b) != EOF) {
		printf("%d\n", !check(a, b));
	}
}
*/





/*
#include<string>
#include<iostream>
#include<cstring>
#include<vector>
#include<queue>
#include<stack>
using namespace std;

class Solution {
public:
	int kthSmallest(vector<vector<int>>& mat, int k) {
		int m, n, nxt;
		m = mat.size();
		n = mat[0].size();
		nxt = 0;
		priority_queue<int> q;
		vector<int> temp;
		stack<vector<int>> pre;
		temp.push_back(0);
		pre.push(temp);

		while (pre.top().size() < k) {
			int tpsz = pre.top().size();
			temp.clear();
			for (int i = 0; i < tpsz; i++)
				for (int j = 0; j < n; j++)
					temp.push_back(pre.top()[i] + mat[nxt][j]);
			
			pre.push(temp);
			nxt++;
		}
		int sz = temp.size();
		for (int i = 0; i < sz; i++) {
			q.push(temp[i]);
		}
		for (int i = 0; i < sz; i++) {
			temp[sz - 1 - i] = q.top();
			q.pop();
		}

		while (nxt < m) {
			int sz = temp.size();
			for (int i = 0; i < sz; i++) {
				q.push(temp[i] + mat[nxt][0]);
			}
			for (int i = 0; i < sz; i++) {
				for (int j = 1; j < n; j++) {
					if (temp[i] + mat[nxt][j] < q.top()) {
						q.pop();
						q.push(temp[i] + mat[nxt][j]);
					}
				}
			}

			for (int i = 0; i < sz; i++) {
				temp[sz - 1 - i] = q.top();
				q.pop();
			}
			nxt++;
		}
		return temp[k - 1];


	}
};
*/

/*
#include<iostream>
#include<fstream>
#include<map>
#include<stack>
#include<algorithm>
#include<string>

using namespace std;



struct bdata {
	int timestamp;
	int quantity;
	int price;
	string symbol;
	void parse(string s);
};

class Tp {
public:
	int maxtimegap;
	int volume;
	int wap;	// average price
	int maxprice;

	int sum;
	int lasttimestamp;
	string symbol;

	Tp() {}
	Tp(bdata& b) {
		maxtimegap = 0;
		volume = b.quantity;
		maxprice = b.price;
		sum = b.price * b.quantity;
		lasttimestamp = b.timestamp;
		symbol = b.symbol;
		wap = sum / volume;
		
	}

	Tp& operator=(const Tp& t) {
		maxtimegap = t.maxtimegap;
		volume = t.volume;
		maxprice = t.maxprice;
		wap = t.wap;
		sum = t.sum;
		lasttimestamp = t.lasttimestamp;

		symbol = t.symbol;
		return *this;
	}

	int proc(bdata& b) {
		if (symbol != b.symbol)
			return -1;

		int tgp = b.timestamp - lasttimestamp;
		if (tgp > maxtimegap)
			maxtimegap = tgp;

		volume += b.quantity;

		if (b.price > maxprice)
			maxprice = b.price;

		sum += b.price * b.quantity;
		lasttimestamp = b.timestamp;

		wap = sum / volume;
	}
};

void bdata::parse(string s) {
	int p = 0;

	timestamp = 0;
	quantity = 0;
	price = 0;
	symbol.clear();

	for (int i = 0; i < s.size(); i++) {
		if (s[i] == ',') {
			p++;
		}
		else if (p == 0) {
			timestamp *= 10;
			timestamp += s[i] - '0';
		}
		else if (p == 1) {
			symbol.push_back(s[i]);
		}
		else if (p == 2) {
			quantity *= 10;
			quantity += s[i] - '0';
		}
		else {
			price *= 10;
			price += s[i] - '0';
		}
	}
}

map<string, Tp> f;
Tp tp;

int main() {
	ifstream fp("input.csv");
	string ts;
	bdata tb;

	

	while (getline(fp, ts)) {
		
		tb.parse(ts);
		
		map<string, Tp>::iterator it = f.find(tb.symbol);
		if (it == f.end()) {
			f.insert(make_pair(tb.symbol, Tp(tb)));
		}
		else {
			it->second.proc(tb);
		}
	}

	

	ofstream ofp("output.csv");

	
	for (map<string, Tp>::iterator it = f.begin(); it != f.end(); it++) {
		ofp << it->second.symbol << "," << it->second.maxtimegap << ","
			<< it->second.volume << "," << it->second.wap << ","
			<< it->second.maxprice << endl;
	}
	

	fp.close();
	ofp.close();
	return 0;

}
*/


/*
#include<iostream>
#include<fstream>
#include<map>
#include<stack>
#include<algorithm>
#include<string>

using namespace std;



struct bdata {
	int timestamp;
	int quantity;
	int price;
	string symbol;
	void parse(string s);
};

class Tp {
public:
	int maxtimegap;
	int volume;
	int wap;	// average price
	int maxprice;

	int sum;
	int lasttimestamp;
	string symbol;

	Tp() {}
	Tp(bdata& b) {
		maxtimegap = 0;
		volume = b.quantity;
		maxprice = b.price;
		sum = b.price * b.quantity;
		lasttimestamp = b.timestamp;
		symbol = b.symbol;
		wap = sum / volume;

	}

	Tp& operator=(const Tp& t) {
		maxtimegap = t.maxtimegap;
		volume = t.volume;
		maxprice = t.maxprice;
		wap = t.wap;
		sum = t.sum;
		lasttimestamp = t.lasttimestamp;

		symbol = t.symbol;
		return *this;
	}

	int proc(bdata& b) {
		if (symbol != b.symbol)
			return -1;

		int tgp = b.timestamp - lasttimestamp;
		if (tgp > maxtimegap)
			maxtimegap = tgp;

		volume += b.quantity;

		if (b.price > maxprice)
			maxprice = b.price;

		sum += b.price * b.quantity;
		lasttimestamp = b.timestamp;

		wap = sum / volume;
	}
};

void bdata::parse(string s) {
	int p = 0;

	timestamp = 0;
	quantity = 0;
	price = 0;
	symbol.clear();

	for (int i = 0; i < s.size(); i++) {
		if (s[i] == ',') {
			p++;
		}
		else if (p == 0) {
			timestamp *= 10;
			timestamp += s[i] - '0';
		}
		else if (p == 1) {
			symbol.push_back(s[i]);
		}
		else if (p == 2) {
			quantity *= 10;
			quantity += s[i] - '0';
		}
		else {
			price *= 10;
			price += s[i] - '0';
		}
	}
}

map<string, Tp> f;
Tp tp;

int main(int argc, char** argv) {
	// that's useless
	char *cn1 = "input.csv", *cn2 = "output.csv";
	char *inp, *otp;
	if (argc > 1) {
		inp = argv[1];
		otp = argv[2];
	}
	else {
		inp = cn1;
		otp = cn2;
	}

	ifstream fp(inp);
	string ts;
	bdata tb;



	while (getline(fp, ts)) {

		tb.parse(ts);

		map<string, Tp>::iterator it = f.find(tb.symbol);
		if (it == f.end()) {
			f.insert(make_pair(tb.symbol, Tp(tb)));
		}
		else {
			it->second.proc(tb);
		}
	}



	ofstream ofp(otp);


	for (map<string, Tp>::iterator it = f.begin(); it != f.end(); it++) {
		ofp << it->second.symbol << "," << it->second.maxtimegap << ","
			<< it->second.volume << "," << it->second.wap << ","
			<< it->second.maxprice << endl;
	}


	fp.close();
	ofp.close();
	return 0;

}
*/

/**
 * Definition for a binary tree node.
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode() : val(0), left(nullptr), right(nullptr) {}
 *     TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
 *     TreeNode(int x, TreeNode *left, TreeNode *right)
 *       : val(x), left(left), right(right) {}
 * };
 **/


/**
 * 问题本质为在有向图中寻找负权环
 **/

/*
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
			res+=dfs(dest);
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
*/

/*
#include<stdio.h>
#include<string.h>

int m, n;
char maze[100][100];
char w[105];

int dir[4][2] = { 0,-1,0,1,-1,0,1,0 };

bool dfs(int i, int j, int pos) {
	
	if (w[pos + 1] == 0)
		return 1;

	int res = 0;
	char temp = maze[i][j];
	maze[i][j] = 0;
	for (int k = 0; k < 4; k++) {
		
		int x = i + dir[k][0];
		int y = j + dir[k][1];
		
		if ((x >= 0) && (x < n) && (y >= 0) && (y < m)) {
			if (maze[x][y] == w[pos + 1]) {
				res = (res || dfs(x, y, pos + 1));
			}
		}
	}
	maze[i][j] = temp;
	return res;
}

int main() {
	scanf("%d%d", &n, &m);
	scanf("%s", w);

	for (int i = 0; i < n; i++) {
		scanf("%s", maze[i]);
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			if (maze[i][j] == w[0])
				if (dfs(i, j, 0)) {
					printf("YES\n");
					return 0;
				}
	}

	printf("NO\n");
}
*/

#include<iostream>
#include<vector>
#include<string>
#include<algorithm>

using namespace std;

struct C {
	int left;
	int right;
	bool operator<(const C& c) {
		if (left == c.left)
			return right < c.right;
		return left < c.left;
	}
};

C par(string& s) {
	C res;
	int i, j, k;

	res.left = 0;
	res.right = 0;
	for (i = 0; s[i] != ','; i++) {
		res.left = res.left * 10;
		res.left += s[i] - '0';
	}
	i++;
	for (; i < s.size(); i++) {
		res.right = res.right * 10;
		res.right += s[i] - '0';
	}

	return res;
}

int main() {
	string s;
	vector<C> v;
	while (cin >> s) {
		v.push_back(par(s));
	}

	
	sort(v.begin(),v.end());

	int l = v[0].left;
	int r = v[0].right;
	for (int i = 1; i < v.size(); i++) {
		// printf("%d,%d ", v[i].left, v[i].right);
		if (v[i].left > r) {
			printf("%d,%d ", l, r);
			l = v[i].left;
			r = v[i].right;
		}
		else {
			if (v[i].right > r)
				r = v[i].right;
		}
	}
	printf("%d,%d\n", l, r);
	return 0;
}