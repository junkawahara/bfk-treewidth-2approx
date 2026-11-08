// gen.cpp : test graph generator (PACE .gr on stdout)
//   gen ktree N K P SEED   random partial k-tree (tw <= K), each k-tree edge kept with prob. P
//   gen grid R C           R x C grid (tw = min(R,C))
//   gen gnm N M SEED       uniform random graph
#include <bits/stdc++.h>
using namespace std;

static void out(int n, set<pair<int, int>>& E, const char* comment) {
    printf("c %s\n", comment);
    printf("p tw %d %zu\n", n, E.size());
    for (auto [u, v] : E) printf("%d %d\n", u + 1, v + 1);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: gen ktree N K P SEED | grid R C | gnm N M SEED\n");
        return 2;
    }
    string t = argv[1];
    set<pair<int, int>> E;
    auto add = [&](int u, int v) {
        if (u == v) return;
        if (u > v) swap(u, v);
        E.insert({u, v});
    };
    char buf[256];
    if (t == "ktree" && argc == 6) {
        int n = atoi(argv[2]), k = atoi(argv[3]);
        double p = atof(argv[4]);
        mt19937_64 rng(atoll(argv[5]));
        uniform_real_distribution<double> U(0, 1);
        vector<int> perm(n);
        iota(perm.begin(), perm.end(), 0);
        shuffle(perm.begin(), perm.end(), rng);
        vector<vector<int>> cliques;  // k-cliques
        int k0 = min(n, k + 1);
        vector<pair<int, int>> kedges;
        for (int i = 0; i < k0; ++i)
            for (int j = i + 1; j < k0; ++j) kedges.push_back({i, j});
        if (n > k) {
            for (int drop = 0; drop < k0; ++drop) {
                vector<int> c;
                for (int i = 0; i < k0; ++i)
                    if (i != drop) c.push_back(i);
                cliques.push_back(c);
            }
        }
        for (int v = k0; v < n; ++v) {
            vector<int> c = cliques[rng() % cliques.size()];
            for (int u : c) kedges.push_back({u, v});
            for (int drop = 0; drop < k; ++drop) {
                vector<int> d = c;
                d[drop] = v;
                cliques.push_back(d);
            }
        }
        for (auto [u, v] : kedges)
            if (U(rng) < p) add(perm[u], perm[v]);
        snprintf(buf, sizeof buf, "partial %d-tree n=%d p=%.2f seed=%s", k, n, p, argv[5]);
        out(n, E, buf);
    } else if (t == "grid" && argc == 4) {
        int R = atoi(argv[2]), C = atoi(argv[3]);
        for (int i = 0; i < R; ++i)
            for (int j = 0; j < C; ++j) {
                int v = i * C + j;
                if (i + 1 < R) add(v, v + C);
                if (j + 1 < C) add(v, v + 1);
            }
        snprintf(buf, sizeof buf, "grid %dx%d", R, C);
        out(R * C, E, buf);
    } else if (t == "gnm" && argc == 5) {
        int n = atoi(argv[2]);
        long long m = atoll(argv[3]);
        mt19937_64 rng(atoll(argv[4]));
        long long maxm = (long long)n * (n - 1) / 2;
        m = min(m, maxm);
        while ((long long)E.size() < m) add(rng() % n, rng() % n);
        snprintf(buf, sizeof buf, "gnm n=%d m=%lld seed=%s", n, m, argv[4]);
        out(n, E, buf);
    } else {
        fprintf(stderr, "bad arguments\n");
        return 2;
    }
    return 0;
}
