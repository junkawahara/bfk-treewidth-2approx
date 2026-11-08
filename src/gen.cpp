// SPDX-License-Identifier: MIT
// gen.cpp : test graph generator (PACE .gr on stdout)
//   gen ktree N K P SEED   random partial k-tree (tw <= K), each k-tree edge kept with prob. P
//   gen grid R C           R x C grid (tw = min(R,C), except tw = 0 for the 1 x 1 grid)
//   gen gnm N M SEED       uniform random graph
#include <bits/stdc++.h>
using namespace std;

// strict argument parsing: the whole string, in [lo, hi] (and finite for reals); otherwise exit 2
static long long arg_int(const char* s, long long lo, long long hi, const char* name) {
    errno = 0;
    char* end;
    long long v = strtoll(s, &end, 10);
    if (!*s || isspace((unsigned char)*s) || errno || *end || v < lo || v > hi) {
        fprintf(stderr, "error: %s must be an integer in [%lld, %lld], got '%s'\n", name, lo, hi, s);
        exit(2);
    }
    return v;
}
static double arg_real(const char* s, double lo, double hi, const char* name) {
    errno = 0;
    char* end;
    double v = strtod(s, &end);
    if (!*s || isspace((unsigned char)*s) || errno || *end || !isfinite(v) || v < lo || v > hi) {
        fprintf(stderr, "error: %s must be a number in [%g, %g], got '%s'\n", name, lo, hi, s);
        exit(2);
    }
    return v;
}

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
        int n = (int)arg_int(argv[2], 0, INT_MAX, "N");
        int k = (int)arg_int(argv[3], 0, INT_MAX - 1, "K");
        double p = arg_real(argv[4], 0, 1, "P");
        mt19937_64 rng(arg_int(argv[5], 0, LLONG_MAX, "SEED"));
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
        int R = (int)arg_int(argv[2], 1, INT_MAX, "R"), C = (int)arg_int(argv[3], 1, INT_MAX, "C");
        if ((long long)R * C > INT_MAX) {
            fprintf(stderr, "error: R*C must be at most %d\n", INT_MAX);
            return 2;
        }
        for (int i = 0; i < R; ++i)
            for (int j = 0; j < C; ++j) {
                int v = i * C + j;
                if (i + 1 < R) add(v, v + C);
                if (j + 1 < C) add(v, v + 1);
            }
        snprintf(buf, sizeof buf, "grid %dx%d", R, C);
        out(R * C, E, buf);
    } else if (t == "gnm" && argc == 5) {
        int n = (int)arg_int(argv[2], 0, INT_MAX, "N");
        long long m = arg_int(argv[3], 0, LLONG_MAX, "M");
        mt19937_64 rng(arg_int(argv[4], 0, LLONG_MAX, "SEED"));
        long long maxm = (long long)n * (n - 1) / 2;
        m = min(m, maxm);
        while ((long long)E.size() < m) add(rng() % n, rng() % n);
        snprintf(buf, sizeof buf, "gnm n=%d m=%lld seed=%s", n, m, argv[4]);
        out(n, E, buf);
    } else {
        fprintf(stderr, "bad arguments\n");
        return 2;
    }
    if (fflush(stdout) != 0 || ferror(stdout)) {  // e.g. a full disk: do not hand on a truncated graph
        fprintf(stderr, "error: writing the graph failed\n");
        return 1;
    }
    return 0;
}
