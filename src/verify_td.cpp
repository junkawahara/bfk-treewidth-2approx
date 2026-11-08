// SPDX-License-Identifier: MIT
// verify_td.cpp : check that a PACE .td file is a valid tree decomposition of a .gr graph
// usage: verify_td graph.gr decomposition.td      -> prints "OK width=W" or the error
#include <bits/stdc++.h>
using namespace std;

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: verify_td graph.gr decomposition.td\n");
        return 2;
    }
    ifstream fg(argv[1]), ft(argv[2]);
    if (!fg || !ft) {
        fprintf(stderr, "cannot open input\n");
        return 2;
    }
    int n = -1;
    vector<pair<int, int>> edges;
    string line;
    while (getline(fg, line)) {
        if (line.empty() || line[0] == 'c') continue;
        istringstream ss(line);
        if (line[0] == 'p') {
            string a, b;
            long long m;
            ss >> a >> b >> n >> m;
            continue;
        }
        int u, v;
        if (ss >> u >> v)
            if (u != v) edges.push_back({u - 1, v - 1});
    }
    int nb = -1, tn = -1;
    vector<vector<int>> bags;
    vector<vector<int>> tadj;
    long long tedges = 0;
    while (getline(ft, line)) {
        if (line.empty() || line[0] == 'c') continue;
        istringstream ss(line);
        if (line[0] == 's') {
            string s, td;
            int mb;
            ss >> s >> td >> nb >> mb >> tn;
            bags.assign(nb, {});
            tadj.assign(nb, {});
            continue;
        }
        if (line[0] == 'b') {
            string b;
            int id, v;
            ss >> b >> id;
            if (id < 1 || id > nb) {
                printf("ERROR bag id\n");
                return 1;
            }
            while (ss >> v) bags[id - 1].push_back(v - 1);
            continue;
        }
        int a, b;
        if (ss >> a >> b) {
            if (a < 1 || b < 1 || a > nb || b > nb) {
                printf("ERROR edge id\n");
                return 1;
            }
            tadj[a - 1].push_back(b - 1);
            tadj[b - 1].push_back(a - 1);
            tedges++;
        }
    }
    if (tn != n) {
        printf("ERROR vertex count %d != %d\n", tn, n);
        return 1;
    }
    // tree check
    if (nb > 0) {
        if (tedges != nb - 1) {
            printf("ERROR not a tree (%lld edges, %d bags)\n", tedges, nb);
            return 1;
        }
        vector<char> seen(nb, 0);
        vector<int> st{0};
        seen[0] = 1;
        int cnt = 1;
        while (!st.empty()) {
            int x = st.back();
            st.pop_back();
            for (int y : tadj[x])
                if (!seen[y]) seen[y] = 1, cnt++, st.push_back(y);
        }
        if (cnt != nb) {
            printf("ERROR tree not connected\n");
            return 1;
        }
    } else if (n > 0) {
        printf("ERROR no bags\n");
        return 1;
    }
    // occurrences
    vector<vector<int>> occ(n);
    int width = -1;
    for (int i = 0; i < nb; ++i) {
        auto b = bags[i];
        sort(b.begin(), b.end());
        if (unique(b.begin(), b.end()) != b.end()) {
            printf("ERROR duplicate vertex in bag %d\n", i + 1);
            return 1;
        }
        width = max(width, (int)b.size() - 1);
        for (int v : b) {
            if (v < 0 || v >= n) {
                printf("ERROR vertex out of range\n");
                return 1;
            }
            occ[v].push_back(i);
        }
        bags[i] = b;
    }
    for (int v = 0; v < n; ++v)
        if (occ[v].empty()) {
            printf("ERROR vertex %d not covered\n", v + 1);
            return 1;
        }
    for (auto [u, v] : edges) {
        bool ok = false;
        const auto& a = occ[u].size() < occ[v].size() ? occ[u] : occ[v];
        int o = occ[u].size() < occ[v].size() ? v : u;
        for (int i : a)
            if (binary_search(bags[i].begin(), bags[i].end(), o)) {
                ok = true;
                break;
            }
        if (!ok) {
            printf("ERROR edge %d-%d not covered\n", u + 1, v + 1);
            return 1;
        }
    }
    // connectivity of occurrences
    vector<int> mark(nb, -1);
    for (int v = 0; v < n; ++v) {
        for (int i : occ[v]) mark[i] = v;
        vector<int> st{occ[v][0]};
        vector<char> vis;
        int cnt = 0;
        unordered_set<int> seen{occ[v][0]};
        while (!st.empty()) {
            int x = st.back();
            st.pop_back();
            cnt++;
            for (int y : tadj[x])
                if (mark[y] == v && !seen.count(y)) seen.insert(y), st.push_back(y);
        }
        if (cnt != (int)occ[v].size()) {
            printf("ERROR occurrences of vertex %d not connected\n", v + 1);
            return 1;
        }
    }
    printf("OK width=%d bags=%d\n", width, nb);
    return 0;
}
