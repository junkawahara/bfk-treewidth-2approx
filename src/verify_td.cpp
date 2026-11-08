// SPDX-License-Identifier: MIT
// verify_td.cpp : check that a PACE .td file is a valid tree decomposition of a .gr graph
// usage: verify_td graph.gr decomposition.td      -> prints "OK width=W bags=N" or "ERROR ..."
// Both files must follow the PACE 2017 formats (https://pacechallenge.org/2017/treewidth/,
// Appendix A: .gr, Appendix B: .td); format errors are reported like invalid decompositions.
// exit code: 0 valid, 1 invalid or malformed input, 2 usage / cannot open / out of memory
#include <bits/stdc++.h>
using namespace std;

struct Invalid : runtime_error {
    using runtime_error::runtime_error;
};

static bool skip_line(const string& s) {  // comment or blank
    return s.find_first_not_of(" \t\r") == string::npos || s[0] == 'c';
}

// the remaining tokens of the line as integers; false if one of them is not an integer
static bool ints(istringstream& ss, vector<long long>& out) {
    out.clear();
    string tok;
    while (ss >> tok) {
        errno = 0;
        char* end;
        long long v = strtoll(tok.c_str(), &end, 10);
        if (errno || *end) return false;
        out.push_back(v);
    }
    return true;
}

static int verify(istream& fg, istream& ft) {
    // ---- graph: "p tw n m" once, before m edge lines "u v" (loops / multiple edges allowed)
    long long n = -1, m = -1, m_read = 0;
    vector<pair<int, int>> edges;
    string line;
    int ln = 0;
    auto bad_g = [&](const string& w) { throw Invalid("graph line " + to_string(ln) + ": " + w); };
    vector<long long> x;
    while (getline(fg, line)) {
        ++ln;
        if (skip_line(line)) continue;
        istringstream ss(line);
        if (line[0] == 'p') {
            string p, tw;
            ss >> p >> tw;
            if (n >= 0) bad_g("second p-line");
            if (p != "p" || tw != "tw" || !ints(ss, x) || x.size() != 2) bad_g("malformed p-line");
            n = x[0], m = x[1];
            if (n < 0 || n > INT_MAX || m < 0) bad_g("vertex or edge count out of range");
            continue;
        }
        if (n < 0) bad_g("edge line before p-line");
        if (!ints(ss, x) || x.size() != 2) bad_g("malformed edge line");
        if (x[0] < 1 || x[0] > n || x[1] < 1 || x[1] > n) bad_g("vertex id out of range");
        ++m_read;
        if (x[0] != x[1]) edges.push_back({(int)x[0] - 1, (int)x[1] - 1});
    }
    if (n < 0) throw Invalid("graph: no p-line");
    if (m_read != m)
        throw Invalid("graph: p-line declares " + to_string(m) + " edges, found " + to_string(m_read));

    // ---- decomposition: "s td N W n" first, then exactly one "b i ..." per bag and tree edges "a b", a < b
    long long nb = -1, mb = -1, tn = -1;
    vector<vector<int>> bags;
    vector<char> defined;
    vector<vector<int>> tadj;
    long long tedges = 0;
    ln = 0;
    auto bad_t = [&](const string& w) { throw Invalid("td line " + to_string(ln) + ": " + w); };
    while (getline(ft, line)) {
        ++ln;
        if (skip_line(line)) continue;
        istringstream ss(line);
        if (nb < 0) {
            string s, td;
            ss >> s >> td;
            if (s != "s" || td != "td" || !ints(ss, x) || x.size() != 3)
                bad_t("the first line must be \"s td N W n\"");
            nb = x[0], mb = x[1], tn = x[2];
            if (nb < 0 || nb > INT_MAX || mb < 0) bad_t("bag count or bag size out of range");
            if (tn != n) bad_t("vertex count " + to_string(tn) + " != " + to_string(n));
            bags.assign(nb, {});
            defined.assign(nb, 0);
            tadj.assign(nb, {});
            continue;
        }
        if (line[0] == 's') bad_t("second s-line");
        if (line[0] == 'b') {
            string b;
            ss >> b;
            if (b != "b" || !ints(ss, x) || x.empty()) bad_t("malformed b-line");
            if (x[0] < 1 || x[0] > nb) bad_t("bag id out of range");
            int id = (int)x[0] - 1;
            if (defined[id]) bad_t("bag " + to_string(id + 1) + " defined twice");
            defined[id] = 1;
            for (size_t i = 1; i < x.size(); ++i) {
                if (x[i] < 1 || x[i] > n) bad_t("vertex out of range");
                bags[id].push_back((int)x[i] - 1);
            }
            continue;
        }
        if (!ints(ss, x) || x.size() != 2) bad_t("malformed tree edge line");
        if (x[0] < 1 || x[1] > nb || x[0] >= x[1])
            bad_t("tree edge must be \"a b\" with 1 <= a < b <= N");
        tadj[x[0] - 1].push_back((int)x[1] - 1);
        tadj[x[1] - 1].push_back((int)x[0] - 1);
        tedges++;
    }
    if (nb < 0) throw Invalid("td: no s-line");
    for (int i = 0; i < nb; ++i)
        if (!defined[i]) throw Invalid("td: bag " + to_string(i + 1) + " not defined");

    const int N = (int)nb, nv = (int)n;
    // ---- tree check
    if (N > 0) {
        if (tedges != N - 1)
            throw Invalid("not a tree (" + to_string(tedges) + " edges, " + to_string(N) + " bags)");
        vector<char> seen(N, 0);
        vector<int> st{0};
        seen[0] = 1;
        int cnt = 1;
        while (!st.empty()) {
            int x = st.back();
            st.pop_back();
            for (int y : tadj[x])
                if (!seen[y]) seen[y] = 1, cnt++, st.push_back(y);
        }
        if (cnt != N) throw Invalid("tree not connected");
    } else if (nv > 0)
        throw Invalid("no bags");
    // ---- occurrences, width and the declared maximum bag size
    vector<vector<int>> occ(nv);
    int width = -1;
    for (int i = 0; i < N; ++i) {
        auto& b = bags[i];
        sort(b.begin(), b.end());
        if (unique(b.begin(), b.end()) != b.end())
            throw Invalid("duplicate vertex in bag " + to_string(i + 1));
        width = max(width, (int)b.size() - 1);
        for (int v : b) occ[v].push_back(i);
    }
    if (width + 1 != mb)
        throw Invalid("s-line declares maximum bag size " + to_string(mb) + ", but it is " +
                      to_string(width + 1));
    for (int v = 0; v < nv; ++v)
        if (occ[v].empty()) throw Invalid("vertex " + to_string(v + 1) + " not covered");
    for (auto [u, v] : edges) {
        bool ok = false;
        const auto& a = occ[u].size() < occ[v].size() ? occ[u] : occ[v];
        int o = occ[u].size() < occ[v].size() ? v : u;
        for (int i : a)
            if (binary_search(bags[i].begin(), bags[i].end(), o)) {
                ok = true;
                break;
            }
        if (!ok) throw Invalid("edge " + to_string(u + 1) + "-" + to_string(v + 1) + " not covered");
    }
    // ---- connectivity of occurrences
    vector<int> mark(N, -1);
    for (int v = 0; v < nv; ++v) {
        for (int i : occ[v]) mark[i] = v;
        vector<int> st{occ[v][0]};
        int cnt = 0;
        unordered_set<int> seen{occ[v][0]};
        while (!st.empty()) {
            int x = st.back();
            st.pop_back();
            cnt++;
            for (int y : tadj[x])
                if (mark[y] == v && !seen.count(y)) seen.insert(y), st.push_back(y);
        }
        if (cnt != (int)occ[v].size())
            throw Invalid("occurrences of vertex " + to_string(v + 1) + " not connected");
    }
    printf("OK width=%d bags=%d\n", width, N);
    return 0;
}

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
    try {
        return verify(fg, ft);
    } catch (Invalid& e) {
        printf("ERROR %s\n", e.what());
        return 1;
    } catch (exception& e) {  // e.g. bad_alloc for a huge declared bag count
        fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
}
