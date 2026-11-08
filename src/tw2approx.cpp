// SPDX-License-Identifier: MIT
// tw2approx.cpp
//
// Implementation of
//   M. Belbasi, M. Fürer, M. Kumar,
//   "Optimized 2-Approximation of Treewidth", arXiv:2411.16918 (v3).
//
// Given a graph G (PACE .gr format) the program computes a tree decomposition
// of width <= 2k+1 or concludes tw(G) > k.  In the default (auto) mode k is
// searched upward, so the output width is <= 2*tw(G)+1.
//
// Main ingredients (section numbers refer to the paper):
//   * grouped tree decomposition (Def. 6, Lemma 7) with main / intersection nodes
//   * DP tables over partitions (W1,W2,W3,X) of each bag (Sec. 2.2, Def. 14)
//     stored at every node and maintained incrementally by Addition /
//     Subtraction / Update (Sec. 3.1)
//   * rounds w = width .. 2k+2, a DFS instance per round (Sec. 3.2)
//   * Move of a maximum bag to the root by edge rotations (Alg. 2, 3)
//   * Merge of nodes that lost their home vertex (Alg. 4)
//   * Split (Def. 11, Alg. 5, 6): choose a good partition at the root,
//     choose the split top-down only inside the editable part, create the
//     copies x_i, re-use the DP table of x for the non-editable children,
//     create the new root with bag S and the forget-chains to the copies r_i.
//   * 3-way partitions (C1,C2,S) while |B_r| >= 3k+4 and 4-way partitions
//     (C1,C2,C3,S) afterwards (Remark 24).
//
// See README.md for the places where the paper is under-specified and the
// choices made here.

#include <bits/stdc++.h>
using namespace std;

typedef int32_t val_t;
static const val_t INF = 1000000000;
static int g_verbose = 0;

// ---------------------------------------------------------------- graph
struct Graph {
    int n = 0;
    long long m = 0;
    vector<vector<int>> adj;  // sorted, no self loops, no duplicates
    bool has_edge(int u, int v) const {
        const auto& a = adj[u];
        return binary_search(a.begin(), a.end(), v);
    }
};

static Graph read_gr(istream& in) {
    Graph G;
    string line;
    bool header = false;
    while (getline(in, line)) {
        if (line.empty() || line[0] == 'c') continue;
        istringstream ss(line);
        if (line[0] == 'p') {
            string p, tw;
            long long n, m;
            ss >> p >> tw >> n >> m;
            G.n = (int)n;
            G.adj.assign(G.n, {});
            header = true;
            continue;
        }
        long long u, v;
        if (!(ss >> u >> v)) continue;
        if (!header) throw runtime_error("edge before p-line");
        if (u < 1 || v < 1 || u > G.n || v > G.n) throw runtime_error("vertex id out of range");
        --u;
        --v;
        if (u == v) continue;
        G.adj[u].push_back((int)v);
        G.adj[v].push_back((int)u);
    }
    if (!header) throw runtime_error("no p-line found");
    for (auto& a : G.adj) {
        sort(a.begin(), a.end());
        a.erase(unique(a.begin(), a.end()), a.end());
        G.m += (long long)a.size();
    }
    G.m /= 2;
    return G;
}

// ------------------------------------------------------- sorted vectors
static vector<int> v_union(const vector<int>& a, const vector<int>& b) {
    vector<int> r;
    r.reserve(a.size() + b.size());
    set_union(a.begin(), a.end(), b.begin(), b.end(), back_inserter(r));
    return r;
}
static vector<int> v_inter(const vector<int>& a, const vector<int>& b) {
    vector<int> r;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(r));
    return r;
}
static vector<int> v_diff(const vector<int>& a, const vector<int>& b) {
    vector<int> r;
    set_difference(a.begin(), a.end(), b.begin(), b.end(), back_inserter(r));
    return r;
}
static bool v_subset(const vector<int>& a, const vector<int>& b) {  // a subset of b
    return includes(b.begin(), b.end(), a.begin(), a.end());
}
static uint64_t bag_hash(const vector<int>& b) {
    uint64_t h = 1469598103934665603ULL ^ b.size();
    for (int v : b) {
        h ^= (uint64_t)v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h *= 1099511628211ULL;
    }
    return h;
}

// ------------------------------------------- initial tree decomposition
// Elimination ordering -> tree decomposition.  Node of v has bag {v} u N+(v)
// and parent = earliest eliminated vertex of N+(v).  This decomposition has
// the unique home property (Def. 4) for free: node v is the home of exactly v.
struct Elim {
    vector<int> order, pos, parent;
    vector<vector<int>> bag;
    int width = -1;
};

// mode: "mindeg", "random", or "mixed" (min-degree, but with probability mix a
// uniformly random remaining vertex is eliminated -> deliberately worse decompositions)
static Elim eliminate(const Graph& G, const string& mode, uint64_t seed, double mix = 0) {
    int n = G.n;
    Elim e;
    e.pos.assign(n, -1);
    e.parent.assign(n, -1);
    e.bag.assign(n, {});
    vector<vector<int>> a = G.adj;
    vector<int> perm;
    bool rnd = (mode == "random"), mixed = (mode == "mixed");
    mt19937_64 rng(seed);
    if (rnd) {
        perm.resize(n);
        iota(perm.begin(), perm.end(), 0);
        shuffle(perm.begin(), perm.end(), rng);
    }
    vector<int> remain, rpos;  // for "mixed": remaining vertices (swap-remove)
    if (mixed) {
        remain.resize(n);
        iota(remain.begin(), remain.end(), 0);
        rpos = remain;
    }
    uniform_real_distribution<double> U(0, 1);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    if (!rnd)
        for (int v = 0; v < n; ++v) pq.push({(int)a[v].size(), v});
    for (int step = 0; step < n; ++step) {
        int v;
        if (rnd) {
            v = perm[step];
        } else if (mixed && U(rng) < mix) {
            v = remain[rng() % remain.size()];
        } else {
            while (true) {
                auto [d, u] = pq.top();
                pq.pop();
                if (e.pos[u] >= 0 || d != (int)a[u].size()) continue;
                v = u;
                break;
            }
        }
        const vector<int> N = a[v];
        vector<int> B = N;
        B.insert(lower_bound(B.begin(), B.end(), v), v);
        e.width = max(e.width, (int)B.size() - 1);
        e.bag[v] = B;
        for (int u : N) {
            vector<int> m;
            m.reserve(a[u].size() + N.size());
            size_t i = 0, j = 0;
            const auto& au = a[u];
            while (i < au.size() || j < N.size()) {
                int x;
                if (j >= N.size() || (i < au.size() && au[i] < N[j]))
                    x = au[i++];
                else if (i >= au.size() || N[j] < au[i])
                    x = N[j++];
                else {
                    x = au[i];
                    ++i;
                    ++j;
                }
                if (x != u && x != v) m.push_back(x);
            }
            a[u].swap(m);
            if (!rnd) pq.push({(int)a[u].size(), u});
        }
        vector<int>().swap(a[v]);
        if (mixed) {
            int last = remain.back();
            remain[rpos[v]] = last;
            rpos[last] = rpos[v];
            remain.pop_back();
        }
        e.pos[v] = step;
        e.order.push_back(v);
    }
    for (int v = 0; v < n; ++v) {
        int best = -1;
        for (int u : e.bag[v])
            if (u != v && (best < 0 || e.pos[u] < e.pos[best])) best = u;
        e.parent[v] = best;
    }
    return e;
}

// core numbers (degeneracy) -- a cheap lower bound: tw >= degeneracy
static vector<int> core_numbers(const Graph& G) {
    int n = G.n, md = 0;
    vector<int> deg(n), core(n, 0);
    for (int v = 0; v < n; ++v) deg[v] = (int)G.adj[v].size(), md = max(md, deg[v]);
    vector<vector<int>> bucket(md + 1);
    for (int v = 0; v < n; ++v) bucket[deg[v]].push_back(v);
    vector<char> done(n, 0);
    int cur = 0;
    for (int d = 0; d <= md; ++d) {
        while (!bucket[d].empty()) {
            int v = bucket[d].back();
            bucket[d].pop_back();
            if (done[v] || deg[v] != d) continue;
            done[v] = 1;
            cur = max(cur, d);
            core[v] = cur;
            for (int u : G.adj[v])
                if (!done[u] && deg[u] > d) {
                    --deg[u];
                    bucket[deg[u]].push_back(u);
                }
        }
    }
    return core;
}

// ------------------------------------------------- grouped decomposition
struct Node {
    vector<int> bag;  // sorted
    vector<int> ch;   // children
    vector<val_t> E;  // DP table, E[P] = s(A[P]) - |X_P|  (INF = illegal)
    int parent = -1, pos = -1;
    bool main = true, alive = true;
    uint8_t st = 0;  // DFS status: 0 unvisited, 1 pre-visited, 2 post-visited
    uint64_t bh = 0;
};

struct Stats {
    long long table_ops = 0, entries = 0, rotations = 0, splits = 0, merges = 0, rounds = 0,
              recomputes = 0, size_violations = 0, dfs_restarts = 0;
};

struct GTD {
    const Graph& G;
    int k = 0;
    int L = 4;  // number of labels: 3 -> (C1,C2,S), 4 -> (C1,C2,C3,S)
    size_t Lp[40];
    deque<Node> nd;
    vector<int> freel;
    int root = -1;
    bool tables_valid = false;
    size_t mem_entries = 0, peak_entries = 0, mem_limit_entries = 0;
    bool check = false;
    Stats st;
    // scratch for split
    vector<int8_t> vlab;
    vector<int> vstamp;
    int stamp = 0;

    explicit GTD(const Graph& g, size_t lim) : G(g), mem_limit_entries(lim) {
        vlab.assign(G.n, 0);
        vstamp.assign(G.n, 0);
        set_L(4);
    }

    void set_L(int l) {
        L = l;
        Lp[0] = 1;
        for (int i = 1; i < 40; ++i) Lp[i] = (Lp[i - 1] > (SIZE_MAX >> 3)) ? SIZE_MAX : Lp[i - 1] * L;
    }

    // ---------------- memory accounting
    void tfree(vector<val_t>& T) {
        mem_entries -= T.size();
        vector<val_t>().swap(T);
    }
    void talloc(vector<val_t>& T, size_t sz) {
        tfree(T);
        if (sz == SIZE_MAX || mem_entries + sz > mem_limit_entries) {
            ostringstream os;
            os << "memory limit exceeded (table of " << sz << " entries, in use " << mem_entries << ")";
            throw runtime_error(os.str());
        }
        T.assign(sz, 0);
        mem_entries += sz;
        peak_entries = max(peak_entries, mem_entries);
    }

    // ---------------- node management
    int new_node(bool main, vector<int> bag) {
        int id;
        if (!freel.empty()) {
            id = freel.back();
            freel.pop_back();
            nd[id] = Node();
        } else {
            id = (int)nd.size();
            nd.emplace_back();
        }
        Node& x = nd[id];
        x.main = main;
        x.bag = std::move(bag);
        x.bh = bag_hash(x.bag);
        x.alive = true;
        return id;
    }
    void del_node(int id) {
        Node& x = nd[id];
        tfree(x.E);
        x.alive = false;
        vector<int>().swap(x.ch);
        vector<int>().swap(x.bag);
        x.parent = x.pos = -1;
        freel.push_back(id);
    }
    void add_child(int p, int c) {
        nd[c].parent = p;
        nd[c].pos = (int)nd[p].ch.size();
        nd[p].ch.push_back(c);
    }
    void remove_child(int p, int c) {
        auto& ch = nd[p].ch;
        int i = nd[c].pos;
        assert(i >= 0 && i < (int)ch.size() && ch[i] == c);
        int last = ch.back();
        ch[i] = last;
        nd[last].pos = i;
        ch.pop_back();
        nd[c].parent = -1;
        nd[c].pos = -1;
    }
    int find_inter(int m, const vector<int>& bag) {
        uint64_t h = bag_hash(bag);
        for (int c : nd[m].ch)
            if (nd[c].bh == h && nd[c].bag == bag) return c;
        return -1;
    }
    static int diffpos(const vector<int>& big, const vector<int>& small) {
        size_t q = 0;
        while (q < small.size() && big[q] == small[q]) ++q;
        return (int)q;
    }
    int home_vertex(int y) {  // non-root main node
        const auto& By = nd[y].bag;
        const auto& Bz = nd[nd[y].parent].bag;
        return By[diffpos(By, Bz)];
    }
    int width() {
        int w = -1;
        for (auto& x : nd)
            if (x.alive && x.main) w = max(w, (int)x.bag.size() - 1);
        return w;
    }

    // ---------------- table operations
    // weights of positions of B in the index space of the sub-bag S (S subset of B)
    void sub_weights(const vector<int>& B, const vector<int>& S, size_t* w) {
        size_t j = 0;
        for (size_t p = 0; p < B.size(); ++p) {
            if (j < S.size() && S[j] == B[p])
                w[p] = Lp[j++];
            else
                w[p] = 0;
        }
        if (j != S.size()) throw logic_error("sub_weights: not a subset");
    }

    // T over B := 0 for legal partitions, INF for illegal ones (Initialization, Sec. 3.1.1)
    void legal_init(vector<val_t>& T, const vector<int>& B) {
        int b = (int)B.size();
        if (b > 30) throw runtime_error("bag too large");
        talloc(T, Lp[b]);
        uint32_t nbr[32] = {0};
        for (int p = 0; p < b; ++p)
            for (int q = p + 1; q < b; ++q)
                if (G.has_edge(B[p], B[q])) nbr[p] |= 1u << q, nbr[q] |= 1u << p;
        st.table_ops++;
        st.entries += (long long)Lp[b];
        uint32_t lm[4] = {0, 0, 0, 0};
        legal_rec(T.data(), b - 1, 0, nbr, lm);
    }
    void legal_rec(val_t* T, int p, size_t base, const uint32_t* nbr, uint32_t* lm) {
        if (p < 0) {
            T[base] = 0;
            return;
        }
        uint32_t any = lm[1] | lm[2] | lm[3];
        for (int l = 0; l < L; ++l) {
            size_t nb = base + (size_t)l * Lp[p];
            if (l > 0 && (nbr[p] & any & ~lm[l])) {
                std::fill(T + nb, T + nb + Lp[p], INF);
                continue;
            }
            lm[l] |= 1u << p;
            legal_rec(T, p - 1, nb, nbr, lm);
            lm[l] &= ~(1u << p);
        }
    }

    // T(B)[P] += sign * S(BS)[P|BS] for all legal P (BS subset of B)
    // (Addition / Subtraction of an intersection child to a main node, Sec. 3.1.2 case 2)
    void proj_add(vector<val_t>& T, const vector<int>& B, const vector<val_t>& S, const vector<int>& BS,
                  int sign) {
        int b = (int)B.size();
        size_t w[32];
        sub_weights(B, BS, w);
        st.table_ops++;
        st.entries += (long long)Lp[b];
        val_t* t = T.data();
        const val_t* s = S.data();
        if (b == 0) {
            if (t[0] != INF) t[0] += sign * s[0];
            return;
        }
        int d[32] = {0};
        size_t si = 0, w0 = w[0], outer = Lp[b - 1];
        size_t idx = 0;
        for (size_t o = 0; o < outer; ++o) {
            for (int l = 0; l < L; ++l, ++idx)
                if (t[idx] != INF) t[idx] += sign * s[si + l * w0];
            for (int p = 1; p < b; ++p) {
                if (d[p] < L - 1) {
                    d[p]++;
                    si += w[p];
                    break;
                }
                si -= (size_t)(L - 1) * w[p];
                d[p] = 0;
            }
        }
    }
    void add_im(int m, int c, int sign) { proj_add(nd[m].E, nd[m].bag, nd[c].E, nd[c].bag, sign); }

    // Addition / Subtraction of a main child y to an intersection node z (Sec. 3.1.2 case 1)
    // E_z[P] += sign * min_l ( E_y[P + v:l] + [l == X] )
    void add_mi(int z, int y, int sign) {
        const auto& Bz = nd[z].bag;
        const auto& By = nd[y].bag;
        if (By.size() != Bz.size() + 1) throw logic_error("add_mi: bag sizes");
        int q = diffpos(By, Bz);
        int bz = (int)Bz.size();
        size_t lo_n = Lp[q], hi_n = Lp[bz - q], hstep = Lp[q + 1];
        val_t* Ez = nd[z].E.data();
        const val_t* Ey = nd[y].E.data();
        st.table_ops++;
        st.entries += (long long)Lp[bz] * L;
        for (size_t hi = 0; hi < hi_n; ++hi) {
            size_t iz = hi * lo_n, by = hi * hstep;
            for (size_t lo = 0; lo < lo_n; ++lo, ++iz) {
                if (Ez[iz] == INF) continue;
                size_t base = by + lo;
                val_t m = Ey[base] + 1;
                for (int l = 1; l < L; ++l) m = min(m, Ey[base + l * lo_n]);
                Ez[iz] += sign * m;
            }
        }
    }

    // T(B') := S(B)[P' with B\B' labelled X]   (B' subset of B)   -- Alg. 5
    void restrict_table(const vector<val_t>& S, const vector<int>& B, vector<val_t>& T, const vector<int>& Bp) {
        int bp = (int)Bp.size();
        talloc(T, Lp[bp]);
        size_t w[32];
        {
            size_t j = 0;
            for (int p = 0; p < bp; ++p) {
                while (B[j] != Bp[p]) ++j;
                w[p] = Lp[j];
            }
        }
        st.table_ops++;
        st.entries += (long long)Lp[bp];
        int d[32] = {0};
        size_t si = 0;
        for (size_t idx = 0; idx < Lp[bp]; ++idx) {
            T[idx] = S[si];
            for (int p = 0; p < bp; ++p) {
                if (d[p] < L - 1) {
                    d[p]++;
                    si += w[p];
                    break;
                }
                si -= (size_t)(L - 1) * w[p];
                d[p] = 0;
            }
        }
    }
    void merge_same_bag(vector<val_t>& T, const vector<val_t>& S) {
        st.table_ops++;
        st.entries += (long long)T.size();
        for (size_t i = 0; i < T.size(); ++i)
            if (T[i] != INF) T[i] += S[i];
    }

    // ---------------- full (re)computation of all tables, bottom-up
    vector<int> preorder() {
        vector<int> ord, stk{root};
        while (!stk.empty()) {
            int x = stk.back();
            stk.pop_back();
            ord.push_back(x);
            for (int c : nd[x].ch) stk.push_back(c);
        }
        return ord;
    }
    void compute_node(int id) {  // Initialization + Addition of all children
        legal_init(nd[id].E, nd[id].bag);
        if (nd[id].main) {
            for (int c : nd[id].ch) add_im(id, c, +1);
        } else {
            for (int y : nd[id].ch) add_mi(id, y, +1);
        }
    }
    size_t estimate_entries() {
        size_t tot = 0;
        for (auto& x : nd)
            if (x.alive) {
                size_t s = Lp[min<size_t>(x.bag.size(), 39)];
                if (s == SIZE_MAX || tot + s < tot) return SIZE_MAX;
                tot += s;
            }
        return tot;
    }
    void recompute_all() {
        for (auto& x : nd)
            if (x.alive) tfree(x.E);
        size_t est = estimate_entries();
        if (est == SIZE_MAX || est > mem_limit_entries) {
            ostringstream os;
            os << "memory limit: DP tables would need " << (est == SIZE_MAX ? -1.0 : est * 4.0 / 1e9)
               << " GB (L=" << L << ")";
            throw runtime_error(os.str());
        }
        vector<int> ord = preorder();
        for (int i = (int)ord.size() - 1; i >= 0; --i) compute_node(ord[i]);
        tables_valid = true;
        st.recomputes++;
    }

    // ---------------- consistency checks (--check)
    void verify_structure() {
        if (nd[root].parent != -1 || !nd[root].main || nd[root].bag.empty())
            throw logic_error("verify: bad root");
        vector<int> ord = preorder();
        for (int x : ord) {
            Node& X = nd[x];
            if (!X.alive) throw logic_error("verify: dead node in tree");
            for (size_t i = 0; i < X.ch.size(); ++i) {
                int c = X.ch[i];
                if (nd[c].parent != x || nd[c].pos != (int)i) throw logic_error("verify: parent/pos");
                if (nd[c].main == X.main) throw logic_error("verify: alternation");
            }
            if (!X.main) {
                if (X.ch.empty()) throw logic_error("verify: inter leaf");
                if (!v_subset(X.bag, nd[X.parent].bag)) throw logic_error("verify: inter not subset of parent");
                for (int y : X.ch) {
                    if (v_inter(nd[y].bag, nd[X.parent].bag) != X.bag)
                        throw logic_error("verify: inter != intersection");
                    if (nd[y].bag.size() != X.bag.size() + 1) throw logic_error("verify: unique home");
                }
            } else {
                for (size_t i = 0; i < X.ch.size(); ++i)
                    for (size_t j = i + 1; j < X.ch.size(); ++j)
                        if (nd[X.ch[i]].bag == nd[X.ch[j]].bag) throw logic_error("verify: sibling inter equal");
            }
        }
    }
    void verify_tables() {
        vector<int> ord = preorder();
        unordered_map<int, vector<val_t>> fresh;
        for (int i = (int)ord.size() - 1; i >= 0; --i) {
            int id = ord[i];
            vector<val_t> E;
            legal_init(E, nd[id].bag);
            if (nd[id].main) {
                for (int c : nd[id].ch) proj_add(E, nd[id].bag, fresh[c], nd[c].bag, +1);
            } else {
                // emulate add_mi using fresh child tables
                const auto& Bz = nd[id].bag;
                for (int y : nd[id].ch) {
                    const auto& By = nd[y].bag;
                    int q = diffpos(By, Bz);
                    size_t lo_n = Lp[q], hi_n = Lp[Bz.size() - q], hstep = Lp[q + 1];
                    const auto& Ey = fresh[y];
                    for (size_t hi = 0; hi < hi_n; ++hi)
                        for (size_t lo = 0; lo < lo_n; ++lo) {
                            size_t iz = hi * lo_n + lo;
                            if (E[iz] == INF) continue;
                            size_t base = hi * hstep + lo;
                            val_t m = Ey[base] + 1;
                            for (int l = 1; l < L; ++l) m = min(m, Ey[base + l * lo_n]);
                            E[iz] += m;
                        }
                }
            }
            if (E != nd[id].E) {
                ostringstream os;
                os << "verify: table mismatch at node " << id << (nd[id].main ? " (main)" : " (inter)");
                throw logic_error(os.str());
            }
            fresh[id] = std::move(E);
        }
        for (auto& kv : fresh) tfree(kv.second);
    }
    void do_checks(const char* where) {
        if (!check) return;
        try {
            verify_structure();
            verify_tables();
        } catch (exception& e) {
            throw logic_error(string(e.what()) + " [after " + where + "]");
        }
    }

    // ---------------- building (Lemma 7)
    void build(const vector<int>& verts, const Elim& e) {
        unordered_map<int, int> id;
        for (int v : verts) id[v] = new_node(true, e.bag[v]);
        unordered_map<int, vector<int>> kids;
        for (int v : verts) {
            if (e.parent[v] < 0)
                root = id[v];
            else
                kids[e.parent[v]].push_back(v);
        }
        for (auto& kv : kids) {
            int P = id[kv.first];
            map<vector<int>, int> groups;  // the "temporary binary tree" of Lemma 7
            for (int v : kv.second) {
                vector<int> I = e.bag[v];
                I.erase(lower_bound(I.begin(), I.end(), v));
                auto it = groups.find(I);
                int z;
                if (it == groups.end()) {
                    z = new_node(false, I);
                    add_child(P, z);
                    groups[I] = z;
                } else
                    z = it->second;
                add_child(z, id[v]);
            }
        }
    }

    // ---------------- DFS instance (Sec. 3.2)
    int dfs_find(int target) {
        struct F {
            int x, ci, yi;
        };
        vector<F> stk;
        if (nd[root].st != 2) stk.push_back({root, -1, 0});
        while (!stk.empty()) {
            F& f = stk.back();
            int x = f.x;
            if (f.ci < 0) {
                if ((int)nd[x].bag.size() == target) return x;
                nd[x].st = 1;
                f.ci = 0;
                f.yi = 0;
            }
            bool pushed = false;
            while (f.ci < (int)nd[x].ch.size()) {
                int c = nd[x].ch[f.ci];
                if (f.yi < (int)nd[c].ch.size()) {
                    int y = nd[c].ch[f.yi++];
                    if (nd[y].st != 2) {
                        stk.push_back({y, -1, 0});
                        pushed = true;
                        break;
                    }
                } else {
                    f.ci++;
                    f.yi = 0;
                }
            }
            if (!pushed) {
                nd[x].st = 2;
                stk.pop_back();
            }
        }
        return -1;
    }

    // ---------------- attach helpers (keep property 5 of Def. 6: sibling inter bags distinct)
    // attach inter node c under main m, updating E_m
    void attach_inter_table(int m, int c) {
        int c2 = find_inter(m, nd[c].bag);
        add_im(m, c, +1);
        if (c2 < 0) {
            add_child(m, c);
        } else {
            merge_same_bag(nd[c2].E, nd[c].E);
            vector<int> chs = nd[c].ch;
            for (int y : chs) add_child(c2, y);
            nd[c].ch.clear();
            del_node(c);
            st.merges++;
        }
    }
    // attach inter node c under main m, contribution already contained in E_m
    void attach_inter_struct(int m, int c) {
        int c2 = find_inter(m, nd[c].bag);
        if (c2 < 0) {
            add_child(m, c);
        } else {
            merge_same_bag(nd[c2].E, nd[c].E);
            vector<int> chs = nd[c].ch;
            for (int y : chs) add_child(c2, y);
            nd[c].ch.clear();
            del_node(c);
            st.merges++;
        }
    }
    // attach main nodes ys (all with B_y cap B_m == I, |B_y \ I| = 1) under main m
    void attach_mains(int m, const vector<int>& I, const vector<int>& ys, vector<int>* created = nullptr) {
        int c = find_inter(m, I);
        if (c < 0) {
            c = new_node(false, I);
            if (created) created->push_back(c);
            legal_init(nd[c].E, I);
            for (int y : ys) {
                add_child(c, y);
                add_mi(c, y, +1);
            }
            add_child(m, c);
            add_im(m, c, +1);
        } else {
            add_im(m, c, -1);
            for (int y : ys) {
                add_child(c, y);
                add_mi(c, y, +1);
            }
            add_im(m, c, +1);
        }
    }

    // forget chain: returns the head main node with bag I u {D[0]} whose
    // descendants lead to 'bottom' (bag I u D).  If |D| == 1 returns bottom.
    int build_chain(int bottom, const vector<int>& I, const vector<int>& D, vector<int>* created = nullptr) {
        int cur = bottom;
        for (int t = (int)D.size() - 1; t >= 1; --t) {
            vector<int> bt = I;
            for (int u = 0; u < t; ++u) bt.push_back(D[u]);
            sort(bt.begin(), bt.end());
            int ib = new_node(false, bt);
            legal_init(nd[ib].E, bt);
            add_child(ib, cur);
            add_mi(ib, cur, +1);
            int cm = new_node(true, bt);
            legal_init(nd[cm].E, bt);
            add_child(cm, ib);
            add_im(cm, ib, +1);
            if (created) created->push_back(ib), created->push_back(cm);
            cur = cm;
        }
        return cur;
    }

    // ---------------- Rotate (Alg. 3): r root, z child of r, y child of z; y becomes root
    void rotate(int r, int z, int y) {
        st.rotations++;
        add_im(r, z, -1);
        add_mi(z, y, -1);
        remove_child(z, y);
        if (nd[z].ch.empty()) {
            remove_child(r, z);
            del_node(z);
        } else {
            add_im(r, z, +1);
        }
        root = y;
        vector<int> D = v_diff(nd[r].bag, nd[y].bag);
        if (D.empty()) {  // Merge r into y (Alg. 4)
            st.merges++;
            vector<int> chs = nd[r].ch;
            for (int c : chs) {
                remove_child(r, c);
                attach_inter_table(y, c);
            }
            del_node(r);
        } else {
            vector<int> I = v_inter(nd[r].bag, nd[y].bag);
            int head = build_chain(r, I, D);
            attach_mains(y, I, {head});
        }
    }

    // ---------------- Move (Alg. 2)
    void move_to_root(int x) {
        vector<int> mains{x}, inters;
        int cur = x;
        while (nd[cur].parent >= 0) {
            int z = nd[cur].parent, m = nd[z].parent;
            inters.push_back(z);
            mains.push_back(m);
            cur = m;
        }
        reverse(mains.begin(), mains.end());
        reverse(inters.begin(), inters.end());
        for (size_t s = 0; s + 1 < mains.size(); ++s) rotate(mains[s], inters[s], mains[s + 1]);
    }

    // ---------------- Split (Def. 11, Alg. 5, Alg. 6)
    struct EInfo {
        int mask = 0;
        vector<int> edi;              // editable intersection children
        vector<pair<int, int>> ned;   // non-editable intersection children, assigned label
        vector<int> SV;               // S cap V_x
    };
    struct Item {
        int type;  // 0 = VIRT table, 1 = STRUCT inter node (already counted), 2 = REAL main node
        int id;
        vector<val_t> T;
        vector<int> bag;
    };

    bool split(int r) {
        const vector<int> Br = nd[r].bag;
        int b = (int)Br.size();
        // ---- choose a good partition of the root bag (Def. 10)
        long long best = -1;
        int bestA = INT_MAX, bestX = -1;
        {
            const val_t* E = nd[r].E.data();
            int d[32] = {0}, cnt[4] = {b, 0, 0, 0};
            for (size_t idx = 0; idx < Lp[b]; ++idx) {
                val_t e = E[idx];
                if (e != INF) {
                    int A = e + cnt[0];
                    if (A <= k + 1 && cnt[1] + A < b && cnt[2] + A < b && cnt[3] + A < b) {
                        if (A < bestA || (A == bestA && cnt[0] > bestX)) {
                            bestA = A;
                            bestX = cnt[0];
                            best = (long long)idx;
                        }
                    }
                }
                for (int p = 0; p < b; ++p) {
                    cnt[d[p]]--;
                    if (d[p] < L - 1) {
                        d[p]++;
                        cnt[d[p]]++;
                        break;
                    }
                    d[p] = 0;
                    cnt[0]++;
                }
            }
        }
        if (best < 0) return false;
        st.splits++;
        ++stamp;
        unordered_map<int, EInfo> info;
        {
            long long t = best;
            int mask = 0;
            for (int p = 0; p < b; ++p) {
                int l = (int)(t % L);
                t /= L;
                vlab[Br[p]] = (int8_t)l;
                vstamp[Br[p]] = stamp;
                if (l) mask |= 1 << l;
            }
            info[r].mask = mask;
        }
        auto lab = [&](int v) -> int {
            if (vstamp[v] != stamp) throw logic_error("split: unlabelled vertex");
            return vlab[v];
        };

        // ---- top-down: choose the split inside the editable part only (Lemma 20, 21)
        vector<int> order, stk{r};
        while (!stk.empty()) {
            int x = stk.back();
            stk.pop_back();
            order.push_back(x);
            EInfo& I = info[x];
            for (int c : nd[x].ch) {
                int m = 0;
                for (int v : nd[c].bag) {
                    int l = lab(v);
                    if (l) m |= 1 << l;
                }
                if (__builtin_popcount(m) >= 2) {
                    I.edi.push_back(c);
                    const auto& Bc = nd[c].bag;
                    for (int y : nd[c].ch) {
                        const auto& By = nd[y].bag;
                        int q = diffpos(By, Bc);
                        int v = By[q];
                        size_t base = 0;
                        for (int p = 0; p < (int)By.size(); ++p)
                            if (p != q) base += (size_t)lab(By[p]) * Lp[p];
                        const val_t* Ey = nd[y].E.data();
                        // minimal size first, then maximal |X| (ties -> X); labels restricted to a(z)
                        long long bv = (long long)Ey[base] + 1;
                        int bl = 0;
                        for (int l = 1; l < L; ++l)
                            if ((m >> l) & 1)
                                if ((long long)Ey[base + l * Lp[q]] < bv) bv = Ey[base + l * Lp[q]], bl = l;
                        vlab[v] = (int8_t)bl;
                        vstamp[v] = stamp;
                        info[y].mask = m | (bl ? 1 << bl : 0);
                        stk.push_back(y);
                    }
                } else {
                    I.ned.push_back({c, m ? __builtin_ctz(m) : -1});
                }
            }
        }

        // ---- bottom-up: create copies and their tables (Alg. 5)
        unordered_map<int, array<vector<Item>, 4>> pend;
        for (int oi = (int)order.size() - 1; oi >= 0; --oi) {
            int x = order[oi];
            EInfo& I = info[x];
            const vector<int> Bx = nd[x].bag;
            bool isroot = (x == r);
            vector<int> SV;
            for (int v : Bx)
                if (lab(v) == 0) SV.push_back(v);
            for (int c : I.edi)
                for (int y : nd[c].ch) SV = v_union(SV, info[y].SV);
            I.SV = SV;
            int hl = -1, g = -1;
            if (!isroot) {
                hl = lab(home_vertex(x));
                g = nd[nd[x].parent].parent;
            }
            // remove the editable children from A_x
            for (int c : I.edi) add_im(x, c, -1);
            // non-editable children: all-X ones go to the copy j* with most children
            int cnt[4] = {0, 0, 0, 0};
            for (auto& pr : I.ned)
                if (pr.second > 0) cnt[pr.second]++;
            int jstar = -1;
            for (int l = 1; l < L; ++l)
                if ((I.mask >> l) & 1)
                    if (jstar < 0 || cnt[l] > cnt[jstar]) jstar = l;
            for (auto& pr : I.ned)
                if (pr.second < 0) pr.second = jstar, cnt[jstar]++;
            // A_x now must only contain the children assigned to j*  (correction of Alg. 5, see README)
            for (auto& pr : I.ned)
                if (pr.second != jstar) add_im(x, pr.first, -1);

            for (int i = 1; i < L; ++i) {
                if (!((I.mask >> i) & 1)) continue;
                vector<int> Bp;
                for (int v : Bx)
                    if (lab(v) == 0 || lab(v) == i) Bp.push_back(v);
                bool real = isroot || hl == i;
                vector<val_t> T;
                bool haveT = cnt[i] > 0;
                if (haveT) {
                    if (i == jstar)
                        restrict_table(nd[x].E, Bx, T, Bp);
                    else {
                        legal_init(T, Bp);
                        for (auto& pr : I.ned)
                            if (pr.second == i) proj_add(T, Bp, nd[pr.first].E, nd[pr.first].bag, +1);
                    }
                }
                if (real) {
                    vector<int> Bn = v_union(Bp, SV);
                    if (Bn.size() > Bx.size()) st.size_violations++;
                    int xi = new_node(true, Bn);
                    if (haveT && Bn.size() == Bp.size()) {
                        nd[xi].E.swap(T);
                    } else {
                        legal_init(nd[xi].E, Bn);
                        if (haveT) {
                            proj_add(nd[xi].E, Bn, T, Bp, +1);
                            tfree(T);
                        }
                    }
                    for (auto& pr : I.ned)
                        if (pr.second == i) {
                            remove_child(x, pr.first);
                            add_child(xi, pr.first);
                        }
                    process_pending(xi, pend[x][i]);
                    pend[x][i].clear();
                    if (!isroot)
                        pend[g][i].push_back(Item{2, xi, {}, {}});
                    else
                        rcopy[i] = xi;
                } else {
                    auto& P = pend[g][i];
                    if (haveT) P.push_back(Item{0, -1, std::move(T), Bp});
                    for (auto& pr : I.ned)
                        if (pr.second == i) {
                            remove_child(x, pr.first);
                            P.push_back(Item{1, pr.first, {}, {}});
                        }
                    for (auto& it : pend[x][i]) P.push_back(std::move(it));
                    pend[x][i].clear();
                }
            }
        }
        return finish_split(r, info, order, bestA);
    }

    int rcopy[4];  // copies r_i of the root created by split()

    void process_pending(int xi, vector<Item>& items) {
        for (auto& it : items)
            if (it.type == 0) {
                proj_add(nd[xi].E, nd[xi].bag, it.T, it.bag, +1);
                tfree(it.T);
            }
        for (auto& it : items)
            if (it.type == 1) attach_inter_struct(xi, it.id);
        map<vector<int>, vector<int>> groups;
        for (auto& it : items)
            if (it.type == 2) {
                vector<int> I = v_inter(nd[it.id].bag, nd[xi].bag);
                if (I.size() + 1 != nd[it.id].bag.size()) throw logic_error("split: unique home violated");
                groups[I].push_back(it.id);
            }
        for (auto& kv : groups) attach_mains(xi, kv.first, kv.second);
    }

    bool finish_split(int r, unordered_map<int, EInfo>& info, const vector<int>& order, int bestA) {
        const vector<int>& S = info[r].SV;
        if ((int)S.size() != bestA) throw logic_error("split: |S| differs from DP value");
        int mask = info[r].mask;
        // new root r' with B = S, intersection node S, forget chains to r_i
        int nr = new_node(true, S);
        int is = new_node(false, S);
        vector<int> heads;
        for (int i = 1; i < L; ++i) {
            if (!((mask >> i) & 1)) continue;
            int ri = rcopy[i];
            vector<int> W = v_diff(nd[ri].bag, S);
            heads.push_back(build_chain(ri, S, W));
        }
        legal_init(nd[is].E, S);
        for (int h : heads) {
            add_child(is, h);
            add_mi(is, h, +1);
        }
        legal_init(nd[nr].E, S);
        add_child(nr, is);
        add_im(nr, is, +1);
        // delete the edited nodes
        for (int x : order) {
            for (int c : info[x].edi) del_node(c);
            del_node(x);
        }
        root = nr;
        return true;
    }

    // ---------------- one run of Algorithm 1 for a fixed k
    // returns false iff a maximum bag without split was found (=> tw > k)
    bool run(int kk, long long& splits_guard_limit) {
        k = kk;
        int cw = width();
        for (int w = cw; w >= 2 * k + 2; --w) {
            if (width() < w) continue;
            int newL = (w + 1 >= 3 * k + 4) ? 3 : 4;
            if (newL != L || !tables_valid) {
                set_L(newL);
                recompute_all();
                do_checks("recompute");
            }
            for (auto& x : nd)
                if (x.alive) x.st = 0;
            st.rounds++;
            long long guard = 0;
            if (g_verbose) fprintf(stderr, "c   round w=%d (bags of size %d), L=%d\n", w, w + 1, L);
            while (true) {
                int x = dfs_find(w + 1);
                st.dfs_restarts++;
                if (x < 0) break;
                move_to_root(x);
                do_checks("move");
                for (int i = 0; i < 4; ++i) rcopy[i] = -1;
                if (!split(root)) return false;
                do_checks("split");
                if (++guard > splits_guard_limit) throw runtime_error("too many splits in one round");
            }
            if (width() >= w + 1) throw logic_error("round finished with a bag of size w+1");
        }
        return true;
    }
};

// ------------------------------------------------------------------ main
static void usage() {
    fprintf(stderr,
            "usage: tw2approx [options] [graph.gr]  (reads stdin if no file)\n"
            "  -k K           fixed k: output width <= 2K+1 or report tw > K\n"
            "                 (default: search the smallest k that succeeds)\n"
            "  --init H       initial decomposition: mindeg (default) | random | mixed\n"
            "  --mix Q        for --init mixed: probability of a random elimination step (default 0.05)\n"
            "  --seed S       seed for --init random\n"
            "  --no-lb        do not use the degeneracy lower bound to skip small k\n"
            "  --mem-gb G     memory limit for the DP tables in GB (default 24)\n"
            "  --check        verify structure and all DP tables after every operation (slow)\n"
            "  -o FILE        write the decomposition to FILE (default stdout)\n"
            "  -v             verbose\n");
}

int main(int argc, char** argv) {
    int fixed_k = -1;
    string init = "mindeg", infile, outfile;
    uint64_t seed = 1;
    bool use_lb = true, check = false;
    double mem_gb = 24, mix = 0.05;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        auto need = [&](void) -> string {
            if (i + 1 >= argc) {
                usage();
                exit(2);
            }
            return argv[++i];
        };
        if (a == "-k")
            fixed_k = stoi(need());
        else if (a == "--init")
            init = need();
        else if (a == "--mix")
            mix = stod(need());
        else if (a == "--seed")
            seed = stoull(need());
        else if (a == "--no-lb")
            use_lb = false;
        else if (a == "--mem-gb")
            mem_gb = stod(need());
        else if (a == "--check")
            check = true;
        else if (a == "-o")
            outfile = need();
        else if (a == "-v")
            g_verbose = 1;
        else if (a == "-h" || a == "--help") {
            usage();
            return 0;
        } else if (!a.empty() && a[0] == '-') {
            usage();
            return 2;
        } else
            infile = a;
    }
    if (init != "mindeg" && init != "random" && init != "mixed") {
        usage();
        return 2;
    }
    auto t0 = chrono::steady_clock::now();
    Graph G;
    try {
        if (infile.empty())
            G = read_gr(cin);
        else {
            ifstream f(infile);
            if (!f) throw runtime_error("cannot open " + infile);
            G = read_gr(f);
        }
    } catch (exception& e) {
        fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
    size_t mem_limit_entries = (size_t)(mem_gb * 1e9 / sizeof(val_t));

    Elim e = eliminate(G, init, seed, mix);
    vector<int> core = core_numbers(G);
    int n = G.n;
    vector<int> comp(n, -1);
    int nc = 0;
    for (int i = n - 1; i >= 0; --i) {
        int v = e.order[i];
        comp[v] = e.parent[v] < 0 ? nc++ : comp[e.parent[v]];
    }
    vector<vector<int>> cverts(nc);
    for (int v : e.order) cverts[comp[v]].push_back(v);
    fprintf(stderr, "c graph: n=%d m=%lld components=%d, initial width (%s) = %d\n", n, G.m, nc, init.c_str(),
            e.width);

    vector<vector<int>> out_bags;
    vector<pair<int, int>> out_edges;
    int final_k = -1, lower_bound = 0, out_width = -1;
    Stats tot;
    size_t peak = 0;
    bool failed_fixed = false;
    int prev_root = -1;

    for (int ci = 0; ci < nc; ++ci) {
        const auto& V = cverts[ci];
        int degen = 0;
        for (int v : V) degen = max(degen, core[v]);
        GTD T(G, mem_limit_entries);
        T.check = check;
        T.build(V, e);
        int w0 = T.width();
        int klo, kres = -1, lb = degen;
        long long guard = 50LL * (long long)V.size() + 1000;
        try {
            if (fixed_k >= 0) {
                if (w0 > 4 * fixed_k + 3) {
                    fprintf(stderr,
                            "error: component %d: initial width %d > 4k+3 = %d (Theorem 1 needs a "
                            "4-approximation)\n",
                            ci, w0, 4 * fixed_k + 3);
                    return 2;
                }
                if (2 * fixed_k + 1 >= w0 || T.run(fixed_k, guard))
                    kres = fixed_k;
                else {
                    failed_fixed = true;
                    lb = max(lb, fixed_k + 1);
                    fprintf(stderr, "c component %d: no split found => tw > %d\n", ci, fixed_k);
                }
            } else {
                // Repeated application of Theorem 1.  Invariant: tw >= lb, current width cw.
                // Use k = max(lb, smallest k with cw <= 4k+3).  Success shrinks cw to <= 2k+1
                // (which may allow a smaller k next time), failure proves tw > k.
                // Stops when cw <= 2*lb+1, i.e. a 2-approximation (replaces Bodlaender's
                // reduction used in Corollary 2 to obtain the initial 4-approximation).
                if (!use_lb) lb = 0;
                while (true) {
                    int cw = T.width();
                    klo = max(0, cw / 4);  // smallest k with cw <= 4k+3
                    int kk = max(lb, klo);
                    if (2 * kk + 1 >= cw) {  // only possible for kk == lb
                        kres = kk;
                        break;
                    }
                    if (g_verbose) fprintf(stderr, "c component %d: trying k=%d (width %d)\n", ci, kk, cw);
                    T.tables_valid = false;
                    if (T.run(kk, guard)) {
                        if (kk == lb) {
                            kres = kk;
                            break;
                        }
                    } else {
                        lb = kk + 1;
                        if (g_verbose) fprintf(stderr, "c component %d: tw > %d\n", ci, kk);
                    }
                }
            }
        } catch (exception& ex) {
            fprintf(stderr, "error: component %d: %s\n", ci, ex.what());
            return 3;
        }
        // statistics
        tot.table_ops += T.st.table_ops;
        tot.entries += T.st.entries;
        tot.rotations += T.st.rotations;
        tot.splits += T.st.splits;
        tot.merges += T.st.merges;
        tot.rounds += T.st.rounds;
        tot.recomputes += T.st.recomputes;
        tot.size_violations += T.st.size_violations;
        peak = max(peak, T.peak_entries);
        lower_bound = max(lower_bound, lb);
        if (kres >= 0) final_k = max(final_k, kres);
        int cw = T.width();
        out_width = max(out_width, cw);
        if (g_verbose || nc <= 20)
            fprintf(stderr, "c component %d: |V|=%zu w0=%d -> width %d (k=%d, lb=%d)\n", ci, V.size(), w0, cw,
                    kres, lb);
        // emit main nodes
        unordered_map<int, int> oid;
        vector<int> ord = T.preorder();
        int croot = -1;
        for (int x : ord)
            if (T.nd[x].main) {
                oid[x] = (int)out_bags.size();
                out_bags.push_back(T.nd[x].bag);
                if (x == T.root) croot = oid[x];
            }
        for (int x : ord)
            if (T.nd[x].main && x != T.root) out_edges.push_back({oid[x], oid[T.nd[T.nd[x].parent].parent]});
        if (prev_root >= 0) out_edges.push_back({prev_root, croot});
        prev_root = croot;
    }
    double secs = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    fprintf(stderr,
            "c result: width=%d  k=%d  (=> width <= 2k+1, tw >= %d)  time=%.3fs\n"
            "c stats: rounds=%lld splits=%lld rotations=%lld merges=%lld table_ops=%lld "
            "entries=%.3g recomputes=%lld peak_table_mem=%.3fGB size_violations=%lld\n",
            out_width, final_k, lower_bound, secs, tot.rounds, tot.splits, tot.rotations, tot.merges,
            tot.table_ops, (double)tot.entries, tot.recomputes, peak * 4.0 / 1e9, tot.size_violations);
    if (failed_fixed) {
        printf("c tw > %d\n", fixed_k);
        return 1;
    }
    FILE* out = stdout;
    if (!outfile.empty()) {
        out = fopen(outfile.c_str(), "w");
        if (!out) {
            fprintf(stderr, "error: cannot write %s\n", outfile.c_str());
            return 2;
        }
    }
    size_t maxb = 0;
    for (auto& b : out_bags) maxb = max(maxb, b.size());
    fprintf(out, "c width %d, k %d, tw lower bound %d\n", out_width, final_k, lower_bound);
    fprintf(out, "s td %zu %zu %d\n", out_bags.size(), maxb, n);
    for (size_t i = 0; i < out_bags.size(); ++i) {
        fprintf(out, "b %zu", i + 1);
        for (int v : out_bags[i]) fprintf(out, " %d", v + 1);
        fprintf(out, "\n");
    }
    // PACE .td (Appendix B): the first bag id of a tree edge must be smaller than the second
    for (auto& ed : out_edges)
        fprintf(out, "%d %d\n", min(ed.first, ed.second) + 1, max(ed.first, ed.second) + 1);
    if (out != stdout) fclose(out);
    return 0;
}
