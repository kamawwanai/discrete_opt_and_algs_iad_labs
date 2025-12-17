#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <iomanip>
#include <cstdint>
#include <limits>
using namespace std;


class MaxCliqueTabuSearch
{
public:
    static int GetRandom(int a, int b)
    {
        static mt19937 generator;
        uniform_int_distribution<int> uniform(a, b);
        return uniform(generator);
    }

    void ReadGraphFile(string filename)
    {
        ifstream fin(filename);
        string line;
        int vertices = 0, edges = 0;
        while (getline(fin, line))
        {
            if (line[0] == 'c')
            {
                continue;
            }

            stringstream line_input(line);
            char command;
            if (line[0] == 'p')
            {
                string type;
                line_input >> command >> type >> vertices >> edges;
                neighbour_sets.resize(vertices);
                qco.resize(vertices);
                index.resize(vertices, -1);
                non_neighbours.resize(vertices);
            }
            else
            {
                int start, finish;
                line_input >> command >> start >> finish;
                // Edges in DIMACS file can be repeated, but it is not a problem for our sets
                neighbour_sets[start - 1].insert(finish - 1);
                neighbour_sets[finish - 1].insert(start - 1);
            }
        }
        for (int i = 0; i < vertices; ++i)
        {
            for (int j = 0; j < vertices; ++j)
            {
                if (neighbour_sets[i].count(j) == 0 && i != j)
                    non_neighbours[i].insert(j);
            }
        }
    }

    void RunSearch(int starts, int randomization)
    {
        for (int iter = 0; iter < starts; ++iter)
        {
            ClearClique();
            for (size_t i = 0; i < neighbour_sets.size(); ++i)
            {
                qco[i] = i;
                index[i] = i;
            }
            RunInitialHeuristic(randomization);
            c_border = q_border;
            int swaps = 0;
            while (swaps < 100)
            {
                if (! Move())
                {
                    if (! Swap1To1())
                    {
                        break;
                    }
                    else
                    {
                        ++swaps;
                    }
                }
            }
            if (q_border > best_clique.size())
            {
                best_clique.clear();
                for (int i = 0; i < q_border; ++i)
                    best_clique.insert(qco[i]);
            }
        }
    }

    const unordered_set<int>& GetClique()
    {
        return best_clique;
    }

    bool Check()
    {
        for (int i : best_clique)
        {
            for (int j : best_clique)
            {
                if (i != j && neighbour_sets[i].count(j) == 0)
                {
                    cout << "Returned subgraph is not a clique\n";
                    return false;
                }
            }
        }
        return true;
    }

    void ClearClique()
    {
        best_clique.clear();
        q_border = 0;
        c_border = 0;
    }

private:
    int ComputeTightness(int vertex)
    {
        int tightness = 0;
        for (int i = 0; i < q_border; ++i)
        {
            if (neighbour_sets[qco[i]].count(vertex) == 0)
                ++tightness;
        }
        return tightness;
    }

    void SwapVertices(int vertex, int border)
    {
        int vertex_at_border = qco[border];
        swap(qco[index[vertex]], qco[border]);
        swap(index[vertex], index[vertex_at_border]);
    }

    void InsertToClique(int i)
    {
        for (int j : non_neighbours[i])
        {
            if (ComputeTightness(j) == 0)
            {
                --c_border;
                SwapVertices(j, c_border);
            }
        }
        SwapVertices(i, q_border);
        ++q_border;
    }

    void RemoveFromClique(int k)
    {
        for (int j : non_neighbours[k])
        {
            if (ComputeTightness(j) == 1)
            {
                SwapVertices(j, c_border);
                c_border++;
            }
        }
        --q_border;
        SwapVertices(k, q_border);
    }

    bool Swap1To1()
    {
        int st = GetRandom(0, q_border - 1);
        for (int counter = 0; counter < q_border; ++counter)
        {
            int vertex_index = (counter + st) % q_border;
            int vertex = qco[vertex_index];
            vector<int> L;
            for (int i : non_neighbours[vertex])
            {
                if (ComputeTightness(i) == 1)
                {
                    L.push_back(i);
                }
            }
            if (L.empty())
                continue;
            int index_in_l = GetRandom(0, L.size() - 1);
            int change = L[index_in_l];
            RemoveFromClique(vertex);
            InsertToClique(change);
            return true;
        }
        return false;
    }

    bool Move()
    {
        if (c_border == q_border)
            return false;
        int index_in_qco = GetRandom(q_border, c_border - 1);
        int vertex = qco[index_in_qco];
        InsertToClique(vertex);
        return true;
    }

    void RunInitialHeuristic(int randomization)
    {
        static mt19937 generator;
        vector<int> candidates(neighbour_sets.size());
        for (size_t i = 0; i < neighbour_sets.size(); ++i)
        {
            candidates[i] = i;
        }
        shuffle(candidates.begin(), candidates.end(), generator);
        while (! candidates.empty())
        {
            int last = candidates.size() - 1;
            int rnd = GetRandom(0, min(randomization - 1, last));
            int vertex = candidates[rnd];
            SwapVertices(vertex, q_border);
            ++q_border;
            for (int c = 0; c < candidates.size(); ++c)
            {
                int candidate = candidates[c];
                if (neighbour_sets[vertex].count(candidate) == 0)
                {
                    swap(candidates[c], candidates[candidates.size() - 1]);
                    candidates.pop_back();
                    --c;
                }
            }
            shuffle(candidates.begin(), candidates.end(), generator);
        }
    }

private:
    vector<unordered_set<int>> neighbour_sets;
    vector<unordered_set<int>> non_neighbours;
    unordered_set<int> best_clique;
    vector<int> qco;
    vector<int> index;
    int q_border = 0;
    int c_border = 0;
};


class BnBSolver
{
public:
    void ReadGraphFile(string filename)
    {
        file = filename;
        ifstream fin(filename);
        if (!fin)
            return;

        int vert = 0, edges = 0;
        char tag = 0;
        while (fin >> tag)
        {
            if (tag == 'c')
            {
                fin.ignore(numeric_limits<streamsize>::max(), '\n');
                continue;
            }
            if (tag == 'p')
            {
                string type;
                fin >> type >> vert >> edges;
                n = vert;
                words = (n + 63) / 64;
                adj.assign(static_cast<size_t>(n) * static_cast<size_t>(words), 0ULL);
                degree.assign(n, 0);
                continue;
            }
            if (tag == 'e')
            {
                int st = 0, fn = 0;
                fin >> st >> fn;
                --st; --fn;
                if (st < 0 || fn < 0 || st >= n || fn >= n || st == fn)
                    continue;
                if (!IsAdjacent(st, fn))
                {
                    SetEdge(st, fn);
                    SetEdge(fn, st);
                    ++degree[st];
                    ++degree[fn];
                }
                continue;
            }

            fin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
    }

    void RunBnB()
    {
        // Fast initial lower bound
        InitBestCliqueGreedy(8, 10);

        vector<int> candidates(n);
        for (int i = 0; i < n; ++i) candidates[i] = i;
        sort(candidates.begin(), candidates.end(),
             [&](int a, int b) { return degree[a] > degree[b]; });

        cur_clique.clear();
        Expand(candidates);
    }

    const unordered_set<int>& GetClique()
    {
        return best_clique_set;
    }

    bool Check()
    {
        for (int i : best_clique_set)
        {
            for (int j : best_clique_set)
            {
                if (i != j && !IsAdjacent(i, j))
                {
                    cout << "Returned subgraph is not clique\n";
                    return false;
                }
            }
        }
        return true;
    }

    void ClearClique()
    {
        best_clique_set.clear();
        best_clique_vec.clear();
        cur_clique.clear();
        best_size = 0;
    }

private:
    inline const uint64_t* AdjRow(int u) const
    {
        return adj.data() + static_cast<size_t>(u) * static_cast<size_t>(words);
    }

    inline uint64_t* AdjRow(int u)
    {
        return adj.data() + static_cast<size_t>(u) * static_cast<size_t>(words);
    }

    inline bool IsAdjacent(int u, int v) const
    {
        const uint64_t mask = 1ULL << (v & 63);
        return (AdjRow(u)[v >> 6] & mask) != 0ULL;
    }

    inline void SetEdge(int u, int v)
    {
        AdjRow(u)[v >> 6] |= (1ULL << (v & 63));
    }

    static inline void OrEq(uint64_t* a, const uint64_t* b, int w)
    {
        for (int i = 0; i < w; ++i) a[i] |= b[i];
    }

    static inline void SetBit(uint64_t* bs, int v)
    {
        bs[v >> 6] |= (1ULL << (v & 63));
    }

    static inline bool TestBit(const uint64_t* bs, int v)
    {
        return (bs[v >> 6] & (1ULL << (v & 63))) != 0ULL;
    }

    void ColorSort(const vector<int>& candidates, vector<int>& order, vector<int>& bounds) const
    {
        order.clear();
        bounds.clear();
        order.reserve(candidates.size());
        bounds.reserve(candidates.size());

        vector<int> U = candidates;
        int color = 0;
        vector<uint64_t> forbidden(static_cast<size_t>(words), 0ULL);
        vector<int> nextU;
        nextU.reserve(U.size());
        while (!U.empty())
        {
            ++color;
            fill(forbidden.begin(), forbidden.end(), 0ULL);
            nextU.clear();

            for (int v : U)
            {
                if (!TestBit(forbidden.data(), v))
                {
                    order.push_back(v);
                    bounds.push_back(color);
                    SetBit(forbidden.data(), v);
                    OrEq(forbidden.data(), AdjRow(v), words);
                }
                else
                {
                    nextU.push_back(v);
                }
            }
            U.swap(nextU);
        }
    }

    void Expand(const vector<int>& candidates)
    {
        if (candidates.empty())
        {
            if (static_cast<int>(cur_clique.size()) > best_size)
            {
                best_size = static_cast<int>(cur_clique.size());
                best_clique_vec = cur_clique;
                best_clique_set.clear();
                for (int v : best_clique_vec) best_clique_set.insert(v);
            }
            return;
        }

        vector<int> order;
        vector<int> bounds;
        ColorSort(candidates, order, bounds);

        for (int i = static_cast<int>(order.size()) - 1; i >= 0; --i)
        {
            if (static_cast<int>(cur_clique.size()) + bounds[i] <= best_size)
                return; // since bounds are nondecreasing for prefixes

            const int v = order[i];
            const uint64_t* row_v = AdjRow(v);
            vector<int> new_candidates;
            new_candidates.reserve(static_cast<size_t>(i));
            for (int j = 0; j < i; ++j)
            {
                const int u = order[j];
                if ((row_v[u >> 6] & (1ULL << (u & 63))) != 0ULL) new_candidates.push_back(u);
            }

            cur_clique.push_back(v);
            Expand(new_candidates);
            cur_clique.pop_back();
        }
    }

private:
    void InitBestCliqueGreedy(int starts, int rcl)
    {
        if (n <= 0) return;
        if (starts <= 0) starts = 1;
        if (rcl <= 0) rcl = 1;

        vector<int> base(n);
        for (int i = 0; i < n; ++i) base[i] = i;
        sort(base.begin(), base.end(), [&](int a, int b) { return degree[a] > degree[b]; });

        static mt19937 gen(1234567);

        for (int s = 0; s < starts; ++s)
        {
            vector<int> cand = base;
            vector<int> clique;
            clique.reserve(64);

            while (!cand.empty())
            {
                const int last = static_cast<int>(cand.size()) - 1;
                const int take = min(rcl - 1, last);
                uniform_int_distribution<int> dist(0, take);
                const int pick_idx = dist(gen);
                const int v = cand[pick_idx];
                const uint64_t* row_v = AdjRow(v);

                clique.push_back(v);

                vector<int> next;
                next.reserve(static_cast<size_t>(last));
                for (int i = 0; i < static_cast<int>(cand.size()); ++i)
                {
                    if (i == pick_idx) continue;
                    const int u = cand[i];
                    if ((row_v[u >> 6] & (1ULL << (u & 63))) != 0ULL) next.push_back(u);
                }
                cand.swap(next);
            }

            if (static_cast<int>(clique.size()) > best_size)
            {
                best_size = static_cast<int>(clique.size());
                best_clique_vec = clique;
                best_clique_set.clear();
                for (int v : best_clique_vec) best_clique_set.insert(v);
            }
        }
    }

    int n = 0;
    int words = 0;
    vector<uint64_t> adj;
    vector<int> degree;

    int best_size = 0;
    vector<int> best_clique_vec;
    vector<int> cur_clique;
    unordered_set<int> best_clique_set;

    string file;
};

int main(int argc, char** argv)
{
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    vector<string> files;
    if (argc > 1)
    {
        for (int i = 1; i < argc; ++i) files.emplace_back(argv[i]);
    }
    else
    {
        files = {
            "brock200_1.clq", "brock200_2.clq", "brock200_3.clq", "brock200_4.clq",
            "C125.9.clq", "gen200_p0.9_44.clq", "gen200_p0.9_55.clq",
            "hamming8-4.clq", "johnson16-2-4.clq", "johnson8-2-4.clq",
            "keller4.clq", "MANN_a27.clq", "MANN_a9.clq",
            "p_hat1000-1.clq", "p_hat1500-1.clq",
            "p_hat300-3.clq", "san1000.clq", "sanr200_0.9.clq"
        };
    }
    ofstream fout("clique_bnb.csv");
    fout << "File; Clique; Time (sec)\n";
    for (string file : files)
    {
        BnBSolver problem;
        string filepath = file;
        if (filepath.find('/') == string::npos && filepath.find('\\') == string::npos)
            filepath = "task4_input/" + filepath;
        problem.ReadGraphFile(filepath);
        problem.ClearClique();
        clock_t start = clock();
        problem.RunBnB();
        if (! problem.Check())
        {
            cout << "*** WARNING: incorrect clique ***\n";
            fout << "*** WARNING: incorrect clique ***\n";
        }

        double time_sec = double(clock() - start) / CLOCKS_PER_SEC;

        fout << file << "; " << problem.GetClique().size() << "; " << fixed << setprecision(6) << time_sec << '\n';
        cout << file << ", result - " << problem.GetClique().size() << ", time - " << fixed << setprecision(6) << time_sec << '\n';
    }
    
    return 0;
}
