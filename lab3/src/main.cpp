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
using namespace std;


class MaxCliqueTabuSearch
{
public:
    MaxCliqueTabuSearch() : rng(random_device{}()) {}

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
                non_neighbours.assign(vertices, {});
                degrees.assign(vertices, 0);
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
            degrees[i] = static_cast<int>(neighbour_sets[i].size());

        // build explicit non-neighbour lists (used for O(deg_non) tightness updates)
        for (int i = 0; i < vertices; ++i)
        {
            non_neighbours[i].reserve(vertices - degrees[i] - 1);
            for (int j = 0; j < vertices; ++j)
            {
                if (i != j && neighbour_sets[i].count(j) == 0)
                    non_neighbours[i].push_back(j);
            }
        }
    }

    // iterations: number of tabu steps per instance
    // randomization: tournament size for choosing candidates
    void RunSearch(int iterations, int randomization)
    {
        if (neighbour_sets.empty())
            return;

        randomization = max(1, randomization);
        iterations = max(1, iterations);

        const int n = static_cast<int>(neighbour_sets.size());
        
        const int num_restarts = min(3, max(1, n / 200));
        int global_best_size = 0;
        vector<int> global_best_vector;

        for (int restart = 0; restart < num_restarts; ++restart)
        {
            InitializeState(n);
            
            // initial maximal clique (randomized greedy)
            BuildInitialClique(randomization);
            best_vector = clique;

            int best_size = static_cast<int>(best_vector.size());
            int current_size = best_size;
            int stall = 0;
            const int stall_limit = max(500, 5 * n);
            int base_tenure = 7;
            int last_improvement = 0;
            int intensification_count = 0;
            bool in_intensification = false;

            // adaptive parameters
            double improvement_rate = 1.0;
            int consecutive_swaps = 0;

            for (int it = 1; it <= iterations / num_restarts; ++it)
            {
                // expand to a maximal clique with non-tabu insertions
                ExpandGreedy(randomization, it, best_size, current_size);

                current_size = static_cast<int>(clique.size());

                if (current_size > best_size)
                {
                    best_size = current_size;
                    best_vector = clique;
                    stall = 0;
                    last_improvement = it;
                    in_intensification = false;
                    intensification_count = 0;
                    improvement_rate = 1.0;
                    
                    // reduce when finding improvements
                    base_tenure = max(5, base_tenure - 1);
                }
                else
                {
                    ++stall;
                    improvement_rate = 0.95 * improvement_rate;
                }

                // when we find a good solution, intensify search
                if (current_size >= best_size - 1 && !in_intensification && it - last_improvement < 100)
                {
                    in_intensification = true;
                    intensification_count = 0;
                }

                if (in_intensification)
                {
                    ++intensification_count;
                    if (intensification_count > 200 || current_size < best_size - 2)
                    {
                        in_intensification = false;
                        intensification_count = 0;
                    }
                }

                // trigger earlier if improvement rate is low
                int adaptive_stall_limit = stall_limit;
                if (improvement_rate < 0.3 && it > iterations / (2 * num_restarts))
                {
                    adaptive_stall_limit = stall_limit / 2;
                }

                if (stall >= adaptive_stall_limit)
                {
                    // remove a chunk of vertices, then rebuild
                    Diversify(it, base_tenure);
                    stall = 0;
                    base_tenure = min(30, base_tenure + 2);
                    improvement_rate = 1.0;
                    in_intensification = false;
                    continue;
                }

                // if no swap exists, drop one vertex
                if (HasAdmissibleAdd(it, best_size, current_size))
                {
                    consecutive_swaps = 0;
                    continue;
                }

                if (TrySwap(it, randomization, base_tenure, best_size, current_size))
                {
                    ++consecutive_swaps;
                    // if too many swaps in a row, drop one to diversify
                    if (consecutive_swaps > 50)
                    {
                        DropOne(it, base_tenure);
                        consecutive_swaps = 0;
                    }
                    continue;
                }

                consecutive_swaps = 0;
                DropOne(it, base_tenure);
            }

            if (best_size > global_best_size)
            {
                global_best_size = best_size;
                global_best_vector = best_vector;
            }
        }

        best_vector = global_best_vector;
        best_clique.clear();
        for (int v : best_vector)
            best_clique.insert(v);
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
                    cout << "Returned subgraph is not clique\n";
                    return false;
                }
            }
        }
        return true;
    }

    void ClearClique()
    {
        clique.clear();
    }

private:
    struct FastSet
    {
        vector<int> items;
        vector<int> pos; // -1 means not in set

        void Init(int n)
        {
            items.clear();
            pos.assign(n, -1);
        }

        bool Contains(int v) const { return pos[v] != -1; }
        int Size() const { return static_cast<int>(items.size()); }

        void Add(int v)
        {
            if (pos[v] != -1)
                return;
            pos[v] = static_cast<int>(items.size());
            items.push_back(v);
        }

        void Remove(int v)
        {
            int p = pos[v];
            if (p == -1)
                return;
            int last = items.back();
            items[p] = last;
            pos[last] = p;
            items.pop_back();
            pos[v] = -1;
        }
    };

    vector<unordered_set<int>> neighbour_sets;
    vector<vector<int>> non_neighbours;
    unordered_set<int> best_clique;
    vector<int> degrees;

    // search state
    mt19937 rng;
    vector<int> clique;
    vector<int> best_vector;
    vector<int> pos_in_clique;
    vector<char> in_clique;
    vector<int> tight;      // tight[v] = # of clique vertices non-adjacent to v
    FastSet C0;             // tight==0 and not in clique
    FastSet C1;             // tight==1 and not in clique
    vector<int> tabu_until; // forbids inserting vertex v while it > current_iter
    vector<int> freq;       // insertion frequency

    int RandInt(int a, int b)
    {
        uniform_int_distribution<int> dist(a, b);
        return dist(rng);
    }

    void InitializeState(int n)
    {
        clique.clear();
        best_vector.clear();
        pos_in_clique.assign(n, -1);
        in_clique.assign(n, 0);
        tight.assign(n, 0);
        C0.Init(n);
        C1.Init(n);
        tabu_until.assign(n, 0);
        freq.assign(n, 0);

        // empty clique => every vertex is feasible to add
        for (int v = 0; v < n; ++v)
            C0.Add(v);
    }

    void AddToClique(int v)
    {
        C0.Remove(v);
        C1.Remove(v);

        in_clique[v] = 1;
        pos_in_clique[v] = static_cast<int>(clique.size());
        clique.push_back(v);

        ++freq[v];
        tight[v] = 0;

        for (int u : non_neighbours[v])
        {
            if (in_clique[u])
                continue;
            int old = tight[u];
            ++tight[u];
            if (old == 0)
            {
                C0.Remove(u);
                C1.Add(u);
            }
            else if (old == 1)
            {
                C1.Remove(u);
            }
        }
    }

    void RemoveFromClique(int v)
    {
        int p = pos_in_clique[v];
        if (p == -1)
            return;

        int last = clique.back();
        clique[p] = last;
        pos_in_clique[last] = p;
        clique.pop_back();
        pos_in_clique[v] = -1;
        in_clique[v] = 0;

        // v was in the clique, so now its tightness becomes 0
        tight[v] = 0;
        C0.Add(v);
        C1.Remove(v);

        for (int u : non_neighbours[v])
        {
            if (in_clique[u])
                continue;
            int old = tight[u];
            --tight[u];
            if (old == 1)
            {
                C1.Remove(u);
                C0.Add(u);
            }
            else if (old == 2)
            {
                C1.Add(u);
            }
        }
    }

    int CurrentTenure(int base_tenure) const
    {
        int t = base_tenure + static_cast<int>(clique.size()) / 8;
        return min(60, max(5, t));
    }

    bool IsTabuToInsert(int v, int iter) const
    {
        return tabu_until[v] > iter;
    }

    // improved aspiration: allow tabu moves that lead to significant improvements
    bool IsAspiration(int v, int iter, int best_size, int current_size) const
    {
        const int target = static_cast<int>(clique.size()) + 1;
        // improves best
        if (target > best_size)
            return true;
        // significantly improves current (within 1 of best)
        if (target >= best_size && current_size < best_size - 1)
            return true;
        return false;
    }

    bool HasAdmissibleAdd(int iter, int best_size, int current_size) const
    {
        if (C0.Size() == 0)
            return false;
        for (int v : C0.items)
        {
            if (!IsTabuToInsert(v, iter) || IsAspiration(v, iter, best_size, current_size))
                return true;
        }
        return false;
    }

    int ChooseFromC0Tournament(int k, int iter, int best_size, int current_size) 
    {
        if (C0.Size() == 0)
            return -1;

        k = min(k, C0.Size());
        // if too small, just return first feasible
        int tries = min(C0.Size(), max(50, k * 15));
        int best_v = -1;
        double best_score = -1e100;

        auto neighbors_in_C0 = [&](int v) -> int
        {
            int cnt = 0;
            for (int u : neighbour_sets[v])
            {
                if (!in_clique[u] && tight[u] == 0)
                    ++cnt;
            }
            return cnt;
        };

        // tournament: sample multiple random candidates from C0
        for (int i = 0; i < tries; ++i)
        {
            int idx = RandInt(0, C0.Size() - 1);
            int v = C0.items[idx];
            bool is_tabu = IsTabuToInsert(v, iter);
            if (is_tabu && !IsAspiration(v, iter, best_size, current_size))
                continue;

            // prefer vertices with high C0 degree and high degree, penalize frequency
            int c0deg = neighbors_in_C0(v);
            
            // improved scoring: more weight to C0 degree, less to frequency
            double score = 2000.0 * static_cast<double>(c0deg) 
                         + 1.5 * static_cast<double>(degrees[v]) 
                         - 0.1 * static_cast<double>(freq[v]);
            
            // slight penalty for tabu moves (even if allowed by aspiration)
            if (is_tabu)
                score *= 0.95;
            
            // small random tie-breaker
            score += 1e-6 * RandInt(0, 1000000);

            if (score > best_score)
            {
                best_score = score;
                best_v = v;
            }
        }
        return best_v;
    }

    void ExpandGreedy(int randomization, int iter, int best_size, int current_size)
    {
        int max_expansions = 1000;
        int expansions = 0;
        while (expansions < max_expansions)
        {
            int v = ChooseFromC0Tournament(randomization, iter, best_size, current_size);
            if (v == -1)
                return;
            AddToClique(v);
            ++expansions;
            current_size = static_cast<int>(clique.size());
        }
    }

    void BuildInitialClique(int randomization)
    {
        while (C0.Size() > 0)
        {
            int k = min(randomization, C0.Size());
            int tries = min(C0.Size(), max(60, k * 12));
            int best_v = -1;
            double best_score = -1e100;

            auto neighbors_in_C0 = [&](int v) -> int
            {
                int cnt = 0;
                for (int u : neighbour_sets[v])
                {
                    if (!in_clique[u] && tight[u] == 0)
                        ++cnt;
                }
                return cnt;
            };

            for (int i = 0; i < tries; ++i)
            {
                int idx = RandInt(0, C0.Size() - 1);
                int v = C0.items[idx];
                int c0deg = neighbors_in_C0(v);
                
                // improved scoring for initial construction
                double score = 1500.0 * static_cast<double>(c0deg) 
                             + 1.2 * static_cast<double>(degrees[v]) 
                             + 1e-6 * RandInt(0, 1000000);
                
                if (score > best_score)
                {
                    best_score = score;
                    best_v = v;
                }
            }
            if (best_v == -1)
                break;
            AddToClique(best_v);
        }
    }

    int FindConflictVertexForC1(int v) const
    {
        // v has tight[v]==1, so there is exactly one vertex in the clique not adjacent to v
        for (int u : clique)
        {
            if (neighbour_sets[v].count(u) == 0)
                return u;
        }
        return -1;
    }

    int SwapDeltaC0(int remove_u, int add_v) const
    {
        int gain = 0;
        for (int w : non_neighbours[remove_u])
        {
            if (!in_clique[w] && tight[w] == 1)
                ++gain;
        }

        int loss = 0;
        for (int w : non_neighbours[add_v])
        {
            if (!in_clique[w] && tight[w] == 0)
                ++loss;
        }
        return gain - loss;
    }

    double EvaluateSwap(int remove_u, int add_v, int iter, int best_size, int current_size) const
    {
        int delta = SwapDeltaC0(remove_u, add_v);
        
        // base score from C0 delta
        double score = 100.0 * static_cast<double>(delta);
        
        // prefer swaps that maintain or improve solution quality
        const int new_size = static_cast<int>(clique.size()); // same size after swap
        if (new_size >= best_size - 1)
            score += 50.0;
        
        // diversification: prefer adding vertices with low frequency and removing high frequency
        score -= 0.5 * static_cast<double>(freq[add_v]);
        score -= 0.2 * static_cast<double>(freq[remove_u]);
        
        // prefer removing vertices with lower degree (easier to replace)
        score += 0.1 * static_cast<double>(degrees[remove_u]);
        
        // prefer adding vertices with higher degree (more connections)
        score += 0.15 * static_cast<double>(degrees[add_v]);
        
        return score;
    }

    bool TrySwap(int iter, int randomization, int base_tenure, int best_size, int current_size)
    {
        if (C1.Size() == 0 || clique.empty())
            return false;

        int k = min(randomization, C1.Size());
        int tries = min(C1.Size(), max(80, k * 15));
        int best_v = -1;
        int best_u = -1;
        double best_score = -1e100;

        for (int i = 0; i < tries; ++i)
        {
            int idx = RandInt(0, C1.Size() - 1);
            int v = C1.items[idx];
            bool is_tabu = IsTabuToInsert(v, iter);
            if (is_tabu && !IsAspiration(v, iter, best_size, current_size))
                continue;

            int u = FindConflictVertexForC1(v);
            if (u == -1)
                continue;

            double score = EvaluateSwap(u, v, iter, best_size, current_size);
            
            if (is_tabu)
                score *= 0.9;
            
            score += 1e-6 * RandInt(0, 1000000);

            if (score > best_score)
            {
                best_score = score;
                best_v = v;
                best_u = u;
            }
        }

        if (best_v == -1)
            return false;

        RemoveFromClique(best_u);
        tabu_until[best_u] = iter + CurrentTenure(base_tenure) + RandInt(0, base_tenure);
        AddToClique(best_v);
        return true;
    }

    void DropOne(int iter, int base_tenure)
    {
        if (clique.empty())
            return;

        int k = min(max(10, static_cast<int>(clique.size()) / 2), static_cast<int>(clique.size()));
        int best_u = clique[RandInt(0, static_cast<int>(clique.size()) - 1)];
        double best_score = -1e100;

        for (int i = 0; i < k; ++i)
        {
            int u = clique[RandInt(0, static_cast<int>(clique.size()) - 1)];
            int gain = 0;
            for (int w : non_neighbours[u])
            {
                if (!in_clique[w] && tight[w] == 1)
                    ++gain;
            }

            double score = 10.0 * static_cast<double>(gain) 
                         + 0.3 * static_cast<double>(freq[u])
                         - 0.05 * static_cast<double>(degrees[u]);
            
            score += 1e-6 * RandInt(0, 1000000);
            
            if (score > best_score)
            {
                best_score = score;
                best_u = u;
            }
        }

        RemoveFromClique(best_u);
        tabu_until[best_u] = iter + CurrentTenure(base_tenure) + RandInt(0, base_tenure);
    }

    void Diversify(int iter, int base_tenure)
    {
        if (clique.empty())
            return;

        int remove_cnt = max(1, static_cast<int>(clique.size()) / 2);
        remove_cnt = min(remove_cnt, static_cast<int>(clique.size()));

        vector<pair<int, int>> candidates; // (freq, vertex)
        for (int v : clique)
        {
            candidates.push_back({freq[v], v});
        }
        sort(candidates.begin(), candidates.end(), greater<pair<int, int>>());

        for (int i = 0; i < remove_cnt && i < static_cast<int>(candidates.size()); ++i)
        {
            // Prefer removing high-frequency vertices to escape well-trodden areas
            int idx = i;
            if (i < remove_cnt - 1 && RandInt(0, 100) < 30) // 30% chance to randomize
            {
                idx = RandInt(i, min(i + 3, static_cast<int>(candidates.size()) - 1));
            }
            
            int best_u = candidates[idx].second;

            RemoveFromClique(best_u);
            // longer tabu tenure for diversification
            tabu_until[best_u] = iter + 3 * CurrentTenure(base_tenure) + RandInt(0, 2 * base_tenure);

            if (clique.empty())
                break;
        }
    }
};

int main()
{
    int iterations;
    cout << "Number of iterations (tabu steps): ";
    cin >> iterations;
    int randomization;
    cout << "Randomization (tournament size): ";
    cin >> randomization;
    vector<string> files = { 
        "brock200_1.clq", "brock200_2.clq", "brock200_3.clq", "brock200_4.clq",
        "brock400_1.clq", "brock400_2.clq", "brock400_3.clq", "brock400_4.clq",
        "C125.9.clq", "gen200_p0.9_44.clq", "gen200_p0.9_55.clq",
        "hamming8-4.clq", "johnson16-2-4.clq", "johnson8-2-4.clq",
        "keller4.clq", "MANN_a27.clq", "MANN_a9.clq",
        "p_hat1000-1.clq", "p_hat1000-2.clq", "p_hat1500-1.clq",
        "p_hat300-3.clq", "p_hat500-3.clq", "san1000.clq",
        "sanr200_0.9.clq", "sanr400_0.7.clq"
    };
    ofstream fout("clique_local.csv");
    fout << "File; Clique; Time (sec)\n";
    for (const string& file : files)
    {
        MaxCliqueTabuSearch problem;
        string filepath = "task3_input/" + file;
        problem.ReadGraphFile(filepath);
        clock_t start = clock();
        problem.RunSearch(iterations, randomization);
        if (!problem.Check())
        {
            cout << "*** WARNING: incorrect clique ***\n";
            fout << "*** WARNING: incorrect clique ***\n";
        }
        
        double time_sec = double(clock() - start) / CLOCKS_PER_SEC;

        fout << file << "; " << problem.GetClique().size() << "; " << fixed << setprecision(6) << time_sec << '\n';
        cout << file << ", result - " << problem.GetClique().size() << ", time - " << fixed << setprecision(6) << time_sec << '\n';
    }
    fout.close();
    return 0;
}
