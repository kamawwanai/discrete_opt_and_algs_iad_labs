#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <climits>
#include <iomanip>
#include <numeric>
using namespace std;


class MaxCliqueProblem
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
        if (!fin.is_open())
        {
            cerr << "Error: Cannot open file '" << filename << "'\n";
            return;
        }
        
        string line;
        int vertices = 0, edges = 0;
        int edges_read = 0;
        bool header_found = false;
        
        while (getline(fin, line))
        {
            if (line.empty())
            {
                continue;
            }
            
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
                if (line_input.fail())
                {
                    cerr << "Error: Failed to parse header line in file '" << filename << "'\n";
                    fin.close();
                    return;
                }
                neighbour_sets.resize(vertices);
                header_found = true;
            }
            else if (line[0] == 'e')
            {
                int start, finish;
                line_input >> command >> start >> finish;
                if (line_input.fail())
                {
                    cerr << "Error: Failed to parse edge line in file '" << filename << "': " << line << "\n";
                    continue;
                }
                if (start < 1 || start > vertices || finish < 1 || finish > vertices)
                {
                    cerr << "Warning: Invalid vertex index in file '" << filename << "': " << start << " or " << finish << "\n";
                    continue;
                }
                // Edges in DIMACS file can be repeated, but it is not a problem for our sets
                neighbour_sets[start - 1].insert(finish - 1);
                neighbour_sets[finish - 1].insert(start - 1);
                edges_read++;
            }
        }
        
        fin.close();
        
        if (!header_found)
        {
            cerr << "Error: Header line (starting with 'p') not found in file '" << filename << "'\n";
            return;
        }
        
        if (edges_read != edges)
        {
            cout << "Warning: Expected " << edges << " edges, but read " << edges_read << " edges from file '" << filename << "'\n";
        }
    }

    // GRASP algorithm for Maximum Clique Problem
    // randomization: size of Restricted Candidate List (RCL)
    // iterations: number of GRASP iterations
    void FindClique(int randomization, int iterations)
    {
        static mt19937 generator;
        for (int iteration = 0; iteration < iterations; ++iteration)
        {
            vector<int> clique;
            vector<int> candidates;
            candidates.reserve(neighbour_sets.size());
            for (int i = 0; i < neighbour_sets.size(); ++i)
            {
                candidates.push_back(i);
            }
            
            // use vector<bool> for faster candidate lookup
            vector<bool> is_candidate(neighbour_sets.size(), false);
            for (int c : candidates)
            {
                is_candidate[c] = true;
            }
            
            // diversity: sometimes start with a high-degree vertex to explore different regions
            // 25% of iterations start with a pre-selected high-degree vertex
            if (GetRandom(0, 100) < 25 && !candidates.empty())
            {
                // find vertex with maximum degree in current candidate set
                int max_degree_vertex = candidates[0];
                int max_degree = 0;
                for (int v : candidates)
                {
                    int degree = neighbour_sets[v].size();
                    if (degree > max_degree)
                    {
                        max_degree = degree;
                        max_degree_vertex = v;
                    }
                }
                
                clique.push_back(max_degree_vertex);
                is_candidate[max_degree_vertex] = false;
                
                vector<int> new_candidates;
                new_candidates.reserve(candidates.size());
                const auto& max_vertex_neighbors = neighbour_sets[max_degree_vertex];
                for (int c : candidates)
                {
                    if (c != max_degree_vertex && max_vertex_neighbors.count(c) > 0)
                    {
                        new_candidates.push_back(c);
                    }
                    else
                    {
                        is_candidate[c] = false;
                    }
                }
                candidates = move(new_candidates);
            }
            
            // build clique iteratively until no candidates remain
            // reuse vectors to avoid allocations
            static vector<int> candidate_degrees_static;
            candidate_degrees_static.assign(neighbour_sets.size(), 0);
            
            while (!candidates.empty())
            {
                int max_degree = 0;
                for (int v : candidates)
                {
                    int degree = 0;
                    const auto& neighbors_v = neighbour_sets[v];
                    for (int neighbor : neighbors_v)
                    {
                        degree += is_candidate[neighbor] ? 1 : 0;
                    }
                    candidate_degrees_static[v] = degree;
                    if (degree > max_degree)
                    {
                        max_degree = degree;
                    }
                }
                
                // compute improved scores for all candidates
                static vector<pair<int, int>> candidate_scores;
                candidate_scores.clear();
                candidate_scores.reserve(candidates.size());
                
                // precompute critical_threshold once (same for all vertices)
                int critical_threshold = max(2, max_degree / 4);
                
                for (int v : candidates)
                {
                    int degree = candidate_degrees_static[v];
                    
                    // analyze neighbors to compute additional metrics
                    int critical_neighbors = 0;      // neighbors with very low degree
                    int total_neighbor_degree = 0;    // sum of neighbor degrees
                    int min_neighbor_degree = INT_MAX; // minimum neighbor degree
                    int potential = 0;                // look-ahead
                    int neighbor_count = 0;
                    int potential_limit = min(15, degree); // precompute limit
                    
                    const auto& neighbors_v = neighbour_sets[v];  // use reference to avoid repeated lookups
                    for (int neighbor : neighbors_v)
                    {
                        if (is_candidate[neighbor])
                        {
                            int n_degree = candidate_degrees_static[neighbor];
                            total_neighbor_degree += n_degree;
                            if (n_degree < min_neighbor_degree)
                            {
                                min_neighbor_degree = n_degree;
                            }
                            
                            // Look-ahead
                            if (neighbor_count < potential_limit)
                            {
                                potential += n_degree;
                            }
                            

                            if (n_degree <= critical_threshold)
                            {
                                critical_neighbors += (critical_threshold - n_degree + 1);
                            }
                            
                            neighbor_count++;
                        }
                    }
                    
                    int avg_neighbor_degree = (neighbor_count > 0) ? (total_neighbor_degree / neighbor_count) : 0;
                    if (min_neighbor_degree == INT_MAX) min_neighbor_degree = 0;
                    

                    int score = degree * 1000 
                              + critical_neighbors * 100 
                              + avg_neighbor_degree * 10
                              + min_neighbor_degree * 20
                              + potential * 5;
                    
                    candidate_scores.push_back({v, score});
                }
                
                int rcl_size = min(randomization, (int)candidate_scores.size());
                if (rcl_size == 0) rcl_size = 1;
                
                // use partial_sort
                if (rcl_size < candidate_scores.size())
                {
                    // partial_sort is often faster for small k
                    partial_sort(candidate_scores.begin(), candidate_scores.begin() + rcl_size,
                                 candidate_scores.end(),
                                 [](const pair<int, int>& a, const pair<int, int>& b) {
                                     return a.second > b.second;
                                 });
                }
                else
                {
                    // if rcl_size == size, just sort all
                    sort(candidate_scores.begin(), candidate_scores.end(),
                         [](const pair<int, int>& a, const pair<int, int>& b) {
                             return a.second > b.second;
                         });
                }
                
                // expand RCL to include all candidates with same score
                int min_score_in_rcl = candidate_scores[rcl_size - 1].second;
                int actual_rcl_size = rcl_size;
                for (int i = rcl_size; i < candidate_scores.size(); ++i)
                {
                    if (candidate_scores[i].second == min_score_in_rcl)
                    {
                        actual_rcl_size++;
                    }
                    else
                    {
                        break;
                    }
                }
                
                // weighted random selection from RCL
                int selected_idx;
                if (actual_rcl_size <= 3)
                {
                    // for small RCL, uniform random is fine
                    selected_idx = GetRandom(0, actual_rcl_size - 1);
                }
                else
                {
                    // use linear weighting for simplicity and speed
                    int max_score = candidate_scores[0].second;
                    int min_score_in_rcl = candidate_scores[actual_rcl_size - 1].second;
                    int score_range = max_score - min_score_in_rcl;
                    
                    if (score_range > 0)
                    {
                        // precompute division for efficiency
                        int score_range_div = max(1, score_range);
                        
                        // calculate total weight and cumulative weights in one pass
                        int total_weight = 0;
                        static vector<int> cumulative_weights;
                        cumulative_weights.clear();
                        cumulative_weights.reserve(actual_rcl_size);
                        
                        for (int i = 0; i < actual_rcl_size; ++i)
                        {
                            int normalized_score = candidate_scores[i].second - min_score_in_rcl;
                            int weight = 1 + (normalized_score * 10) / score_range_div;
                            total_weight += weight;
                            cumulative_weights.push_back(total_weight);
                        }
                        
                        // select based on cumulative weights
                        int random_val = GetRandom(0, total_weight - 1);
                        for (int i = 0; i < actual_rcl_size; ++i)
                        {
                            if (random_val < cumulative_weights[i])
                            {
                                selected_idx = i;
                                break;
                            }
                        }
                    }
                    else
                    {
                        // all scores are equal, use uniform random
                        selected_idx = GetRandom(0, actual_rcl_size - 1);
                    }
                }
                
                int selected_vertex = candidate_scores[selected_idx].first;
                
                // add selected vertex to clique
                clique.push_back(selected_vertex);
                
                // update candidates: keep only neighbors of selected vertex (optimized)
                is_candidate[selected_vertex] = false;
                static vector<int> new_candidates;
                new_candidates.clear();
                
                // this is faster when selected vertex has fewer neighbors than candidates
                const auto& selected_neighbors = neighbour_sets[selected_vertex];
                if (selected_neighbors.size() < candidates.size())
                {
                    new_candidates.reserve(selected_neighbors.size());
                    for (int neighbor : selected_neighbors)
                    {
                        if (is_candidate[neighbor])
                        {
                            new_candidates.push_back(neighbor);
                        }
                    }
                    // mark all old candidates as false, then restore true for new ones
                    for (int c : candidates)
                    {
                        is_candidate[c] = false;
                    }
                    for (int c : new_candidates)
                    {
                        is_candidate[c] = true;
                    }
                }
                else
                {
                    new_candidates.reserve(candidates.size());
                    for (int c : candidates)
                    {
                        if (c != selected_vertex && selected_neighbors.count(c) > 0)
                        {
                            new_candidates.push_back(c);
                            // is_candidate[c] remains true
                        }
                        else
                        {
                            is_candidate[c] = false;
                        }
                    }
                }
                candidates = move(new_candidates);
            }

            if (clique.size() > best_clique.size())
            {
                best_clique = clique;
            }
        }
    }
    

    const vector<int>& GetClique()
    {
        return best_clique;
    }

    bool IsGraphValid()
    {
        return !neighbour_sets.empty();
    }

    int GetVertexCount()
    {
        return neighbour_sets.size();
    }

    bool Check()
    {
        if (unique(best_clique.begin(), best_clique.end()) != best_clique.end())
        {
            cout << "Duplicated vertices in the clique\n";
            return false;
        }
        for (int i : best_clique)
        {
            const auto& neighbors_i = neighbour_sets[i];
            for (int j : best_clique)
            {
                if (i != j && neighbors_i.count(j) == 0)
                {
                    cout << "Returned subgraph is not a clique\n";
                    return false;
                }
            }
        }
        return true;
    }

private:
    vector<unordered_set<int>> neighbour_sets;
    vector<int> best_clique;
};

int main()
{
    int iterations;
    cout << "Number of iterations: ";
    cin >> iterations;
    int randomization;
    cout << "Randomization: ";
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
    ofstream fout("clique.csv");
    fout << "File; Clique; Time (sec)\n";
    for (string file : files)
    {
        MaxCliqueProblem problem;
        string filepath = "task2_input/" + file;
        problem.ReadGraphFile(filepath);
        
        if (!problem.IsGraphValid())
        {
            cerr << "Error: Failed to read graph from file '" << filepath << "'. Skipping...\n";
            fout << file << "; ERROR: Failed to read file; N/A\n";
            continue;
        }
        
        // I use different iterations and randomization for specific files
        int current_iterations = iterations;
        int current_randomization = randomization;
        if (file.find("MANN_a27") != string::npos)
        {
            current_iterations = 100;
        }
        else if (file.find("p_hat") != string::npos || file.find("san") != string::npos)
        {
            current_iterations = 300;
            current_randomization = 10;
        }
        
        clock_t start = clock();
        problem.FindClique(current_randomization, current_iterations);
        if (! problem.Check())
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
