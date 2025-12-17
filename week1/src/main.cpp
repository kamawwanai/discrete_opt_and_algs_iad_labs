#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <climits>
#include <cstdint>
#include <time.h>
#include <filesystem>

using namespace std;


class ColoringProblem
{
public:
    int GetRandom(int a, int b)
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
            colors.assign(vertices, 0);
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
    }

    void GreedyGraphColoring()
    {
        const int n = static_cast<int>(neighbour_sets.size());
        colors.assign(n, 0);
        maxcolor = 1;
        
        vector<int> saturation(n, 0);
        vector<bool> used_colors;
        used_colors.reserve(n + 1);
        
        unordered_set<int> uncolored_vertices;
        for (int i = 0; i < n; ++i)
            uncolored_vertices.insert(i);

        while (!uncolored_vertices.empty())
        {
            // Select vertex with highest saturation degree
            int max_saturation = -1;
            int max_degree = -1;
            int max_uncolored_neighbors = -1;
            long long max_neighbor_degree_sum = -1;
            vector<int> candidate_list;

            for (int vertex : uncolored_vertices)
            {
                int sat = saturation[vertex];
                int degree = static_cast<int>(neighbour_sets[vertex].size());
                
                int uncolored_neighbors = 0;
                long long neighbor_degree_sum = 0;
                for (int neighbour : neighbour_sets[vertex])
                {
                    if (colors[neighbour] == 0)
                    {
                        uncolored_neighbors++;
                        neighbor_degree_sum += static_cast<int>(neighbour_sets[neighbour].size());
                    }
                }

                if (sat > max_saturation)
                {
                    max_saturation = sat;
                    max_degree = degree;
                    max_uncolored_neighbors = uncolored_neighbors;
                    max_neighbor_degree_sum = neighbor_degree_sum;
                    candidate_list.clear();
                    candidate_list.push_back(vertex);
                }
                else if (sat == max_saturation)
                {
                    if (degree > max_degree)
                    {
                        max_degree = degree;
                        max_uncolored_neighbors = uncolored_neighbors;
                        max_neighbor_degree_sum = neighbor_degree_sum;
                        candidate_list.clear();
                        candidate_list.push_back(vertex);
                    }
                    else if (degree == max_degree)
                    {
                        if (uncolored_neighbors > max_uncolored_neighbors)
                        {
                            max_uncolored_neighbors = uncolored_neighbors;
                            max_neighbor_degree_sum = neighbor_degree_sum;
                            candidate_list.clear();
                            candidate_list.push_back(vertex);
                        }
                        else if (uncolored_neighbors == max_uncolored_neighbors)
                        {
                            if (neighbor_degree_sum > max_neighbor_degree_sum)
                            {
                                max_neighbor_degree_sum = neighbor_degree_sum;
                                candidate_list.clear();
                                candidate_list.push_back(vertex);
                            }
                            else if (neighbor_degree_sum == max_neighbor_degree_sum)
                            {
                                candidate_list.push_back(vertex);
                            }
                        }
                    }
                }
            }

            int vertex = candidate_list[0];
            
            // Find the smallest available color
            used_colors.assign(maxcolor + 2, false);
            for (int neighbour : neighbour_sets[vertex])
            {
                if (colors[neighbour] != 0)
                {
                    used_colors[colors[neighbour]] = true;
                }
            }
            
            // prefer colors that dont increase neighbor saturation
            // among those, prefer largest existing color
            int best_color = 1;
            int min_saturation_increase = INT_MAX;
            bool found_zero_increase = false;
            
            for (int color = 1; color <= maxcolor; ++color)
            {
                if (used_colors[color])
                    continue;
                
                int saturation_increase = 0;
                for (int neighbour : neighbour_sets[vertex])
                {
                    if (colors[neighbour] == 0)
                    {
                        bool has_color = false;
                        for (int other_neighbour : neighbour_sets[neighbour])
                        {
                            if (colors[other_neighbour] == color)
                            {
                                has_color = true;
                                break;
                            }
                        }
                        if (!has_color)
                            saturation_increase++;
                    }
                }
                
                if (saturation_increase == 0)
                {
                    if (!found_zero_increase || color > best_color)
                    {
                        found_zero_increase = true;
                        best_color = color;
                        min_saturation_increase = 0;
                    }
                }
                else if (!found_zero_increase)
                {
                    if (saturation_increase < min_saturation_increase)
                    {
                        min_saturation_increase = saturation_increase;
                        best_color = color;
                    }
                    else if (saturation_increase == min_saturation_increase && color > best_color)
                    {
                        // Among colors with same saturation increase, prefer larger (better balance)
                        best_color = color;
                    }
                }
            }
            
            if (!found_zero_increase && min_saturation_increase == INT_MAX)
            {
                best_color = maxcolor + 1;
                maxcolor = best_color;
            }

            colors[vertex] = best_color;

            // faster than recalculating for all uncolored vertices
            for (int neighbour : neighbour_sets[vertex])
            {
                if (colors[neighbour] == 0)
                {
                    unordered_set<int> neighbor_colors;
                    for (int other_neighbour : neighbour_sets[neighbour])
                    {
                        if (colors[other_neighbour] != 0)
                        {
                            neighbor_colors.insert(colors[other_neighbour]);
                        }
                    }
                    saturation[neighbour] = static_cast<int>(neighbor_colors.size());
                }
            }

            uncolored_vertices.erase(vertex);
        }
    }

    bool Check()
    {
        for (size_t i = 0; i < neighbour_sets.size(); ++i)
        {
            if (colors[i] == 0)
            {
                cout << "Vertex " << i + 1 << " is not colored\n";
                return false;
            }
            for (int neighbour : neighbour_sets[i])
            {
                if (colors[neighbour] == colors[i])
                {
                    cout << "Neighbour vertices " << i + 1 << ", " << neighbour + 1 <<  " have the same color\n";
                    return false;
                }
            }
        }
        return true;
    }

    int GetNumberOfColors()
    {
        return maxcolor;
    }

    const vector<int>& GetColors()
    {
        return colors;
    }

private:
    vector<int> colors;
    int maxcolor = 1;
    vector<unordered_set<int>> neighbour_sets;
};

int main()
{
    string folder_path = "task1_files";
    vector<string> file_names = { "myciel3.col", "myciel7.col", "school1.col", "school1_nsh.col",
        "anna.col","miles1000.col", "miles1500.col","le450_5a.col",
        "le450_15b.col", "queen11_11.col" };
    vector<string> files;
    for (const auto& name : file_names) {
        files.push_back(folder_path + "/" + name);
    }
    ofstream fout("color.csv");
    fout << "Instance; Colors; Time (sec)\n";
    cout << "Instance; Colors; Time (sec)\n";
    for (string file : files)
    {
        ColoringProblem problem;
        problem.ReadGraphFile(file);
        clock_t start = clock();
        problem.GreedyGraphColoring();
        if (! problem.Check())
        {
            fout << "*** WARNING: incorrect coloring: ***\n";
            cout << "*** WARNING: incorrect coloring: ***\n";
        }
        fout << file << "; " << problem.GetNumberOfColors() << "; " 
            << fixed << setprecision(6) << double(clock() - start) / CLOCKS_PER_SEC << '\n';
        cout << file << "; " << problem.GetNumberOfColors() << "; " 
            << fixed << setprecision(6) << double(clock() - start) / CLOCKS_PER_SEC << '\n';
            }
    fout.close();
    return 0;
}