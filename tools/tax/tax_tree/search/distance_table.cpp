#include <iostream>
#include <chrono>
#include <vector>
#include <thread>
#include <omp.h>
#include "file_list_loader.h"
#include "config_distance_table.h"

using namespace std;
using namespace std::chrono;

typedef vector<double> DistanceTo;

const double INF_DISTANCE = 1000000000;

FileListLoader::Files files; // todo: remove

struct MergePair
{
    int a, b;
    double distance;
    MergePair(int a, int b, double distance = INF_DISTANCE) : a(a), b(b), distance(distance){}     

    bool valid() const { return distance < INF_DISTANCE; }
};

struct DistanceTable
{
    struct Node
    {
        DistanceTo distance_to;
        MergePair merge;
        bool removed;

        Node() : merge(-1, -1), removed(false) {}

    };

    vector<Node> nodes;
    DistanceTable(int nodes_count = 0): nodes(nodes_count){}
        
    void set_distance(int node, DistanceTo &d)
    {
        if (!valid_node(node))
            throw std::runtime_error("set_distance invalid node");

        if (d.empty()) // todo: remove - wourkaround for non processed
        {
            throw std::runtime_error("no distance info!");
            nodes[node].distance_to.resize(nodes.size());
            std::fill(nodes[node].distance_to.begin(), nodes[node].distance_to.end(), INF_DISTANCE);
            nodes[node].removed = true;

            return;
        }
            

//        cout << d.size() << " " << nodes.size() << endl;
        if (d.size() != nodes.size())
            throw std::runtime_error("d.size() != distance_to.size()");

        nodes[node].distance_to = d;
    }

    bool valid_node(int node) const
    {
        return node >=0 && node < int(nodes.size());
    }

    double distance(int a, int b) const
    {
        if (!valid_node(a) || !valid_node(b))
            throw "distance !valid(a) || !valid(b)";

        if (nodes[a].distance_to.size() <= b)
            throw "distance:: nodes[a].distance_to.size() <= b";

        return nodes[a].distance_to[b];
    }

    void add_node(const MergePair &pair)
    {
        Node n;
        n.merge = pair;

        nodes[pair.a].removed = true;
        nodes[pair.b].removed = true;

        add_distance_column(n);
        calculate_distance_row(n);
        nodes.push_back(n);
    }

    void add_distance_column(const Node &new_node)
    {
        for (int a = 0; a < int(nodes.size()); a++)
        {
            Node &node = nodes[a];
            if (node.removed) // todo: comment ?
                continue;

            node.distance_to.push_back(calculate_column_distance(a, new_node));
        }
    }

    double calculate_column_distance(int a, const Node &new_node) // todo: really calculate
    {
        auto a_dist = distance(a, new_node.merge.a);
        auto b_dist = distance(a, new_node.merge.b);

        return middle(a_dist, b_dist);
    }

    double middle(double a_dist, double b_dist)
    {
        if (a_dist == INF_DISTANCE || b_dist == INF_DISTANCE)
            return INF_DISTANCE;

        return (a_dist + b_dist)/2;
    }

    void calculate_distance_row(Node &node)
    {
        for (int to = 0; to < int(nodes.size()); to++)
        {
            auto a_dist = distance(node.merge.a, to);
            auto b_dist = distance(node.merge.b, to);

            node.distance_to.push_back(middle(a_dist, b_dist)); // todo: skip some
        }

        node.distance_to.push_back(INF_DISTANCE);
    }

};

void print_tree(const DistanceTable &t, int node, int indent = 0)
{
    if (node < int(files.size()))
    {
        for (int i=0; i<indent; i++)
            cout << ' ';
        cout << files[node].filename << endl;
    }

    if (t.nodes[node].merge.a >= 0)
    {
        print_tree(t, t.nodes[node].merge.a, indent + 1);
        print_tree(t, t.nodes[node].merge.b, indent + 1);
    }
}

void print(const DistanceTable &t)
{
    return;
    for (int a = 0; a < int(t.nodes.size()); a++)
    {
        if (t.nodes[a].removed)
            cout << "x  ";
        else
            cout << "   ";

        for (int b = 0; b < int(t.nodes[a].distance_to.size()); b++)
            cout << t.nodes[a].distance_to[b] << '\t';
        cout << endl;
    }
}

DistanceTo load_distance(const string &filename)
{
    DistanceTo d;

    ifstream f(filename);
	if (f.fail() || f.bad())
        return d;

    int count = 0;
    f >> count;
    d.resize(count);
    for (int i = 0; i < count; i++)
    {
        f >> d[i];
        if (!d[i])
            throw "load_distance:: !d[i]";
    }

    return d;
}

void load_table(const string &file_list_name, DistanceTable &table)
{
	FileListLoader file_list(file_list_name);
    int file_number = 0;

    files = file_list.files; // global. todo: fix
    table = DistanceTable(file_list.files.size());

    const int THREADS = 16;

	#pragma omp parallel num_threads(THREADS)
    for (int file_number = omp_get_thread_num(); file_number < int(file_list.files.size()); file_number += THREADS)
    {
        auto &file = file_list.files[file_number];
//        cerr << file_number << " of " << file_list.files.size() << " loading file " << file.filename << endl;
        auto distance = load_distance(file.filename + ".distance.13");
        table.set_distance(file_number, distance);
    }
}

MergePair find_closest(const DistanceTable &table, int to_node)
{
    MergePair best_pair(-1, -1);
   
    for (int b = 0; b < table.nodes.size(); b++)
        if (b != to_node && !table.nodes[b].removed)
        {
            auto dist = std::max(table.distance(to_node, b), table.distance(b, to_node));
            if (dist < best_pair.distance)
                best_pair = MergePair(to_node, b, dist);
        }

    return best_pair;
}

MergePair find_merge_pair(const DistanceTable &table)
{
    MergePair best_pair(-1, -1);
    for (int a = 0; a < table.nodes.size(); a++)
        if (!table.nodes[a].removed)
        {
            auto pair = find_closest(table, a);
//            cout << a << " closest " << pair.a << " " << pair.b << " " << pair.distance << endl;
            if (pair.valid() && pair.distance < best_pair.distance)
                best_pair = pair;
        }

    return best_pair;
}

void merge_table(DistanceTable &table)
{
    while (true)
    {
        auto merge_pair = find_merge_pair(table);
        if (!merge_pair.valid())
            break;

        cout << "merge: " << merge_pair.a << " " << merge_pair.b << " " << merge_pair.distance << endl;
        print_tree(table, merge_pair.a);
        print_tree(table, merge_pair.b);
        table.add_node(merge_pair);
        print(table);

        //static int count = 0;
        //count ++;
        //if (count >=2)
        //    return;
    }
}

int main(int argc, char const *argv[])
{
    try
    {
		Config config(argc, argv);

		auto before = high_resolution_clock::now();

        DistanceTable table;
        load_table(config.file_list, table);
        print(table);
        merge_table(table);

		cerr << "total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;
        return 0;
    }
    catch ( exception & x )
    {
        cerr << x.what() << endl;
		cerr << "exit 3" << endl;
		return 3;
    }
  //  catch ( string & x )
  //  {
  //      cerr << x << endl;
		//cerr << "exit 4" << endl;
		//return 4;
  //  }
    catch ( const char * x )
    {
        cerr << x << endl;
		cerr << "exit 5" << endl;
		return 5;
    }
    catch ( char const * x )
    {
        cerr << x << endl;
		cerr << "exit 5.1" << endl;
		return 5;
    }
    catch ( ... )
    {
        cerr << "unknown exception" << endl;
//		cerr << "exit 6" << endl;
		return 6;
    }
}
