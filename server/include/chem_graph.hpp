#ifndef CHGRAPH_CHEM_GRAPH_HPP
#define CHGRAPH_CHEM_GRAPH_HPP

#include <pqxx/pqxx>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>


class ChemGraph
{
public:
    using graph_t = std::unordered_map<int32_t, std::vector<int32_t>>;

private:
    pqxx::connection *conn{nullptr};

    graph_t graph;

    void setup_graph();

public:
    ChemGraph(std::string db_name, std::string user="", std::string password="");
    ~ChemGraph();
};

#endif // CHGRAPH_CHEM_GRAPH_HPP