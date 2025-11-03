#ifndef CHEMS_GRAPH_HPP
#define CHEMS_GRAPH_HPP

#include <pqxx/pqxx>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

#include "fuzzy_map.hpp"

namespace chm
{
    class App
    {
    public:
        using cid_t = int32_t;
        using graph_t = std::unordered_map<cid_t, std::vector<cid_t>>;

    private:
        pqxx::connection *conn{nullptr};

        graph_t graph;

        FuzzyMap<cid_t> fuzzy;

        void setup_graph();
        void setup_fuzzy();

    public:
        App(std::string db_name, std::string user="", std::string password="");
        ~App();
    };
}

#endif // CHEMS_GRAPH_HPP