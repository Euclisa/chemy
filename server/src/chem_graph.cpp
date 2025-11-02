#include "chem_graph.hpp"
#include <iostream>
#include <unordered_set>


ChemGraph::ChemGraph(std::string db_name, std::string user, std::string password)
{
    std::string pq_params = "dbname=" + db_name;
    if(user.size() != 0)
    {
        pq_params += " user=" + user;
        pq_params += " password=" + password;
    }
    
    try
    {
        this->conn = new pqxx::connection(pq_params);
    }
    catch(const pqxx::broken_connection& e)
    {
        std::cerr << "Failed to connect to database:\n" << e.what() << '\n';
        throw e;
    }

    this->setup_graph();
}


ChemGraph::~ChemGraph()
{
    if(this->conn)
        delete this->conn;
}


void ChemGraph::setup_graph()
{
    pqxx::work tx(*this->conn);

    std::string sql_pairs = 
        "SELECT r.cid, p.cid "
        "FROM reaction_reactants r "
        "JOIN reaction_products p "
        "ON r.rid = p.rid";
    
    pqxx::result pair_rows = tx.exec(sql_pairs);

    std::unordered_set<int32_t> unique_cids;
    for(const auto& row : pair_rows)
    {
        int32_t r_cid = row[0].as<int32_t>();
        int32_t p_cid = row[1].as<int32_t>();

        unique_cids.insert({r_cid, p_cid});

        graph[r_cid].push_back(p_cid);
    }

    std::cout << graph.size() << ' ' << unique_cids.size() << '\n';
}



