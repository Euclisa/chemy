#include "app.hpp"
#include <iostream>
#include <unordered_set>


chm::App::App(std::string db_name, std::string user, std::string password)
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
    this->setup_fuzzy();
}


chm::App::~App()
{
    if(this->conn)
        delete this->conn;
}


void chm::App::setup_graph()
{
    pqxx::work tx(*this->conn);

    std::string sql_pairs = 
        "SELECT r.cid, p.cid "
        "FROM reaction_reactants r "
        "JOIN reaction_products p "
        "ON r.rid = p.rid";
    
    pqxx::result pair_rows = tx.exec(sql_pairs);

    std::unordered_set<cid_t> unique_cids;
    for(const auto& row : pair_rows)
    {
        cid_t r_cid = row[0].as<cid_t>();
        cid_t p_cid = row[1].as<cid_t>();

        unique_cids.insert({r_cid, p_cid});

        graph[r_cid].push_back(p_cid);
    }

    std::cout << graph.size() << ' ' << unique_cids.size() << '\n';
}


void chm::App::setup_fuzzy()
{
    pqxx::work tx(*this->conn);

    std::string sql_synonyms =
        "SELECT synonym, cid "
        "FROM compound_synonyms";
    
    pqxx::result synonym_rows = tx.exec(sql_synonyms);
    std::vector<std::pair<std::string, cid_t>> syn_entries;
    for(const auto& row : synonym_rows)
    {
        std::string synonym = row[0].as<std::string>();
        cid_t cid = row[1].as<cid_t>();
        syn_entries.push_back(std::make_pair(synonym, cid));
    }

    this->fuzzy.insert_entries(syn_entries);
    std::cout << "Inserted\n";

    auto res = this->fuzzy.search("amphetamine");
    for(const auto& entry : res)
    {
        std::cout << entry.first << " -> " << entry.second << '\n';
    }
    std::cout << res.size() << '\n';
}


