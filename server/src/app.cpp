#include "app.hpp"
#include <iostream>
#include <unordered_set>
#include <cstdarg>


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


pqxx::result chm::App::run_pqxx_request(std::string sql)
{
    pqxx::work tx(*this->conn);
    return tx.exec(sql);
}


void chm::App::setup_graph()
{
    std::string sql_pairs = 
        "SELECT r.cid, p.cid, r.rid "
        "FROM reaction_reactants r "
        "JOIN reaction_products p "
        "ON r.rid = p.rid";
    
    pqxx::result pair_rows = this->run_pqxx_request(sql_pairs);

    std::unordered_set<cid_t> unique_cids;
    for(const auto& row : pair_rows)
    {
        cid_t r_cid = row[0].as<cid_t>();
        cid_t p_cid = row[1].as<cid_t>();
        std::string rid = row[2].as<std::string>();

        unique_cids.insert({r_cid, p_cid});

        this->graph[r_cid][p_cid].push_back(rid);
        this->graph_reverse[p_cid][r_cid].push_back(rid);
    }

    this->retrieve_fingerprints(unique_cids);
    std::cout << this->compound_infos[222].dump(4) << '\n';
    std::cout << "Starting" << '\n';
    auto paths = this->find_paths({222}, {313}, 4, 100);
    for(const auto& path : paths)
    {
        for(cid_t node : path)
            std::cout << node << ' ';
        std::cout << '\n';
    }

    std::cout << this->graph.size() << ' ' << unique_cids.size() << "\n\n";
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


