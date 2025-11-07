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

    this->setup_sql();
    this->setup_sortings();
    this->setup_graph();
    this->setup_fuzzy();
}


chm::App::~App()
{
    if(this->conn)
        delete this->conn;
}


pqxx::result chm::App::run_pqxx_request(const std::string& sql)
{
    pqxx::work tx(*this->conn);
    return tx.exec(sql);
}


void chm::App::setup_sql()
{
    std::string sql_fp =
        "SELECT cid, bits, popcount "
        "FROM compound_fingerprints "
        "WHERE cid = ANY($1)";
    this->conn->prepare("compound_fingerprints", sql_fp);

    std::string sql_compound_properties =
        "SELECT cid, property_name, property_value, rank "
        "FROM compound_properties "
        "WHERE cid = ANY($1)";
    this->conn->prepare("compound_properties", sql_compound_properties);

    std::string sql_compound_pictograms =
        "SELECT cid, pictogram "
        "FROM compound_hazard_pictograms "
        "WHERE cid = ANY($1)";
    this->conn->prepare("compound_pictograms", sql_compound_pictograms);

    std::string sql_compound_nfpa =
        "SELECT cid, health, flammability, instability "
        "FROM compound_nfpa "
        "WHERE cid = ANY($1)";
    this->conn->prepare("compound_nfpa", sql_compound_nfpa);

    std::string sql_reaction_details =
        "SELECT r.rid, r.balanced, r.complexity, r.source, d.description, d.patent, d.doi "
        "FROM reactions r "
        "JOIN reaction_details d "
        "ON r.rid = d.rid "
        "WHERE r.rid = ANY($1)";
    this->conn->prepare("reaction_details", sql_reaction_details);

    std::string sql_reaction_reactants =
        "SELECT rid, cid, coeff "
        "FROM reaction_reactants "
        "WHERE rid = ANY($1)";
    this->conn->prepare("reaction_reactants", sql_reaction_reactants);

    std::string sql_reaction_products =
        "SELECT rid, cid, coeff "
        "FROM reaction_products "
        "WHERE rid = ANY($1)";
    this->conn->prepare("reaction_products", sql_reaction_products);
}


void chm::App::setup_sortings()
{
    auto populate_sorting = [&](std::string table, sorting_t& sorting)
    {
        std::string sql =
        "SELECT cid, rank "
        "FROM " + table;

        pqxx::result rows = this->run_pqxx_request(sql);
        for(const auto& row : rows)
        {
            cid_t cid = row[0].as<cid_t>();
            uint32_t rank = row[1].as<uint32_t>();
            sorting[cid] = rank;
        }
    };

    populate_sorting("compound_complexity_sorting", this->complexity_sorting);
    populate_sorting("compound_commonness_sorting", this->commonness_sorting);
    populate_sorting("compound_curiosity_sorting", this->curiosity_sorting);
}


void chm::App::setup_graph()
{
    std::string sql_pairs = 
        "SELECT r.cid, p.cid, r.rid "
        "FROM reaction_reactants r "
        "LEFT JOIN reaction_products p "
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
    std::cout << this->retrieve_compound_info_single(222).dump(4) << '\n';
    this->retrieve_reaction_info_single("+hcR5qo9GYZLBXiylUkB+Q==");
    std::cout << this->reaction_infos["+hcR5qo9GYZLBXiylUkB+Q=="].dump(4) << '\n';
    std::cout << "Starting" << '\n';
    std::cout << this->build_graph({222}, {313}, 4, 10, false).dump(4) << '\n';

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


