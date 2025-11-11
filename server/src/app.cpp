#include "app.hpp"

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <cstdarg>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>



chm::App::App(const fs::path& data_dir)
{
    this->DATA_DIR_PATH = data_dir;

    this->DB_INFO_FULL_PATH = this->DATA_DIR_PATH / this->DB_INFO_REL_PATH;
    if(!fs::exists(this->DB_INFO_FULL_PATH))
        throw std::invalid_argument(
            fmt::format("Failed to find db info file in provided data directory: {}",
                (std::string)this->DB_INFO_FULL_PATH));

    std::ifstream dbfin(this->DB_INFO_FULL_PATH);
    if(!dbfin.is_open())
        throw std::invalid_argument(
            fmt::format("Failed to open db info file: {}",
                (std::string)this->DB_INFO_FULL_PATH));
    
    std::stringstream db_info_buffer;
    db_info_buffer << dbfin.rdbuf();
    std::string db_info_str = db_info_buffer.str(); 
    
    nlohmann::json db_info = nlohmann::json::parse(db_info_str);

    this->UI_DIR_PATH = this->DATA_DIR_PATH / "ui";
    this->STRUCTURES_DIR_UI_PATH = "assets/structures";

    std::string pq_params = "dbname=" + (std::string)db_info["dbname"];
    if(db_info.contains("user"))
    {
        pq_params += " user=" + (std::string)db_info["user"];
        
        std::string password = db_info.contains("password") ? (std::string)db_info["password"] : std::string("");
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
    std::lock_guard<std::mutex> lock(this->conn_mutex);
    pqxx::work tx(*this->conn);
    pqxx::result result = tx.exec(sql);
    tx.commit();
    return result;
}


void chm::App::setup_sql()
{
    std::string sql_fp =
        "SELECT cid, bits, popcount "
        "FROM compound_fingerprints "
        "WHERE cid = ANY($1)";
    this->conn->prepare("compound_fingerprints", sql_fp);

    std::string sql_compounds =
        "SELECT c.cid, c.name, c.smiles, w.wiki "
        "FROM compounds c "
        "LEFT JOIN compound_wiki w "
        "ON c.cid = w.cid "
        "WHERE c.cid = ANY($1)";
    this->conn->prepare("compounds", sql_compounds);

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

    std::string sql_compound_descriptions =
        "SELECT cid, description "
        "FROM compound_descriptions "
        "WHERE cid = ANY($1)";
    this->conn->prepare("compound_descriptions", sql_compound_descriptions);

    std::string sql_reaction_details =
        "SELECT r.rid, r.balanced, r.complexity, r.source, d.description, d.patent, d.doi "
        "FROM reactions r "
        "LEFT JOIN reaction_details d "
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
    auto populate_sorting = [this](std::string table, std::string key)
    {
        std::string sql =
        "SELECT cid, rank "
        "FROM " + table;

        pqxx::result rows = this->run_pqxx_request(sql);
        sorting_map_t *sorting_map = &this->sorting[key].first;
        std::vector<cid_t> *sorted_cids = &this->sorting[key].second;
        sorted_cids->resize(rows.size());
        for(const auto& row : rows)
        {
            cid_t cid = row[0].as<cid_t>();
            uint32_t rank = row[1].as<uint32_t>();
            (*sorting_map)[cid] = rank;
            (*sorted_cids)[rank] = cid;
        }
    };

    populate_sorting("compound_complexity_sorting", "complexity");
    populate_sorting("compound_commonness_sorting", "commonness");
    populate_sorting("compound_curiosity_sorting", "curiosity");

    this->sorting["none"].second = this->sorting["commonness"].second;
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



void chm::App::generate_compound_structure_svg(cid_t cid)
{
    fs::path svg_path = this->UI_DIR_PATH / this->STRUCTURES_DIR_UI_PATH / fmt::format("{}.svg", cid);
    if(fs::exists(svg_path))
        return;

    std::string smiles = this->compound_infos[cid]["smiles"];
    RDKit::RWMol* mol = RDKit::SmilesToMol(smiles);

    if (!mol) return;

    RDDepict::compute2DCoords(*mol);

    RDKit::MolDraw2DSVG drawer(300, 200);
    drawer.drawMolecule(*mol);
    drawer.finishDrawing();

    std::ofstream svg_fout(svg_path);
    if(!svg_fout.is_open())
        throw std::runtime_error(fmt::format("Failed to open '{}' to save structure", svg_path.string()));

    svg_fout << drawer.getDrawingText();

    svg_fout.close();

    delete mol;
}