#ifndef CHEMS_GRAPH_HPP
#define CHEMS_GRAPH_HPP

#include <pqxx/pqxx>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <array>
#include <filesystem>
#include <fmt/format.h>
#include <fmt/color.h>
#include <fmt/core.h>
#include <nlohmann/json.hpp>
#include <mutex>

#include "fuzzy_map.hpp"
#include "misc.hpp"


namespace fs = std::filesystem;


namespace chm
{
    class App
    {
    public:
        static constexpr double epsilon = 1e-6;
        
        using cid_t = int32_t;

        static constexpr uint16_t ecfp4_chunks_num = 16;
        using ecfp4_chunk_t = uint64_t;
        using ecfp4_bits_t = std::array<ecfp4_chunk_t, ecfp4_chunks_num>;
        using ecfp4_pc_t = uint32_t;
        using ecfp4_t = std::pair<ecfp4_bits_t, ecfp4_pc_t>;
        static constexpr uint32_t fp_bits_num = ecfp4_chunks_num*sizeof(ecfp4_chunk_t)*8;

        inline static const fs::path DB_INFO_REL_PATH = "db.json";

    protected:
        fs::path DATA_DIR_PATH;
        fs::path DB_INFO_FULL_PATH;
        fs::path UI_DIR_PATH, STRUCTURES_DIR_UI_PATH;

        uint32_t PAGE_SIZE = 100;

        using sorting_map_t = std::unordered_map<cid_t, uint32_t>;
        sorting_map_t complexity_sorting, commonness_sorting, curiosity_sorting;
        std::vector<cid_t> complexity_sorted_cids, commonness_sorted_cids, curiosity_sorted_cids;
        std::unordered_map<std::string, std::pair<sorting_map_t, std::vector<cid_t>>> sorting;

        using graph_t = std::unordered_map<cid_t, std::unordered_map<cid_t, std::vector<std::string>>>;
        graph_t graph, graph_reverse;

        nlohmann::json convert_paths_to_graph(const std::vector<std::vector<cid_t>>& paths, bool primary_only);
    
    private:
        pqxx::connection *conn{nullptr};
        std::mutex conn_mutex;

        using fingerprints_t = std::unordered_map<cid_t, ecfp4_t>;
        fingerprints_t fingerprints;

        using compound_infos_t = std::unordered_map<cid_t, nlohmann::json>;
        compound_infos_t compound_infos;

        using reaction_infos_t = std::unordered_map<std::string, nlohmann::json>;
        reaction_infos_t reaction_infos;

        using reaction_participants_t = std::unordered_map<std::string, std::pair<std::vector<cid_t>, std::vector<cid_t>>>;
        reaction_participants_t reaction_participants;

        FuzzyMap<cid_t> fuzzy;

        pqxx::result run_pqxx_request(const std::string& sql);
        template<typename... Args>
        pqxx::result run_pqxx_request_params(const std::string& sql, Args&&... args);
        template<typename... Args>
        pqxx::result run_pqxx_request_params_prepared(const std::string& prepared, Args&&... args);

        void setup_sql();
        void setup_sortings();
        void setup_graph();
        void setup_fuzzy();

        
        ecfp4_t parse_pqxx_row_fingerprint(const pqxx::row& row);
        template<typename Iterable>
        void retrieve_fingerprints(Iterable cids);
        ecfp4_t retrieve_fingerprint_single(cid_t cid);

        double compute_tanimoto(const ecfp4_t& a, const ecfp4_t& b);
        inline void compute_max_targets_tanimoto(cid_t cid, const std::vector<ecfp4_t>& targets_fp, std::unordered_map<cid_t, double>& cid_dist);
        bool max_targets_tanimoto_cmp(const cid_t a, const cid_t b, const std::vector<ecfp4_t>& targets_fp, std::unordered_map<cid_t, double>& cid_dist);
        
        std::vector<std::vector<cid_t>> find_paths(const std::vector<cid_t>& sources, const std::vector<cid_t>& targets, uint16_t max_cost, uint16_t max_paths);
        std::vector<std::vector<cid_t>> find_paths_both_lists(const std::vector<cid_t>& sources, const std::vector<cid_t>& targets, uint16_t max_cost, uint16_t max_paths);
        std::vector<std::vector<cid_t>> find_paths_single_list(const std::vector<cid_t>& cids, graph_t& graph, uint16_t max_cost, uint16_t max_paths);
        std::vector<std::vector<cid_t>> find_paths_sources_only(const std::vector<cid_t>& sources, uint16_t max_cost, uint16_t max_paths);
        std::vector<std::vector<cid_t>> find_paths_targets_only(const std::vector<cid_t>& targets, uint16_t max_cost, uint16_t max_paths);

    public:
        App(const fs::path& data_dir);
        ~App();


        template<typename Iterable>
        void retrieve_compound_infos(Iterable cids);
        template<typename Iterable>
        nlohmann::json retrieve_compound_infos_json(Iterable cids);
        nlohmann::json retrieve_compound_info_single(cid_t cid);

        template<typename Iterable>
        void retrieve_reaction_infos(Iterable rids);
        template<typename Iterable>
        nlohmann::json retrieve_reaction_infos_json(Iterable rids);
        nlohmann::json retrieve_reaction_info_single(const std::string& rid);

        nlohmann::json build_graph(const std::vector<cid_t>& sources, const std::vector<cid_t>& targets, uint16_t max_cost, uint16_t max_paths, bool primary_only);

        std::vector<cid_t> search_compounds(const std::string& query);

        void generate_compound_structure_svg(cid_t cid);
    };


    template<typename... Args>
    pqxx::result App::run_pqxx_request_params(const std::string& sql, Args&&... args)
    {
        std::lock_guard<std::mutex> lock(this->conn_mutex);
        pqxx::work tx(*this->conn);
        pqxx::result result = tx.exec_params(sql, std::forward<Args>(args)...);
        tx.commit();
        return result;
    }


    template<typename... Args>
    pqxx::result App::run_pqxx_request_params_prepared(const std::string& prepared, Args&&... args)
    {
        std::lock_guard<std::mutex> lock(this->conn_mutex);
        pqxx::work tx(*this->conn);
        pqxx::result result = tx.exec_prepared(prepared, std::forward<Args>(args)...);;
        tx.commit();
        return result;
    }


    template<typename Iterable>
    void chm::App::retrieve_fingerprints(Iterable cids)
    {    
        std::string cids_str = entries_to_sql_list(cids, this->fingerprints);
        pqxx::result res{this->run_pqxx_request_params_prepared("compound_fingerprints", cids_str)};

        for(const auto& row : res)
        {
            cid_t cid = row[0].as<cid_t>();
            this->fingerprints[cid] = this->parse_pqxx_row_fingerprint(row);
        }
    }


    template<typename Iterable>
    void chm::App::retrieve_compound_infos(Iterable cids)
    {
        std::string cids_str = entries_to_sql_list(cids, this->compound_infos);

        compound_infos_t compound_infos;

        pqxx::result cmpd_rows{this->run_pqxx_request_params_prepared("compounds", cids_str)};
        for(const auto& row : cmpd_rows)
        {
            cid_t cid = row[0].as<cid_t>();
            compound_infos[cid]["cid"] = cid;
            compound_infos[cid]["name"] = row[1].as<std::string>();
            compound_infos[cid]["smiles"] = row[2].as<std::string>();
            compound_infos[cid]["organic"] = row[3].as<bool>();
            compound_infos[cid]["wiki"] = row[4].is_null() ? "" : row[4].as<std::string>();
        }

        pqxx::result prop_rows{this->run_pqxx_request_params_prepared("compound_properties", cids_str)};
        for(const auto& row : prop_rows)
        {
            cid_t cid = row[0].as<cid_t>();
            std::string prop_name = row[1].as<std::string>();
            std::string prop_value = row[2].as<std::string>();
            compound_infos[cid]["properties"].push_back({{"property_name", prop_name}, {"property_value", prop_value}});
        }

        pqxx::result picts_rows{this->run_pqxx_request_params_prepared("compound_pictograms", cids_str)};
        for(const auto& row : picts_rows)
        {
            cid_t cid = row[0].as<cid_t>();
            compound_infos[cid]["pictograms"].push_back(row[1].as<std::string>());
        }
        
        pqxx::result nfpa_rows{this->run_pqxx_request_params_prepared("compound_nfpa", cids_str)};
        nlohmann::json nfpa;
        for(const auto& row : nfpa_rows)
        {
            cid_t cid = row[0].as<cid_t>();
            if(!row[1].is_null())
                compound_infos[cid]["nfpa"]["health"] = row[1].as<std::string>();
            if(!row[2].is_null())
                compound_infos[cid]["nfpa"]["flammability"] = row[2].as<std::string>();
            if(!row[3].is_null())
                compound_infos[cid]["nfpa"]["instability"] = row[3].as<std::string>();
        }

        pqxx::result descr_rows{this->run_pqxx_request_params_prepared("compound_descriptions", cids_str)};
        for(const auto& row : descr_rows)
        {
            cid_t cid = row[0].as<cid_t>();
            compound_infos[cid]["description"] = row[1].as<std::string>();
        }

        for(auto& [cid, info] : compound_infos)
            this->compound_infos[cid] = std::move(info);
    }


    template<typename Iterable>
    nlohmann::json chm::App::retrieve_compound_infos_json(Iterable cids)
    {
        this->retrieve_compound_infos(cids);

        nlohmann::json compound_infos;
        for(cid_t cid : cids)
            compound_infos.push_back(this->compound_infos[cid]);
        
        return compound_infos;
    }


    template<typename Iterable>
    void chm::App::retrieve_reaction_infos(Iterable rids)
    {
        std::string rids_str = entries_to_sql_list(rids, this->reaction_infos);

        reaction_infos_t reaction_infos;
        reaction_participants_t reaction_participants;

        pqxx::result details_rows{this->run_pqxx_request_params_prepared("reaction_details", rids_str)};
        for(const auto& row : details_rows)
        {
            std::string rid = row[0].as<std::string>();
            reaction_infos[rid]["balanced"] = row[1].as<bool>();
            reaction_infos[rid]["complexity"] = row[2].as<float>();
            reaction_infos[rid]["source"] = row[3].as<std::string>();
            if(!row[4].is_null())
                reaction_infos[rid]["description"] = row[4].as<std::string>();
            if(!row[5].is_null())
                reaction_infos[rid]["patent"] = row[5].as<std::string>();
            if(!row[6].is_null())
                reaction_infos[rid]["doi"] = row[6].as<std::string>();   
        }


        pqxx::result reactants_rows{this->run_pqxx_request_params_prepared("reaction_reactants", rids_str)};
        for (const auto& row : reactants_rows)
        {
            std::string rid = row[0].as<std::string>();
            nlohmann::json entry;
            cid_t cid = row[1].as<cid_t>();
            entry["cid"] = cid;
            if (!row[2].is_null())
                entry["coeff"] = row[2].as<int16_t>();
            else
                entry["coeff"] = nullptr;
            reaction_infos[rid]["reactants"].push_back(entry);
            reaction_participants[rid].first.push_back(cid);
        }

        pqxx::result products_rows{this->run_pqxx_request_params_prepared("reaction_products", rids_str)};
        for (const auto& row : products_rows)
        {
            std::string rid = row[0].as<std::string>();
            nlohmann::json entry;
            cid_t cid = row[1].as<cid_t>();
            entry["cid"] = cid;
            if (!row[2].is_null())
                entry["coeff"] = row[2].as<int16_t>();
            else
                entry["coeff"] = nullptr;
            reaction_infos[rid]["products"].push_back(entry);
            reaction_participants[rid].second.push_back(cid);
        }

        for(auto& [rid, info] : reaction_infos)
            this->reaction_infos[rid] = std::move(info);
        for(auto& [rid, participants] : reaction_participants)
            this->reaction_participants[rid] = std::move(participants);
    }

    template<typename Iterable>
    nlohmann::json chm::App::retrieve_reaction_infos_json(Iterable rids)
    {
        this->retrieve_reaction_infos(rids);

        nlohmann::json reaction_infos;
        for(const std::string& rid : rids)
            reaction_infos.push_back(this->reaction_infos[rid]);
        
        return reaction_infos;
    }
}

#endif // CHEMS_GRAPH_HPP