#ifndef CHEMS_GRAPH_HPP
#define CHEMS_GRAPH_HPP

#include <pqxx/pqxx>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <array>
#include <fmt/format.h>
#include <nlohmann/json.hpp>

#include "fuzzy_map.hpp"

namespace chm
{
    class App
    {
    public:
        static constexpr double epsilon = 1e-6;
        
        using cid_t = int32_t;
        using graph_t = std::unordered_map<cid_t, std::unordered_map<cid_t, std::vector<std::string>>>;

        static constexpr uint16_t ecfp4_chunks_num = 16;
        using ecfp4_chunk_t = uint64_t;
        using ecfp4_bits_t = std::array<ecfp4_chunk_t, ecfp4_chunks_num>;
        using ecfp4_pc_t = uint32_t;
        using ecfp4_t = std::pair<ecfp4_bits_t, ecfp4_pc_t>;
        static constexpr uint32_t fp_bits_num = ecfp4_chunks_num*sizeof(ecfp4_chunk_t)*8;

    private:
        pqxx::connection *conn{nullptr};

        graph_t graph, graph_reverse;

        std::unordered_map<cid_t, ecfp4_t> fingerprints;
        std::unordered_map<cid_t, nlohmann::json> compound_infos;

        FuzzyMap<cid_t> fuzzy;

        pqxx::result run_pqxx_request(std::string sql);
        template<typename... Args>
        pqxx::result run_pqxx_request_params(std::string sql, Args&&... args);

        void setup_graph();
        void setup_fuzzy();


        template<typename Iterable, typename CidsBlacklist>
        std::string cids_to_sql_list(Iterable cids, const CidsBlacklist& target_map);
        
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


        template<typename Iterable>
        void retrieve_compound_infos(Iterable cids);
        nlohmann::json retrieve_compound_info_single(cid_t cid);

    public:
        App(std::string db_name, std::string user="", std::string password="");
        ~App();
    };


    template<typename... Args>
    pqxx::result App::run_pqxx_request_params(std::string sql, Args&&... args)
    {
        pqxx::work tx(*this->conn);
        return tx.exec_params(sql, std::forward<Args>(args)...);
    }


    template<typename Iterable, typename CidsBlacklist>
    std::string chm::App::cids_to_sql_list(Iterable cids, const CidsBlacklist& blacklist)
    {
        std::ostringstream cids_stream;
        cids_stream << '{';
        size_t added = 0;
        for(cid_t cid : cids)
        {
            if(blacklist.find(cid) != blacklist.end())
                continue;

            if(added != 0)
                cids_stream << ',';
            cids_stream << cid;
            ++added;
        }
        cids_stream << '}';

        return cids_stream.str();
    }

    template<typename Iterable>
    void chm::App::retrieve_fingerprints(Iterable cids)
    {    
        std::string cids_str = this->cids_to_sql_list(cids, this->fingerprints);

        std::string sql_fp =
            "SELECT cid, bits, popcount "
            "FROM compound_fingerprints "
            "WHERE cid = ANY($1)";
        
        pqxx::result res = this->run_pqxx_request_params(sql_fp, cids_str);

        for(const auto& row : res)
        {
            cid_t cid = row[0].as<cid_t>();
            this->fingerprints[cid] = this->parse_pqxx_row_fingerprint(row);
        }
    }


    template<typename Iterable>
    void chm::App::retrieve_compound_infos(Iterable cids)
    {
        std::string cids_str = this->cids_to_sql_list(cids, this->compound_infos);

        std::unordered_map<cid_t, nlohmann::json> compound_infos;
        
        std::string sql_compound_properties =
            "SELECT cid, property_name, property_value, rank "
            "FROM compound_properties "
            "WHERE cid = ANY($1)";
        
        pqxx::result prop_rows{this->run_pqxx_request_params(sql_compound_properties, cids_str)};
        for(const auto& row : prop_rows)
        {
            cid_t cid = row[0].as<cid_t>();
            std::string prop_name = row[1].as<std::string>();
            std::string prop_value = row[2].as<std::string>();
            compound_infos[cid]["properties"].push_back({{"property_name", prop_name}, {"property_value", prop_value}});
        }
        
        std::string sql_compound_pictograms =
            "SELECT cid, pictogram "
            "FROM compound_hazard_pictograms "
            "WHERE cid = ANY($1)";
        
        pqxx::result picts_rows{this->run_pqxx_request_params(sql_compound_pictograms, cids_str)};
        for(const auto& row : picts_rows)
        {
            cid_t cid = row[0].as<cid_t>();
            compound_infos[cid]["hazard_pictograms"].push_back(row[1].as<std::string>());
        }
        
        std::string sql_compound_nfpa =
            "SELECT cid, health, flammability, instability "
            "FROM compound_nfpa "
            "WHERE cid = ANY($1)";
        
        pqxx::result nfpa_rows{this->run_pqxx_request_params(sql_compound_nfpa, cids_str)};
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

        for(auto& entry : compound_infos)
            this->compound_infos[entry.first] = std::move(entry.second);
    }
}

#endif // CHEMS_GRAPH_HPP