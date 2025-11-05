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

        static constexpr uint16_t ecfp4_items_num = 16;
        using ecfp4_chunk_t = uint64_t;
        using ecfp4_bits_t = std::array<ecfp4_chunk_t, ecfp4_items_num>;
        using ecfp4_pc_t = uint32_t;
        using ecfp4_t = std::pair<ecfp4_bits_t, ecfp4_pc_t>;

    private:
        pqxx::connection *conn{nullptr};

        graph_t graph, graph_reverse;

        std::unordered_map<cid_t, ecfp4_t> fingerprints;

        FuzzyMap<cid_t> fuzzy;

        pqxx::result run_pqxx_request(std::string sql);
        template<typename... Args>
        pqxx::result run_pqxx_request_params(std::string sql, Args&&... args);

        void setup_graph();
        void setup_fuzzy();
        
        ecfp4_t parse_pqxx_row_fingerprint(const pqxx::row& row);
        template<typename InputIt>
        void retrieve_fingerprints(InputIt cids_begin, InputIt cids_end);
        ecfp4_t retrieve_fingerprint_single(cid_t cid);

        double compute_tanimoto(const ecfp4_t& a, const ecfp4_t& b);
        inline void compute_min_targets_tanimoto(cid_t cid, const std::vector<ecfp4_t>& targets_fp, std::unordered_map<cid_t, double>& cid_dist);
        bool min_targets_tanimoto_cmp(const cid_t a, const cid_t b, const std::vector<ecfp4_t>& targets_fp, std::unordered_map<cid_t, double>& cid_dist);

        std::vector<std::vector<cid_t>> find_paths(const std::vector<cid_t>& sources, const std::vector<cid_t>& targets, uint16_t max_cost, uint16_t max_paths);

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

    template<typename InputIt>
    void chm::App::retrieve_fingerprints(InputIt cids_begin, InputIt cids_end)
    {    
        std::ostringstream cids_stream;
        cids_stream << '{';
        for(InputIt it = cids_begin; it != cids_end; ++it)
        {
            cid_t cid = *it;
            if(this->fingerprints.find(cid) != this->fingerprints.end())
                continue;

            if(it != cids_begin)
                cids_stream << ',';

            cids_stream << cid;
        }
        cids_stream << '}';

        std::string sql_fp =
            "SELECT cid, bits, popcount "
            "FROM compound_fingerprints "
            "WHERE cid = ANY($1)";
        
        pqxx::result res = this->run_pqxx_request_params(sql_fp, (std::string)cids_stream.str());

        for(const auto& row : res)
        {
            cid_t cid = row[0].as<cid_t>();
            this->fingerprints[cid] = this->parse_pqxx_row_fingerprint(row);
        }
    }
}

#endif // CHEMS_GRAPH_HPP