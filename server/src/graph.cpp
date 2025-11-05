#include <queue>
#include <bit>
#include <algorithm>
#include <execution>
#include <tuple>

#include "app.hpp"
#include "misc.hpp"


chm::App::ecfp4_t chm::App::parse_pqxx_row_fingerprint(const pqxx::row& row)
{
    pqxx::array_parser fp_parser = row[1].as_array();
    fp_parser.get_next();

    ecfp4_bits_t fp_bits;
    for(size_t i = 0; i < ecfp4_items_num; ++i)
    {
        auto parsed = fp_parser.get_next();
        if(parsed.first != fp_parser.string_value)
        {
            cid_t cid = row[0].as<cid_t>();
            std::string fp_str = row[1].as<std::string>();
            throw std::invalid_argument(fmt::format("(juncture={}) Invalid fingerprint for cid {}: '{}'", parsed.first, cid, fp_str));
        }
        
        int64_t fp_chunk = str_to_numeric<int64_t>(parsed.second);
        fp_bits[i] = static_cast<ecfp4_chunk_t>(fp_chunk);
    }

    ecfp4_pc_t popcount = row[2].as<ecfp4_pc_t>();

    return std::make_pair(std::move(fp_bits), popcount);
}


chm::App::ecfp4_t chm::App::retrieve_fingerprint_single(cid_t cid)
{

    if(this->fingerprints.find(cid) != this->fingerprints.end())
        return this->fingerprints[cid];

    std::string sql_fp =
        "SELECT cid, bits, popcount "
        "FROM compound_fingerprints "
        "WHERE cid = " + numeric_to_str(cid);
    
    pqxx::result res = this->run_pqxx_request(sql_fp);

    if(res.size() != 1)
        std::invalid_argument(fmt::format("Invalid number of entries retrieved for cid {}. Got {}, expected 1", cid, res.size()));

    pqxx::row row = *res.begin();

    ecfp4_t fp = this->parse_pqxx_row_fingerprint(row);
    this->fingerprints[cid] = fp;

    return fp;
}


double chm::App::compute_tanimoto(const ecfp4_t& a, const ecfp4_t& b)
{
    ecfp4_pc_t and_pc = 0;
    for(size_t i = 0; i < ecfp4_items_num; ++i)
        and_pc += std::popcount(a.first[i] & a.first[i]);
    
    ecfp4_pc_t or_pc = a.second + b.second - and_pc;

    return or_pc == 0 ? 0 : (double)and_pc / or_pc;
}


inline void chm::App::compute_min_targets_tanimoto(cid_t cid, const std::vector<ecfp4_t>& targets_fp, std::unordered_map<cid_t, double>& cids_sim)
{
    if(cids_sim.find(cid) != cids_sim.end())
        return;

    ecfp4_t fp = this->retrieve_fingerprint_single(cid);
    std::vector<double> source_tanimoto;
    source_tanimoto.reserve(targets_fp.size());
    std::transform(targets_fp.begin(), targets_fp.end(), std::back_inserter(source_tanimoto),
        [this, &fp](ecfp4_t src_fp) { return this->compute_tanimoto(fp, src_fp); });
    cids_sim[cid] = *std::min_element(source_tanimoto.begin(), source_tanimoto.end());
}


bool chm::App::min_targets_tanimoto_cmp(const cid_t a, const cid_t b, const std::vector<ecfp4_t>& targets_fp, std::unordered_map<cid_t, double>& cids_sim)
{
    this->compute_min_targets_tanimoto(a, targets_fp, cids_sim);
    this->compute_min_targets_tanimoto(b, targets_fp, cids_sim);
    return cids_sim[a] > cids_sim[b];
}


std::vector<std::vector<chm::App::cid_t>> chm::App::find_paths(const std::vector<cid_t>& sources, const std::vector<cid_t>& targets, uint16_t max_cost, uint16_t max_paths)
{
    std::vector<ecfp4_t> targets_fp;
    targets_fp.reserve(targets.size());
    std::transform(targets.begin(), targets.end(), std::back_inserter(targets_fp),
        [this](cid_t cid) { return this->retrieve_fingerprint_single(cid); });
    
    std::unordered_map<cid_t, double> cids_sim;
    
    using traverse_trie_t = Trie<int, cid_t>;
    traverse_trie_t paths_trie(0);
    auto cmp = [this, &targets_fp, &cids_sim](const traverse_trie_t *a, const traverse_trie_t *b)
        { return this->min_targets_tanimoto_cmp(a->get_key(), b->get_key(), targets_fp, cids_sim); };
    
    std::priority_queue<traverse_trie_t*, std::vector<traverse_trie_t*>, decltype(cmp)> pqueue(cmp);
    for(cid_t cid : sources)
    {
        traverse_trie_t *src_node = paths_trie.insert_key(cid);
        pqueue.push(src_node);
    }

    auto get_visited_cids = [](const traverse_trie_t *trie, std::unordered_set<cid_t>& visited)
    {
        const traverse_trie_t *curr_trie = trie;
        while(curr_trie->get_parent())
        {
            visited.insert(curr_trie->get_key());
            curr_trie = curr_trie->get_parent();
        }
    };
    
    std::unordered_set<cid_t> targets_set(targets.begin(), targets.end());

    std::vector<std::vector<cid_t>> result_paths;
    while(pqueue.size())
    {
        traverse_trie_t *curr_trie = pqueue.top();
        pqueue.pop();

        cid_t cid = curr_trie->get_key();

        std::unordered_set<cid_t> visited;
        get_visited_cids(curr_trie, visited);

        for(const auto& neigh : this->graph[cid])
        {
            cid_t neigh_cid = neigh.first;
            if(visited.find(neigh_cid) != visited.end())
                continue;

            traverse_trie_t *neigh_trie = curr_trie->insert_key(neigh_cid);
            if(targets_set.find(neigh_cid) != targets_set.end())
            {
                result_paths.push_back(neigh_trie->get_path());
                if(result_paths.size() == max_paths)
                    goto return_paths;
            }
            else if(curr_trie->get_depth() < max_cost)
                pqueue.push(neigh_trie);
        }
    }

    return_paths:

    return result_paths;
}
