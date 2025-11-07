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
    for(size_t i = 0; i < ecfp4_chunks_num; ++i)
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

    if(this->fingerprints.find(cid) == this->fingerprints.end())    
        this->retrieve_fingerprints(std::vector<cid_t>{cid});
    
    return this->fingerprints[cid];
}


double chm::App::compute_tanimoto(const ecfp4_t& a, const ecfp4_t& b)
{
    ecfp4_pc_t and_pc = 0;
    for(size_t i = 0; i < ecfp4_chunks_num; ++i)
        and_pc += std::popcount(a.first[i] & a.first[i]);
    
    ecfp4_pc_t or_pc = a.second + b.second - and_pc;

    return or_pc == 0 ? 0 : (double)and_pc / or_pc;
}


inline void chm::App::compute_max_targets_tanimoto(cid_t cid, const std::vector<ecfp4_t>& targets_fp, std::unordered_map<cid_t, double>& cids_sim)
{
    if(cids_sim.find(cid) != cids_sim.end())
        return;

    ecfp4_t fp = this->retrieve_fingerprint_single(cid);
    double max_tanimoto = 0.0;
    for(const ecfp4_t& targ_fp : targets_fp)
        max_tanimoto = std::max(this->compute_tanimoto(fp, targ_fp), max_tanimoto);
    cids_sim[cid] = max_tanimoto;
}


bool chm::App::max_targets_tanimoto_cmp(const cid_t a, const cid_t b, const std::vector<ecfp4_t>& targets_fp, std::unordered_map<cid_t, double>& cids_sim)
{
    this->compute_max_targets_tanimoto(a, targets_fp, cids_sim);
    this->compute_max_targets_tanimoto(b, targets_fp, cids_sim);
    return cids_sim[a] < cids_sim[b];
}


std::vector<std::vector<chm::App::cid_t>> chm::App::find_paths(const std::vector<cid_t>& sources, const std::vector<cid_t>& targets, uint16_t max_cost, uint16_t max_paths)
{
    if(sources.size() == 0 && targets.size() == 0)
        return std::vector<std::vector<cid_t>>();

    if(sources.size() == 0)
        return this->find_paths_targets_only(targets, max_cost, max_paths);
    if(targets.size() == 0)
        return this->find_paths_sources_only(sources, max_cost, max_paths);
    
    return this->find_paths_both_lists(sources, targets, max_cost, max_paths);
}


std::vector<std::vector<chm::App::cid_t>> chm::App::find_paths_both_lists(const std::vector<cid_t>& sources, const std::vector<cid_t>& targets, uint16_t max_cost, uint16_t max_paths)
{
    std::vector<ecfp4_t> targets_fp;
    targets_fp.reserve(targets.size());
    std::transform(targets.begin(), targets.end(), std::back_inserter(targets_fp),
        [this](cid_t cid) { return this->retrieve_fingerprint_single(cid); });
    
    std::unordered_map<cid_t, double> cids_sim;
    
    using traverse_trie_t = Trie<int, cid_t>;
    traverse_trie_t paths_trie(0);
    auto cmp =
    [this, &targets_fp, &cids_sim](const traverse_trie_t *a, const traverse_trie_t *b)
    { return this->max_targets_tanimoto_cmp(a->get_key(), b->get_key(), targets_fp, cids_sim); };
    
    std::priority_queue<traverse_trie_t*, std::vector<traverse_trie_t*>, decltype(cmp)> pqueue(cmp);
    for(cid_t cid : sources)
    {
        traverse_trie_t *src_node = paths_trie.insert_key(cid);
        pqueue.push(src_node);
    }

    auto get_visited_cids =
    [](const traverse_trie_t *trie, std::unordered_set<cid_t>& visited)
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



std::vector<std::vector<chm::App::cid_t>> chm::App::find_paths_single_list(const std::vector<cid_t>& cids, graph_t& graph, uint16_t max_cost, uint16_t max_paths)
{
    std::vector<ecfp4_t> cids_fp;
    cids_fp.reserve(cids.size());
    std::transform(cids.begin(), cids.end(), std::back_inserter(cids_fp),
        [this](cid_t cid) { return this->retrieve_fingerprint_single(cid); });
    
    std::unordered_map<cid_t, double> cids_sim;
    
    using traverse_trie_t = Trie<int, cid_t>;
    traverse_trie_t paths_trie(0);
    auto cmp =
    [this, &cids_fp, &cids_sim](const traverse_trie_t *a, const traverse_trie_t *b)
    { return this->max_targets_tanimoto_cmp(a->get_key(), b->get_key(), cids_fp, cids_sim); };
    
    std::priority_queue<traverse_trie_t*, std::vector<traverse_trie_t*>, decltype(cmp)> pqueue(cmp);
    for(cid_t cid : cids)
    {
        traverse_trie_t *src_node = paths_trie.insert_key(cid);
        pqueue.push(src_node);
    }

    auto get_visited_cids = 
    [](const traverse_trie_t *trie, std::unordered_set<cid_t>& visited)
    {
        const traverse_trie_t *curr_trie = trie;
        while(curr_trie->get_parent())
        {
            visited.insert(curr_trie->get_key());
            curr_trie = curr_trie->get_parent();
        }
    };

    std::vector<std::vector<cid_t>> result_paths;

    auto save_path_check_exit =
    [&result_paths, max_paths] (const traverse_trie_t *trie)
    {
        result_paths.push_back(trie->get_path());
        return result_paths.size() == max_paths;
    };

    while(pqueue.size())
    {
        traverse_trie_t *curr_trie = pqueue.top();
        pqueue.pop();

        cid_t cid = curr_trie->get_key();

        std::unordered_set<cid_t> visited;
        get_visited_cids(curr_trie, visited);

        for(const auto& neigh : graph[cid])
        {
            cid_t neigh_cid = neigh.first;
            if(visited.find(neigh_cid) != visited.end())
                continue;

            traverse_trie_t *neigh_trie = curr_trie->insert_key(neigh_cid);
            if(curr_trie->get_depth() == max_cost + 1U)
            {
                if(save_path_check_exit(neigh_trie))
                    goto return_paths;
            }
            else
                pqueue.push(neigh_trie);
        }

        if(curr_trie->children.size() == 0)
        {
            if(save_path_check_exit(curr_trie))
                goto return_paths;
        }
    }

    return_paths:

    return result_paths;
}


std::vector<std::vector<chm::App::cid_t>> chm::App::find_paths_sources_only(const std::vector<cid_t>& sources, uint16_t max_cost, uint16_t max_paths)
{
    return this->find_paths_single_list(sources, this->graph, max_cost, max_paths);
}

std::vector<std::vector<chm::App::cid_t>> chm::App::find_paths_targets_only(const std::vector<cid_t>& targets, uint16_t max_cost, uint16_t max_paths)
{
    return this->find_paths_single_list(targets, this->graph_reverse, max_cost, max_paths);
}



std::pair<std::unordered_set<chm::App::cid_t>, chm::App::graph_t> chm::App::convert_paths_to_graph(const std::vector<std::vector<cid_t>>& paths, bool primary_only)
{
    graph_t paths_graph;
    std::unordered_set<cid_t> secondary_cids;
    for(const auto& path : paths)
    {
        for(size_t i = 0; i < path.size()-1; ++i)
        {
            cid_t src = path[i];
            cid_t tgt = path[i + 1];

            auto src_it = this->graph.find(src);
            if(src_it == this->graph.end())
                throw std::invalid_argument("Invalid path: missing source compound");

            auto tgt_it = src_it->second.find(tgt);
            if(tgt_it == src_it->second.end())
                throw std::invalid_argument("Invalid path: missing target compound");

            paths_graph[src][tgt] = tgt_it->second;
        }
    }
    
    if(!primary_only)
    {
        graph_t secondary_graph;
        for(auto& [source, targets] : paths_graph)
        {
            for(const auto& [target, reactions] : targets)
            {
                this->retrieve_reaction_infos(reactions);

                for(const std::string& rid : reactions)
                {
                    for(cid_t src_cid : this->reaction_participants[rid].first)
                    {
                        if(paths_graph.find(src_cid) != paths_graph.end())
                            continue;
                        secondary_graph[src_cid][target].push_back(rid);
                        secondary_cids.insert(src_cid);
                    }
                }
            }
        }
        paths_graph.insert(secondary_graph.begin(), secondary_graph.end());
    }

    return std::make_pair(secondary_cids, paths_graph);
}
