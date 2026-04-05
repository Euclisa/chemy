#include "app.hpp"
#include "misc.hpp"


nlohmann::json chm::App::retrieve_compound_info_single(cid_t cid)
{
    {
        std::lock_guard<std::mutex> lock(this->compound_infos_mutex);
        auto info_it = this->compound_infos.find(cid);
        if(info_it != this->compound_infos.end())
            return info_it->second;
    }

    this->retrieve_compound_infos(std::vector<cid_t>{cid});

    std::lock_guard<std::mutex> lock(this->compound_infos_mutex);
    auto info_it = this->compound_infos.find(cid);
    if(info_it == this->compound_infos.end())
        throw std::invalid_argument(fmt::format("Unknown compound cid {}", cid));

    return info_it->second;
}



nlohmann::json chm::App::retrieve_reaction_info_single(const std::string& rid)
{
    {
        std::lock_guard<std::mutex> lock(this->reaction_cache_mutex);
        auto info_it = this->reaction_infos.find(rid);
        if(info_it != this->reaction_infos.end())
            return info_it->second;
    }

    this->retrieve_reaction_infos(std::vector<std::string>{rid});

    std::lock_guard<std::mutex> lock(this->reaction_cache_mutex);
    auto info_it = this->reaction_infos.find(rid);
    if(info_it == this->reaction_infos.end())
        throw std::invalid_argument(fmt::format("Unknown reaction rid '{}'", rid));

    return info_it->second;
}


nlohmann::json chm::App::build_graph(const std::vector<cid_t>& sources, const std::vector<cid_t>& targets, uint16_t max_cost, uint16_t max_paths, bool primary_only)
{
    auto paths{this->find_paths(sources, targets, max_cost, max_paths)};
    nlohmann::json graph_json = this->convert_paths_to_graph(paths, primary_only);

    return {
        {"params", {{"max_cost", max_cost}, {"max_paths", max_paths}, {"primary_only", primary_only}}},
        {"graph", graph_json}
    };
}


std::vector<chm::App::cid_t> chm::App::search_compounds(const std::string& query)
{
    auto results = this->fuzzy->search(query);
    std::vector<cid_t> result_cids(results.size());
    std::transform(results.begin(), results.end(), result_cids.begin(), [](const std::pair<cid_t, double>& x) { return x.first; });

    return result_cids;
}
