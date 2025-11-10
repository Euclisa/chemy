#include "app.hpp"
#include "misc.hpp"


nlohmann::json chm::App::retrieve_compound_info_single(cid_t cid)
{
    if(this->compound_infos.find(cid) == this->compound_infos.end())
        this->retrieve_compound_infos(std::vector<cid_t>{cid});

    return this->compound_infos[cid];
}



nlohmann::json chm::App::retrieve_reaction_info_single(const std::string& rid)
{
    if(this->reaction_infos.find(rid) == this->reaction_infos.end())
        this->retrieve_reaction_infos(std::vector<std::string>{rid});

    return this->reaction_infos[rid];
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


std::pair<std::vector<chm::App::cid_t>, bool> chm::App::search_compounds(const std::string& query, uint32_t offset)
{
    auto results = this->fuzzy.search(query);
    uint32_t results_size = results.size();
    
    if(results_size <= offset)
        throw std::invalid_argument(fmt::format("Invalid offset {} for results list with size {}", offset, results_size));

    uint32_t results_on_page = std::min(results_size - offset, this->PAGE_SIZE);
    std::vector<cid_t> result_cids(results_on_page);
    std::transform(results.begin() + offset, results.begin() + offset + results_on_page, result_cids.begin(), [](const std::pair<cid_t, double>& x) { return x.first; });

    bool is_end = results_size == (offset + results_on_page);

    return std::make_pair(result_cids, is_end);
}