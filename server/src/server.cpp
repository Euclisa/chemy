#include "server.hpp"
#include <fstream>


void chm::Server::get_compound_infos_request(const httplib::Request& req, httplib::Response& res)
{
    std::string body = req.body;

    nlohmann::json cids_json = nlohmann::json::parse(body);;
    std::vector<cid_t> cids;

    if(!cids_json.is_array())
        throw std::runtime_error("Expected JSON array of cids");

    for(const auto& item : cids_json)
    {

        if (!item.is_number_integer())
            throw std::runtime_error("Expected numeric cid values");

        cids.push_back(item.get<cid_t>());
    }

    res.set_content(this->retrieve_compound_infos_json(cids).dump(), "application/json");
}


void chm::Server::build_graph_request(const httplib::Request& req, httplib::Response& res)
{
    std::string body = req.body;

    nlohmann::json j = nlohmann::json::parse(body);

    std::vector<cid_t> targets;
    std::vector<cid_t> sources;

    if (j.contains("targets") && j["targets"].is_array())
    {
        for (const auto& item : j["targets"])
        {
            if (!item.is_number_integer())
                throw std::invalid_argument("Entries of 'targets' must be numeric");
            targets.push_back(item.get<cid_t>());
        }
    }
    else
        throw std::invalid_argument("Invalid 'targets' field");


    if (j.contains("sources") && j["sources"].is_array())
    {
        for (const auto& item : j["sources"])
        {
            if (!item.is_number_integer())
                throw std::invalid_argument("Entries of 'sources' must be numeric");
            sources.push_back(item.get<cid_t>());
        }
    }
    else
        throw std::invalid_argument("Invalid 'sources' field");


    bool primary_only;
    if (j.contains("primary_only"))
    {
        if (!j["primary_only"].is_boolean())
            throw std::runtime_error("'primary_only' must be a boolean");
        primary_only = j["primary_only"].get<bool>();
    }
    else
        throw std::invalid_argument("'primary_only' field must be provided");


    uint16_t max_paths ;
    if (j.contains("max_paths")) 
    {
        if (!j["max_paths"].is_number_integer())
            throw std::runtime_error("'max_paths' must be an integer");
        max_paths = static_cast<uint16_t>(j["max_paths"].get<int>());
    }
    else
        throw std::invalid_argument("'max_paths' field must be provided");


    uint16_t max_len;
    if (j.contains("max_len")) 
    {
        if (!j["max_len"].is_number_integer())
            throw std::runtime_error("'max_len' must be an integer");
        max_len = static_cast<uint16_t>(j["max_len"].get<int>());
    }
    else
        throw std::invalid_argument("'max_len' field must be provided");

    res.set_content(this->build_graph(targets, sources, primary_only, max_paths, max_len).dump(), "application/json");
}


void chm::Server::get_reaction_infos_request(const httplib::Request& req, httplib::Response& res)
{
    std::string body = req.body;
    nlohmann::json rids_json = nlohmann::json::parse(body);;
    std::vector<std::string> rids;

    if(!rids_json.is_array())
        throw std::runtime_error("Expected JSON array of cids");

    for(const auto& item : rids_json)
        rids.push_back(item.get<std::string>());

    res.set_content(this->retrieve_reaction_infos_json(rids).dump(), "application/json");
}


void chm::Server::get_structure_request(const httplib::Request& req, httplib::Response& res)
{
    std::string cid_str = req.matches[1];
    cid_t cid = str_to_numeric<cid_t>(cid_str);
    fs::path svg_path = this->UI_DIR_PATH / this->STRUCTURES_DIR_UI_PATH / (cid_str + ".svg");

    if (!fs::exists(svg_path)) {
        try {
            this->generate_compound_structure_svg(cid); 
        } catch (const std::exception& e) {
            res.status = 500;
            res.set_content(std::string("Failed to generate SVG: ") + e.what(), "text/plain");
            return;
        }
    }

    std::ifstream file(svg_path);
    if (!file.is_open()) {
        res.status = 500;
        res.set_content("Could not open SVG file after generation.", "text/plain");
        return;
    }

    std::string svg_content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    res.set_content(svg_content, "image/svg+xml");
}


void chm::Server::search_request(const httplib::Request& req, httplib::Response& res)
{
    std::string page_str = req.has_param("page") ? req.get_param_value("page") : "1";
    std::string sorting_order_str = req.has_param("sort") ? req.get_param_value("sort") : "complexity";
    uint32_t page = str_to_numeric<uint32_t>(page_str);

    if(page == 0)
        throw std::invalid_argument("Page index must be >= 1");

    uint32_t offset = (page - 1) * this->PAGE_SIZE;

    if(this->sorting.find(sorting_order_str) == this->sorting.end())
        throw std::invalid_argument(fmt::format("Invalid sorting order '{}'", sorting_order_str));
    
    const auto& [sorting_map, sorted_cids] = this->sorting[sorting_order_str];

    std::vector<cid_t> query_result;
    if (req.has_param("q")) 
    {
        std::string query = req.get_param_value("q");
        query_result = this->search_compounds(query, offset);
        
        std::sort(query_result.begin(), query_result.end(),
            [&sorting_map](cid_t a, cid_t b) { return sorting_map.at(a) < sorting_map.at(b); });
    }
    else
    {
        uint32_t cids_size = sorted_cids.size();

        if(offset >= sorted_cids.size())
            throw std::invalid_argument(fmt::format("Invalid offset {} for results list with size {}", offset, cids_size));

        uint32_t results_on_page = std::min(cids_size - offset, this->PAGE_SIZE);
        query_result.resize(results_on_page);
        std::copy(sorted_cids.begin() + offset,
            sorted_cids.begin() + offset + results_on_page,
            query_result.begin());
    }

    res.set_content(nlohmann::json(query_result).dump(), "application/json");
}


void chm::Server::process_request(const std::function<void(const httplib::Request&,httplib::Response&)>& handler, const httplib::Request& req, httplib::Response& res)
{
    try {
        handler(req, res);
    } catch (const std::exception &e) {
        nlohmann::json j;
        j["error"] = e.what();
        res.status = 400;
        res.set_content(j.dump(), "application/json");
    }
}


chm::Server::Server(const fs::path& data_dir, const std::string& listen_addr, uint16_t listen_port) : App(data_dir), listen_addr(listen_addr), listen_port(listen_port)
{
    svr.set_mount_point("/", this->UI_DIR_PATH);

    svr.Post("/api/compound_info",
        [this](const httplib::Request& req, httplib::Response& res)
        {
            this->process_request(
                [this](const httplib::Request& req, httplib::Response& res) { return this->get_compound_infos_request(req, res); },
                req, res
            );
        });
        
    svr.Post("/api/reaction_info",
        [this](const httplib::Request& req, httplib::Response& res)
        {
            this->process_request(
                [this](const httplib::Request& req, httplib::Response& res) { return this->get_reaction_infos_request(req, res); },
                req, res
            );
        });
    
    svr.Post("/api/build_graph",
        [this](const httplib::Request& req, httplib::Response& res)
        {
            this->process_request(
                [this](const httplib::Request& req, httplib::Response& res) { return this->build_graph_request(req, res); },
                req, res
            );
        });
    
    std::string svg_route = this->STRUCTURES_DIR_UI_PATH.string() + R"(\d+)\.svg)";
    svr.Get(svg_route,
        [&](const httplib::Request& req, httplib::Response& res)
        {
            this->process_request(
                [this](const httplib::Request& req, httplib::Response& res) { return this->get_structure_request(req, res); },
                req, res
            );
        });
    
    svr.Get("/search",
        [&](const httplib::Request &req, httplib::Response &res)
        {
            this->process_request(
                [this](const httplib::Request& req, httplib::Response& res) { return this->search_request(req, res); },
                req, res
            );
        });
}



void chm::Server::listen()
{
    this->svr.listen(this->listen_addr, this->listen_port);
}