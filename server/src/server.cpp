#include "server.hpp"
#include <fstream>
#include <limits>


void chm::Server::compound_info_request(const httplib::Request& req, httplib::Response& res)
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
    auto parse_uint16_field = [&j](const char *field_name)
    {
        if (!j.contains(field_name))
            throw std::invalid_argument(fmt::format("'{}' field must be provided", field_name));
        if (!j[field_name].is_number_integer())
            throw std::runtime_error(fmt::format("'{}' must be an integer", field_name));

        int value = j[field_name].get<int>();
        if (value < 0 || value > std::numeric_limits<uint16_t>::max())
        {
            throw std::out_of_range(
                fmt::format("'{}' must be between 0 and {}", field_name, std::numeric_limits<uint16_t>::max()));
        }

        return static_cast<uint16_t>(value);
    };

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


    uint16_t max_paths = parse_uint16_field("max_paths");
    uint16_t max_cost = parse_uint16_field("max_cost");

    res.set_content(this->build_graph(sources, targets, max_cost, max_paths, primary_only).dump(), "application/json");
}


void chm::Server::adjacent_edges_request(const httplib::Request& req, httplib::Response& res)
{
    if(!req.has_param("cid"))
        throw std::invalid_argument("'cid' paramter must be present");

    cid_t cid = str_to_numeric<cid_t>(req.get_param_value("cid"));
    std::vector<std::vector<cid_t>> edges;
    auto targets_it = this->graph.find(cid);
    if(targets_it != this->graph.end())
    {
        for(const auto& [target_cid, reactions] : targets_it->second)
            edges.push_back({cid, target_cid});
    }
    auto sources_it = this->graph_reverse.find(cid);
    if(sources_it != this->graph_reverse.end())
    {
        for(const auto& [source_cid, reactions] : sources_it->second)
            edges.push_back({source_cid, cid});
    }

    res.set_content(this->convert_paths_to_graph(edges, true).dump(), "application/json");
}


void chm::Server::reaction_info_request(const httplib::Request& req, httplib::Response& res)
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
    fs::path svg_path = this->STRUCTURES_DIR_PATH / (cid_str + ".svg");

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
    std::string page_str = req.has_param("page") ? req.get_param_value("page") : "0";
    std::string sorting_order_str = req.has_param("sorting_order") ? req.get_param_value("sorting_order") : "none";
    std::string sorting_direction = req.has_param("sorting_direction") ? req.get_param_value("sorting_direction") : "asc";
    
    uint32_t page = str_to_numeric<uint32_t>(page_str);
    
    if(sorting_direction != "asc" && sorting_direction != "desc")
        throw std::invalid_argument(fmt::format("Invalid value for 'sorting_direction' parameter: '{}'", sorting_direction));
    
    if(this->sorting.find(sorting_order_str) == this->sorting.end())
        throw std::invalid_argument(fmt::format("Invalid sorting order '{}'", sorting_order_str));
    
    uint32_t offset = page * this->PAGE_SIZE;
    const auto& [sorting_map, sorted_cids] = this->sorting[sorting_order_str];
    
    std::vector<cid_t> search_result;
    const std::vector<cid_t>* source_vector;
    uint32_t total_size;
    
    if (req.has_param("q")) 
    {
        std::string query = req.get_param_value("q");
        search_result = this->search_compounds(query);
        
        if(sorting_order_str != "none")
            std::sort(search_result.begin(), search_result.end(),
                [&sorting_map](cid_t a, cid_t b) { return sorting_map.at(a) < sorting_map.at(b); });
        
        source_vector = &search_result;
        total_size = search_result.size();
    }
    else
    {
        source_vector = &sorted_cids;
        total_size = sorted_cids.size();
    }

    if(total_size == 0)
    {
        res.set_content(nlohmann::json({{"cids", nlohmann::json::array()}, {"is_end", true}}).dump(), "application/json");
        return;
    }

    if(offset >= total_size)
        throw std::invalid_argument(fmt::format("Invalid offset {} for results list with size {}", offset, total_size));
    
    uint32_t results_on_page = std::min(total_size - offset, this->PAGE_SIZE);
    bool is_end = (offset + results_on_page) == total_size;
    
    std::vector<cid_t> query_result;
    query_result.reserve(results_on_page);
    
    if(sorting_direction == "asc")
        std::copy(source_vector->begin() + offset,
                 source_vector->begin() + offset + results_on_page,
                 std::back_inserter(query_result));
    else
        std::copy(source_vector->rbegin() + offset,
                 source_vector->rbegin() + offset + results_on_page,
                 std::back_inserter(query_result));
    
    res.set_content(nlohmann::json({{"cids", query_result}, {"is_end", is_end}}).dump(), "application/json");
}

void chm::Server::process_request(const std::function<void(const httplib::Request&,httplib::Response&)>& handler, const httplib::Request& req, httplib::Response& res)
{
    try {
        handler(req, res);
    } catch (const nlohmann::json::exception &e) {
        nlohmann::json j;
        j["error"] = e.what();
        res.status = 400;
        res.set_content(j.dump(), "application/json");
    } catch (const std::invalid_argument &e) {
        nlohmann::json j;
        j["error"] = e.what();
        res.status = 400;
        res.set_content(j.dump(), "application/json");
    } catch (const std::out_of_range &e) {
        nlohmann::json j;
        j["error"] = e.what();
        res.status = 400;
        res.set_content(j.dump(), "application/json");
    } catch (const std::exception &e) {
        nlohmann::json j;
        j["error"] = e.what();
        res.status = 500;
        res.set_content(j.dump(), "application/json");
    }
}


chm::Server::Server(const fs::path& data_dir) : App(data_dir)
{
    svr.set_mount_point("/", this->UI_DIR_PATH);

    svr.Post("/api/compound_info",
        [this](const httplib::Request& req, httplib::Response& res)
        {
            this->process_request(
                [this](const httplib::Request& req, httplib::Response& res) { return this->compound_info_request(req, res); },
                req, res
            );
        });
        
    svr.Post("/api/reaction_info",
        [this](const httplib::Request& req, httplib::Response& res)
        {
            this->process_request(
                [this](const httplib::Request& req, httplib::Response& res) { return this->reaction_info_request(req, res); },
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
    
    std::string svg_route = fmt::format("/{}/(-?\\d+)\\.svg", this->STRUCTURES_DIR_UI_PATH.string());
    svr.Get(svg_route,
        [&](const httplib::Request& req, httplib::Response& res)
        {
            this->process_request(
                [this](const httplib::Request& req, httplib::Response& res) { return this->get_structure_request(req, res); },
                req, res
            );
        });
    
    svr.Get("/api/search",
        [&](const httplib::Request &req, httplib::Response &res)
        {
            this->process_request(
                [this](const httplib::Request& req, httplib::Response& res) { return this->search_request(req, res); },
                req, res
            );
        });
    
    svr.Get("/api/adjacent_edges",
        [&](const httplib::Request &req, httplib::Response &res)
        {
            this->process_request(
                [this](const httplib::Request& req, httplib::Response& res) { return this->adjacent_edges_request(req, res); },
                req, res
            );
        });
    
    svr.set_post_routing_handler([](const httplib::Request& req, httplib::Response& res) {
        using namespace fmt::literals;

        fmt::color status_color;
        if (res.status >= 500) status_color = fmt::color::yellow;
        else if (res.status >= 400) status_color = fmt::color::red;
        else if (res.status >= 300) status_color = fmt::color::gray;
        else if (res.status >= 200) status_color = fmt::color::green;
        else status_color = fmt::color::white;

        std::string log_msg;

        log_msg += fmt::format("Request: ");
        log_msg += fmt::format(fg(fmt::color::cyan) | fmt::emphasis::bold, "{} ", req.path);
        log_msg += fmt::format(fg(fmt::color::white), "({}) ", req.method);
        log_msg += fmt::format(fg(fmt::color::yellow), "from {} ", req.remote_addr);
        log_msg += fmt::format(fg(status_color) | fmt::emphasis::bold, "[{}]", res.status);

        if (auto it = req.headers.find("Content-Length"); it != req.headers.end()) {
            log_msg += fmt::format(fg(fmt::color::gray), " len:{}", it->second);
        }

        if (res.status >= 400 && res.body.size() > 0 && res.body.size() < 1000)
            log_msg += fmt::format(fg(fmt::color::red) | fmt::emphasis::italic, " Body: {}", res.body);

        log_msg += "\n";

        fmt::print("{}", log_msg);
    });
}



void chm::Server::listen(const std::string& listen_addr, uint16_t listen_port)
{
    fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "🚀 Server is up at ");
    fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold, "{}:{}", listen_addr, listen_port);
    fmt::print("\n");
    this->svr.listen(listen_addr, listen_port);
}
