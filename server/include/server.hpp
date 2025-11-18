#ifndef CHEMS_SERVER_HPP
#define CHEMS_SERVER_HPP

#include "httplib.h"
#include <string>

#include "app.hpp"


namespace chm
{
    class Server : App
    {
    private:
        httplib::Server svr;

        std::string listen_addr;
        uint16_t listen_port;

        void compound_info_request(const httplib::Request& req, httplib::Response& res);
        void reaction_info_request(const httplib::Request& req, httplib::Response& res);
        void build_graph_request(const httplib::Request& req, httplib::Response& res);
        void adjacent_edges_request(const httplib::Request& req, httplib::Response& res);
        void get_structure_request(const httplib::Request& req, httplib::Response& res);
        void search_request(const httplib::Request& req, httplib::Response& res);

        void process_request(const std::function<void(const httplib::Request&,httplib::Response&)>& handler, const httplib::Request& req, httplib::Response& res);

    public:

        Server(const fs::path& data_dir);

        void listen(const std::string& listen_addr, uint16_t listen_port);
    };
}

#endif  // CHEMS_SERVER_HPP