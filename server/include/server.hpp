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

        void get_compound_infos_request(const httplib::Request& req, httplib::Response& res);
        void get_reaction_infos_request(const httplib::Request& req, httplib::Response& res);
        void build_graph_request(const httplib::Request& req, httplib::Response& res);
        void get_structure_request(const httplib::Request& req, httplib::Response& res);

        void process_request(const std::function<void(const httplib::Request&,httplib::Response&)>& handler, const httplib::Request& req, httplib::Response& res);

    public:

        Server(const fs::path& data_dir, const std::string& listen_addr, uint16_t listen_port);

        void listen();
    };
}

#endif  // CHEMS_SERVER_HPP