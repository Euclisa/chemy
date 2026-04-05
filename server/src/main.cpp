#include "server.hpp"

#include <iostream>


int main(int argc, char **argv)
{
    std::string data_dir = argc >= 2 ? argv[1] : PROJECT_DIR "/data";
    std::string listen_addr = argc >= 3 ? argv[2] : "0.0.0.0";
    uint16_t listen_port = argc >= 4 ? chm::str_to_numeric<uint16_t>(argv[3]) : 8080;

    chm::Server server(data_dir);
    server.listen(listen_addr, listen_port);

    return 0;
}
