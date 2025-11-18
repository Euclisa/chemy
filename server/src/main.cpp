#include "server.hpp"


int main()
{
    std::string data_dir = PROJECT_DIR "/data";
    chm::Server server(data_dir);
    server.listen("0.0.0.0", 8080);

    return 0;
}