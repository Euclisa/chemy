#include "server.hpp"


int main()
{
    chm::Server server("data", "0.0.0.0", 8080);
    server.listen();

    return 0;
}