#include "app.hpp"
#include "misc.hpp"


nlohmann::json chm::App::retrieve_compound_info_single(cid_t cid)
{
    if(this->compound_infos.find(cid) == this->compound_infos.end())
        this->retrieve_compound_infos(std::vector<cid_t>{cid});

    return this->compound_infos[cid];
}


