#pragma once

#include <vector>

namespace mars
{

class motif_store
{
private:
    std::vector<int> bpseq;
    std::vector<std::pair<int, int>> stemloops;

public:
    explicit motif_store(std::vector<int> bpseq): bpseq(std::move(bpseq))
    {}

    void stem_loop_partition(std::vector<int> plevel);

    std::vector<std::pair<int, int>> get_stemloops()
    {
        return stemloops;
    }
};

} // namespace mars
