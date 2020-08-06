#include <seqan3/std/ranges>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/zip.hpp>

#include "motif_store.hpp"

void mars::motif_store::stem_loop_partition(std::vector<int> plevel)
{
    struct pk_info
    {
        int level;
        bool closing;
        std::pair<int, int> previous;
    };
    std::vector<pk_info> pk_infos{};

    for (auto && [idx, bp, pk] : seqan3::views::zip(std::ranges::views::iota(0), bpseq, plevel))
    {
        if (pk == -1) // skip unpaired
            continue;

        while (pk + 1 > pk_infos.size()) // allocate a new pseudoknot layer
            pk_infos.push_back({0, false, {0, 0}});

        pk_info & status = pk_infos[pk];
        if (bp < idx) // close an interaction
        {
            status.previous = {bp, idx};
            status.closing = true;
            if (--status.level == 0)
                stemloops.push_back(status.previous);
        }
        else if (status.closing) // open an interaction (after closing the previous)
        {
            if (status.level > 0)
                stemloops.push_back(status.previous);
            status.level = 1;
            status.closing = false;
        }
        else // open another interaction
        {
            ++status.level;
        }
    }

    seqan3::debug_stream << stemloops << "\n";
}
