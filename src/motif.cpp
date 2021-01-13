#include <seqan3/std/ranges>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/zip.hpp>

#include "bi_alphabet.hpp"
#include "motif.hpp"

namespace mars
{

stemloop_type detect_stem_loops(std::vector<int> const & bpseq, std::vector<int> const & plevel)
{
    struct pk_info
    {
        int level;
        bool closing;
        std::pair<int, int> previous;
    };
    std::vector<pk_info> pk_infos{};
    stemloop_type stemloops;

    // 0-based indices
    for (auto &&[idx, bp, pk] : seqan3::views::zip(std::ranges::views::iota(0), bpseq, plevel))
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
    return std::move(stemloops);
}

stem_loop_motif analyze_stem_loop(msa_type const & msa, std::vector<int> const & bpseq, std::pair<int, int> const & pos)
{
    stem_loop_motif motif{};
    motif.bounds = pos;

    return std::move(motif);
}

} // namespace mars
