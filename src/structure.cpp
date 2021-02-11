#include <seqan3/std/algorithm>
#include <seqan3/std/iterator>

#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/zip.hpp>

#include "structure.hpp"

namespace mars
{

std::pair<std::vector<int>, std::vector<int>> compute_structure(Msa const & msa)
{
    // Convert names
    std::list<std::string> names{msa.names.size()};
    std::ranges::copy(msa.names, names.begin());

    // Convert sequences
    std::list<std::string> seqs{msa.sequences.size()};
    for (auto && [src, trg] : seqan3::views::zip(msa.sequences | seqan3::views::to_char, seqs))
        std::ranges::copy(src, std::cpp20::back_inserter(trg));

    return std::move(run_ipknot(names, seqs));
}

}
