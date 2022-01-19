#include <seqan3/std/algorithm>
#include <seqan3/std/iterator>

#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/utility/views/zip.hpp>

#include "format_clustal.hpp"
#include "format_stockholm.hpp"
#include "multiple_alignment.hpp"
#include "structure.hpp"

namespace mars
{

Msa read_msa(std::filesystem::path const & filepath)
{
    if (filepath.extension() == std::filesystem::path{".aln"} ||
        filepath.extension() == std::filesystem::path{".msa"})
    {
        Msa msa = std::move(read_clustal_file<typename Msa::Alphabet>(filepath));
        compute_structure(msa);
        return std::move(msa);
    }
    else if (filepath.extension() == std::filesystem::path{".sth"} ||
             filepath.extension() == std::filesystem::path{".stk"} ||
             filepath.extension() == std::filesystem::path{".sto"})
    {
        return std::move(read_stockholm_file<typename Msa::Alphabet>(filepath));
    }
    else
    {
        throw seqan3::file_open_error{"Unknown file extension for the alignment file " + filepath.string() + "."};
    }
}

void compute_structure(Msa & msa)
{
    // Convert names
    std::list<std::string> names{msa.names.size()};
    std::ranges::copy(msa.names, names.begin());

    // Convert sequences
    std::list<std::string> seqs{msa.sequences.size()};
    for (auto && [src, trg] : seqan3::views::zip(msa.sequences | seqan3::views::to_char, seqs))
        std::ranges::copy(src, std::cpp20::back_inserter(trg));

    msa.structure = std::move(run_ipknot(names, seqs));
}

} // namespace mars
