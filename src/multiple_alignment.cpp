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

} // namespace mars
