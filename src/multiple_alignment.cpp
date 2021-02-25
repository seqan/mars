#include "format_clustal.hpp"
#include "multiple_alignment.hpp"

namespace mars
{

Msa read_msa(std::istream & stream)
{
    return std::move(read_clustal_file<typename Msa::Alphabet>(stream));
}

Msa read_msa(std::filesystem::path const & filepath)
{
    return std::move(read_clustal_file<typename Msa::Alphabet>(filepath));
}

} // namespace mars
