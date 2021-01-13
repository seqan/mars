#include "format_clustal.hpp"
#include "input_output.hpp"

mars::msa_type mars::read_msa(std::istream & stream)
{
    return read_clustal_file<seqan3::rna15>(stream);
}

mars::msa_type mars::read_msa(std::filesystem::path const & filepath)
{
    return read_clustal_file<seqan3::rna15>(filepath);
}
