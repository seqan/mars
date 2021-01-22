#include <seqan3/io/sequence_file/input.hpp>

#include "format_clustal.hpp"
#include "input_output.hpp"

namespace mars
{

msa_type read_msa(std::istream & stream)
{
    return read_clustal_file<seqan3::rna15>(stream);
}

msa_type read_msa(std::filesystem::path const & filepath)
{
    return read_clustal_file<seqan3::rna15>(filepath);
}

std::vector<seqan3::dna4_vector> read_genome(std::filesystem::path const & filepath)
{
    std::vector<seqan3::dna4_vector> seqs{};

    struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
    {
        using sequence_alphabet = seqan3::dna4;
        using sequence_legal_alphabet = seqan3::dna15;
    };

    for (auto && record : seqan3::sequence_file_input<dna4_traits>{filepath})
        seqs.push_back(std::move(seqan3::get<seqan3::field::seq>(record)));

    return std::move(seqs);
}

} // namespace mars
