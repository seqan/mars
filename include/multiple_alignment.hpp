#pragma once

#include <string>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>

namespace mars
{

/*!
 * \brief A multiple alignment representation.
 * \tparam alphabet_type The alphabet type of the contained sequences.
 */
template<seqan3::alphabet alphabet_type>
struct multiple_alignment
{
    std::vector<std::vector<seqan3::gapped<alphabet_type>>> sequences;
    std::vector<std::string> names;
};

/*!
 * \brief The multiple alignment type used in MaRs.
 * \relates multiple_alignment
 */
typedef multiple_alignment<seqan3::rna15> msa_type;

} // namespace mars
