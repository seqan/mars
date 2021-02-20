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
 * \tparam Alphabet The alphabet type of the contained sequences.
 */
template<seqan3::alphabet Alphabet>
struct MultipleAlignment
{
    std::vector<std::vector<seqan3::gapped<Alphabet>>> sequences;
    std::vector<std::string> names;
};

/*!
 * \brief The multiple alignment type used in MaRs.
 * \relates multiple_alignment
 */
typedef MultipleAlignment<seqan3::rna15> Msa;

//! \brief Type for sequence id. We expect to have less than 65k sequences.
using SeqNum = uint16_t;

//! \brief Type for positions within a sequence.
using SeqLen = size_t;

} // namespace mars
