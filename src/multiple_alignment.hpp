#pragma once

#include <seqan3/std/filesystem>
#include <istream>
#include <string>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>

namespace mars
{

/*!
 * \brief A multiple alignment representation.
 * \tparam AlphabetType The alphabet type of the contained sequences.
 */
template<seqan3::alphabet AlphabetType>
struct MultipleAlignment
{
    using Alphabet = AlphabetType;

    std::vector<std::vector<seqan3::gapped<Alphabet>>> sequences;
    std::vector<std::string> names;
};

/*!
 * \brief The multiple alignment type used in MaRs.
 * \relates multiple_alignment
 */
typedef MultipleAlignment<seqan3::rna15> Msa;

//! \brief Type for sequence id. We expect to have less than 4G sequences.
using SeqNum = size_t;

//! \brief Type for positions within a sequence.
using SeqLen = size_t;

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \param stream The input stream where the alignment is parsed from.
 * \return The alignment.
 */
Msa read_msa(std::istream & stream);

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \param filepath The file where the alignment is stored.
 * \return The alignment.
 */
Msa read_msa(std::filesystem::path const & filepath);

} // namespace mars
