#pragma once

#include <seqan3/std/filesystem>
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
    //! \brief The underlying alphabet type.
    using Alphabet = AlphabetType;

    //! \brief The gapped sequences.
    std::vector<std::vector<seqan3::gapped<Alphabet>>> sequences;
    //! \brief The sequence names or identifiers.
    std::vector<std::string> names;
    //! \brief The consensus structure of the alignment.
    std::pair<std::vector<int>, std::vector<int>> structure;
};

/*!
 * \brief The multiple alignment type used in MaRs.
 * \relates multiple_alignment
 */
typedef MultipleAlignment<seqan3::rna15> Msa;

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \param filepath The file where the alignment is stored.
 * \return The alignment.
 */
Msa read_msa(std::filesystem::path const & filepath);

/*!
 * \brief Compute the secondary structure of a given multiple structural alignment.
 * \param msa The multiple structural alignment.
 */
void compute_structure(Msa & msa);

} // namespace mars
