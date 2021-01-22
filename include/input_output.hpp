#pragma once

#include <seqan3/std/filesystem>
#include <istream>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>

#include "multiple_alignment.hpp"

namespace mars
{

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \param stream The input stream where the alignment is parsed from.
 * \return The alignment.
 */
msa_type read_msa(std::istream & stream);

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \param filepath The file where the alignment is stored.
 * \return The alignment.
 */
msa_type read_msa(std::filesystem::path const & filepath);

/*!
 * \brief Read a FASTA file of sequences.
 * \param filepath The file path where the sequences can be read from.
 * \return A vector of DNA sequences.
 */
std::vector<seqan3::dna4_vector> read_genome(std::filesystem::path const & filepath);

} // namespace mars
