#pragma once

#include <seqan3/std/filesystem>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

namespace mars
{

//! \brief The type of a bi-directional index over the 4-letter DNA alphabet.
using Index = seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::collection>;

/*!
 * \brief Read a FASTA file of sequences.
 * \param filepath The file path where the sequences can be read from.
 * \return A vector of DNA sequences.
 */
std::vector<seqan3::dna4_vector> read_genome(std::filesystem::path const & filepath);

/*!
 * \brief Archive an index and store it in a file on disk.
 * \param index The index that should be archived.
 * \param indexpath The path of the index output file.
 */
void write_index(Index const & index, std::filesystem::path & indexpath);

/*!
 * \brief Unarchive an index and read it from a file on disk.
 * \param index The index object which is filled with the contents of the file.
 * \param indexpath The path of the index input file.
 * \return whether an index could be parsed.
 */
bool read_index(Index & index, std::filesystem::path const & indexpath);

} // namespace mars
