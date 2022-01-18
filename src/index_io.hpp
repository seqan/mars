#pragma once

#include <seqan3/std/filesystem>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

namespace mars
{

//! \brief The type of a bi-directional index over the 4-letter DNA alphabet.
using Index = seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::collection>;

/*!
 * \brief Read a FASTA file of sequences.
 * \param[out] seqs The object where the sequences can be stored.
 * \param[out] names The object where the sequence names can be stored.
 * \param[in] filepath The file path where the sequences and names can be read from.
 */
void read_genome(std::vector<seqan3::dna4_vector> & seqs,
                 std::vector<std::string> & names,
                 std::filesystem::path const & filepath);

/*!
 * \brief Archive an index and store it in a file on disk.
 * \param[in] index The index that should be archived.
 * \param[in] names The sequence names.
 * \param[in,out] indexpath The path of the index output file.
 * \param[in] compress Use gzip compression for creating the index file.
 */
void write_index(Index const & index, std::vector<std::string> const & names, std::filesystem::path & indexpath,
                 bool compress);

/*!
 * \brief Unarchive an index and read it from a file on disk.
 * \param[out] index The index object which is filled with the contents of the file.
 * \param[out] names The sequence names.
 * \param[in] indexpath The path of the index input file.
 * \return whether an index could be parsed.
 */
bool read_index(Index & index, std::vector<std::string> & names, std::filesystem::path & indexpath);

} // namespace mars
