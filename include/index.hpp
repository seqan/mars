#pragma once

#include <tuple>
#include <vector>

#include <seqan3/std/filesystem>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

namespace mars
{

//! \brief The type of a bi-directional index over the 4-letter DNA alphabet.
using index_type = seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::collection>;

/*!
 * \brief Create an index of a genome from the specified file.
 * \param filepath The filepath to the file.
 * \return a valid index of the genome.
 * \throws seqan3::file_open_error if neither `filepath` nor `filepath.marsindex` exist.
 *
 * \details
 *
 * This function has two modes:
 *
 * 1. If `filepath.marsindex` exists: Read the already created index from this file.
 * 2. Else if `filepath` exists: Read sequences from this file, create an index
 *    and write the index to `filepath.marsindex`.
 */
index_type create_index(std::filesystem::path const & filepath);

//! \brief Provides a bi-directional search step-by-step with backtracking.
class bi_directional_search
{
private:
    //! \brief The index in which the search is performed.
    index_type const & index;
    //! \brief The history of queries (needed for backtracking).
    std::vector<seqan3::dna4_vector> queries{};

public:
    //! \brief The resulting matches of the current search step.
    std::vector<std::pair<size_t, size_t>> matches{};

    /*!
     * \brief Constructor for a bi-directional search.
     * \param bi_dir_index The index in which the search is performed.
     */
    explicit bi_directional_search(index_type const & bi_dir_index): index(bi_dir_index), queries{}, matches{}
    {
        queries.emplace_back();
    }

    /*!
     * \brief Append a character to the 5' (left) side of the query.
     * \param left The character to be added.
     */
    void append_5prime(seqan3::dna4 left);

    /*!
     * \brief Append a character to the 3' (right) side of the query.
     * \param right The character to be added.
     */
    void append_3prime(seqan3::dna4 right);

    /*!
     * \brief Append a character pair at both sides of the query.
     * \param left The character to be added on the 5' (left) side.
     * \param right The character to be added on the 3' (right) side.
     */
    void append_stem(seqan3::dna4 left, seqan3::dna4 right);

    /*!
     * \brief Revert the previous append step, which shrinks the query by one or two characters.
     * \return true if the operation is successful, false if the history of queries is empty.
     */
    bool backtrack();

    /*!
     * \brief Perform the search with the current query and store the result in `matches`.
     * \return The number of matches that have been found.
     */
    size_t compute_matches();
};

} // namespace mars
