#pragma once

#include <tuple>
#include <vector>

#include <seqan3/std/filesystem>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>

#include "bi_alphabet.hpp"
#include "index_io.hpp"

namespace mars
{

//! \brief Provides a bi-directional search step-by-step with backtracking.
class BiDirectionalIndex
{
private:
    //! \brief The index in which the search is performed.
    Index index;

    //! \brief The history of queries (needed for backtracking).
    std::vector<seqan3::dna4_vector> queries;

    //! \brief The history of scores;
    std::vector<float> scores;

    //! \brief The xdrop parameter.
    unsigned char const xdrop_dist;

public:
    //! \brief The resulting matches of the current search step.
    std::vector<std::pair<uint16_t, size_t>> matches{};

    /*!
     * \brief Constructor for a bi-directional search.
     * \param xdrop The xdrop parameter.
     */
    explicit BiDirectionalIndex(unsigned char xdrop):
        index{},
        queries{},
        scores{},
        xdrop_dist{xdrop},
        matches{}
    {
        queries.emplace_back();
        scores.emplace_back(0);
    }

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
    void create(std::filesystem::path const & filepath);

    /*!
     * \brief Append a character to the 5' (left) side of the query.
     * \param item The character to be added.
     * \param left Whether the loop is at the 5' side.
     */
    void append_loop(std::pair<float, seqan3::rna4> item, bool left);

    /*!
     * \brief Append a character pair at both sides of the query.
     * \param stem_item The score and characters to be added.
     */
    void append_stem(std::pair<float, bi_alphabet<seqan3::rna4>> stem_item);

    //! \brief Revert the previous append step, which shrinks the query by one or two characters.
    void backtrack();

    /*!
     * \brief Whether the search should be aborted through the xdrop condition.
     * \return True if the score dropped over the previous x elements, false otherwise.
     */
    [[nodiscard]] bool xdrop() const;

    /*!
     * \brief Get the total score of the current query.
     * \return the score.
     */
    [[nodiscard]] float get_score() const;

    /*!
     * \brief Perform the search with the current query and store the result in `matches`.
     * \return The number of matches that have been found.
     */
    size_t compute_matches();
};

} // namespace mars
