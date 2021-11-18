#pragma once

#include <seqan3/std/filesystem>
#include <istream>
#include <stack>
#include <string>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/exception.hpp>

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

    /*!
     * \brief Read a wuss string and extract the base pair positions and pseudoknot levels for the MSA.
     * \param[in] wuss_string The string of wuss characters to parse the structure from.
     */
    void parse_structure(std::vector<seqan3::wuss51> const & wuss_string)
    {
        structure.first.resize(wuss_string.size(), -1);
        structure.second.resize(wuss_string.size(), -1);

        std::stack<int> brackets[seqan3::max_pseudoknot_depth<seqan3::wuss51>];
        int pos = 0ul;
        for (seqan3::wuss51 symbol: wuss_string)
        {
            int const id = seqan3::pseudoknot_id(symbol).value_or(-1);

            if (symbol.is_pair_open())
            {
                brackets[id].push(pos);
            }
            else if (symbol.is_pair_close())
            {
                if (!brackets[id].empty())
                {
                    structure.first[pos] = brackets[id].top();
                    structure.first[brackets[id].top()] = pos;
                    structure.second[pos] = structure.second[brackets[id].top()] = id;
                    brackets[id].pop();
                }
                else
                {
                    throw seqan3::parse_error{std::string{
                        "Invalid bracket notation: Unpaired closing bracket at position "} + std::to_string(pos) + "."};
                };
            }
            // no actions for unpaired
            ++pos;
        }
        for (std::stack<int> const & stack: brackets)
        {
            if (!stack.empty())
            {
                throw seqan3::parse_error{std::string{"Invalid bracket notation: Unpaired opening bracket at position "}
                                          + std::to_string(stack.top()) + "."};
            }
        }
    }
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
 * \param filepath The file where the alignment is stored.
 * \return The alignment.
 */
Msa read_msa(std::filesystem::path const & filepath);

} // namespace mars
