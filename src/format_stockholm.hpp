#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/filesystem>
#include <fstream>
#include <stack>
#include <string>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/range/detail/misc.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/views/detail/istreambuf_view.hpp>
#include <seqan3/io/views/detail/take_exactly_view.hpp>
#include <seqan3/io/views/detail/take_line_view.hpp>
#include <seqan3/utility/views/type_reduce.hpp>
#include <seqan3/utility/views/zip.hpp>

#include "multiple_alignment.hpp"

namespace mars
{

/*!
 * \brief Read a wuss string and extract the base pair positions and pseudoknot levels for the MSA.
 * \param[out] structure The structure storage of the MSA.
 * \param[in] wuss_string The string of wuss characters to parse the structure from.
 */
void parse_structure(std::pair<std::vector<int>, std::vector<int>> & structure,
                     std::vector<seqan3::wuss51> const & wuss_string)
{
    structure.first.resize(wuss_string.size(), -1);
    structure.second.resize(wuss_string.size(), -1);

    std::stack<int> brackets[seqan3::max_pseudoknot_depth<seqan3::wuss51>];
    int pos = 0;
    for (seqan3::wuss51 symbol: wuss_string)
    {
        int const pkid = seqan3::pseudoknot_id(symbol).value_or(-1);

        if (symbol.is_pair_open())
        {
            brackets[pkid].push(pos);
        }
        else if (symbol.is_pair_close())
        {
            if (!brackets[pkid].empty())
            {
                structure.first[pos] = brackets[pkid].top();
                structure.first[brackets[pkid].top()] = pos;
                int reduced_pk = pkid;
                while (reduced_pk > 0 && brackets[reduced_pk - 1].empty())
                    --reduced_pk;
                structure.second[pos] = structure.second[brackets[pkid].top()] = reduced_pk;
                brackets[pkid].pop();
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

/*!
 * \brief Read a Stockholm file (*.sth) into a multiple alignment representation.
 * \tparam alphabet_type The alphabet type of the sequences.
 * \param stream The input stream where the alignment is parsed from.
 * \return The alignment.
 */
template<seqan3::alphabet alphabet_type>
MultipleAlignment<alphabet_type> read_stockholm_file(std::istream & stream)
{
    MultipleAlignment<alphabet_type> msa;

    // Define a lambda function to decide if a character is legal.
    auto check_legal_alphabet = [] (char const chr)
    {
        using legal_alphabet_type = std::conditional_t<seqan3::nucleotide_alphabet<alphabet_type>,
                                                       seqan3::rna15,
                                                       seqan3::aa27>;
        auto constexpr is_legal_alph = seqan3::char_is_valid_for<seqan3::gapped<legal_alphabet_type>>;
        if (!is_legal_alph(chr))
            throw seqan3::parse_error{"Encountered an unexpected letter: char_is_valid_for<" +
                                      seqan3::detail::type_name_as_string<legal_alphabet_type> +
                                      "> evaluated to false on " + seqan3::detail::make_printable(chr)};
        return chr;
    };

    auto stream_view = seqan3::detail::istreambuf(stream);

    // skip initial whitespace and check if file starts with "# STOCKHOLM 1.0"
    seqan3::detail::consume(stream_view | seqan3::detail::take_until(!seqan3::is_space));
    if (!std::ranges::equal(stream_view | seqan3::detail::take_exactly_or_throw(15), std::string{"# STOCKHOLM 1.0"}))
        throw seqan3::parse_error{"Expected to read '# STOCKHOLM 1.0' in the beginning of the file."};

    // skip rest of the line and go to next line
    seqan3::detail::consume(stream_view | seqan3::detail::take_line);
    seqan3::detail::consume(stream_view | seqan3::detail::take_until(!seqan3::is_space));

    std::size_t idx = 0;
    bool first_block = true;
    std::vector<seqan3::wuss51> wuss_string{};

    do // read line-wise
    {
        if (seqan3::is_char<'#'>(*stream_view.begin())) // skip or parse #= line
        {
            // found secondary structure
            if (std::ranges::equal(stream_view | seqan3::detail::take_exactly(12), std::string{"#=GC SS_cons"}))
            {
                seqan3::detail::consume(stream_view | seqan3::detail::take_until_or_throw(!seqan3::is_blank));
                std::ranges::copy(stream_view | seqan3::detail::take_until_or_throw(seqan3::is_space)
                                  | seqan3::views::char_to<seqan3::wuss51>,
                                  std::cpp20::back_inserter(wuss_string));
                idx = 0;
                first_block = false;
            }
            seqan3::detail::consume(stream_view | seqan3::detail::take_line);
        }
        else if (seqan3::is_space(*stream_view.begin())) // skip empty lines
        {
            seqan3::detail::consume(stream_view | seqan3::detail::take_until(!seqan3::is_space));
        }
        else // parse sequence line
        {
            // parse the sequence name
            std::string name{};
            std::ranges::copy(stream_view | seqan3::detail::take_until_or_throw(seqan3::is_blank),
                              std::cpp20::back_inserter(name));

            if (first_block)
            {
                // add new name and sequence
                msa.names.push_back(name);
                msa.sequences.push_back(std::vector<seqan3::gapped<alphabet_type>>{});
            }
            else if (idx >= msa.names.size()) // check for inconsistencies
            {
                throw seqan3::parse_error{"Inconsistent alignment depth in the input file."};
            }
            else if (!std::ranges::equal(name, msa.names[idx])) // validate the sequence name
            {
                throw seqan3::parse_error{"Expected to read '" + msa.names[idx] + "' in the input file."};
            }

            // go to the beginning of sequence
            seqan3::detail::consume(stream_view | seqan3::detail::take_until_or_throw(!seqan3::is_blank));

            // read the sequence
            std::ranges::copy(stream_view | seqan3::detail::take_until_or_throw(seqan3::is_space)
                              | std::views::filter(!seqan3::is_digit)
                              | std::views::transform(check_legal_alphabet)
                              | seqan3::views::char_to<seqan3::gapped<alphabet_type>>,
                              std::cpp20::back_inserter(msa.sequences[idx]));

            // go to the next line
            seqan3::detail::consume(stream_view | seqan3::detail::take_line);
            ++idx;
        }

    }
    while (!seqan3::is_char<'/'>(*stream_view.begin())); // end of stockholm record

    parse_structure(msa.structure, wuss_string);
    return std::move(msa);
}

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \tparam alphabet_type The alphabet type of the sequences.
 * \param filepath The file where the alignment is stored.
 * \return The alignment.
 */
template<seqan3::alphabet alphabet_type>
MultipleAlignment<alphabet_type> read_stockholm_file(std::filesystem::path const & filepath)
{
    // Open filepath as stream.
    std::ifstream stream(filepath, std::ios_base::in | std::ios::binary);
    if (!stream.good())
        throw seqan3::file_open_error{"Could not open file " + filepath.string() + " for reading."};

    auto result = read_stockholm_file<alphabet_type>(stream);
    stream.close();
    return std::move(result);
}

} // namespace mars
