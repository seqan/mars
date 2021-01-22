#pragma once

#include <variant>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>

#include "bi_alphabet.hpp"
#include "multiple_alignment.hpp"
#include "profile_char.hpp"

namespace mars
{

//!\brief Store the type of statistics values.
typedef uint32_t value_type;

//!\brief Store {start, end} positions of an inverval.
typedef std::pair<value_type, value_type> interval_type;

//!\brief Store {min, median, max} of a distribution.
typedef std::tuple<value_type, value_type, value_type> stat_type;

struct loop_element
{
    stat_type length;
    std::vector<profile_char<seqan3::rna4>> profile;
    std::vector<std::unordered_map<uint16_t, uint16_t>> gaps;
    bool is_5prime;
};

struct stem_element
{
    stat_type length;
    std::vector<profile_char<bi_alphabet<seqan3::rna4>>> profile;
    std::vector<std::unordered_map<uint16_t, uint16_t>> gaps;
};

struct stem_loop_motif
{
    using element_type = std::variant<loop_element, stem_element>;

    interval_type bounds;
    stat_type length;
    std::vector<element_type> elements;

    stem_element & new_stem()
    {
        elements.emplace_back<stem_element>({});
        return std::get<stem_element>(elements.back());
    }

    loop_element & new_loop()
    {
        elements.emplace_back<loop_element>({});
        return std::get<loop_element>(elements.back());
    }
};

/*!
 * \brief Stream a representation of a stem loop motif.
 * \tparam ostream_type The type of the output stream.
 * \param os The output stream.
 * \param motif The motif to be printed.
 * \return The stream with the motif representation appended.
 */
template <typename ostream_type>
inline ostream_type & operator<<(ostream_type & os, stem_loop_motif const & motif)
{
    for (auto const & el : motif.elements)
    {
        std::visit([&os] (auto element)
        {
            if constexpr (std::is_same_v<decltype(element), mars::loop_element>)
               os << "Loop " << (element.is_5prime ? "5' " : "3' ");
            else
               os << "Stem ";

            for (auto const & profile_char : element.profile)
                os << profile_char << ' ';
            os << "\nGaps: ";
            for (int idx = 0; idx < element.gaps.size(); ++idx)
                if (!element.gaps[idx].empty())
                {
                    os << "\t" << idx << ": ";
                    for (auto && [key, val] : element.gaps[idx])
                        os << "(" << key << "," << val << ")";
                }
            os << "\n";
        }, el);
    }
    return os;
}

using stemloop_type = std::vector<std::pair<int, int>>;
stemloop_type detect_stem_loops(std::vector<int> const & bpseq, std::vector<int> const & plevel);

stem_loop_motif analyze_stem_loop(msa_type const & msa,
                                  std::vector<int> const & bpseq,
                                  std::pair<int, int> const & pos);

} // namespace mars
