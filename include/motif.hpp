#pragma once

#include <variant>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>

#include "bi_alphabet.hpp"
#include "multiple_alignment.hpp"
#include "profile_sequence.hpp"

namespace mars
{

//!\brief Store the type of statistics values.
typedef uint32_t value_type;

//!\brief Store {start, end} positions of an inverval.
typedef std::pair<value_type, value_type> interval_type;

//!\brief Store {min, median, max} of a distribution.
typedef std::tuple<value_type, value_type, value_type> stat_type;

struct secondary_structure
{
    using stem_t = profile_sequence<bi_alphabet<seqan3::rna4>>;
    using loop_t = profile_sequence<seqan3::rna4>;
    bool loop_left;
    std::variant<stem_t, loop_t> profile;
    stat_type length;
    // TODO: gaps
};

struct stem_loop_motif
{
    interval_type bounds;
    stat_type length;
    std::vector<secondary_structure> elements;
};

using stemloop_type = std::vector<std::pair<int, int>>;
stemloop_type detect_stem_loops(std::vector<int> const & bpseq, std::vector<int> const & plevel);

stem_loop_motif analyze_stem_loop(msa_type const & msa, std::vector<int> const & bpseq, std::pair<int, int> const & pos);

} // namespace mars
