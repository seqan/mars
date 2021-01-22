#pragma once

#include <seqan3/std/filesystem>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

namespace mars
{

using index_type = seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::collection>;

index_type create_index(std::filesystem::path const & filepath);

} // namespace mars
