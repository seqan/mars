#pragma once

#include <tuple>
#include <vector>

#include <seqan3/std/filesystem>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

namespace mars
{

using index_type = seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::collection>;

index_type create_index(std::filesystem::path const & filepath);

class bi_directional_search
{
private:
    index_type const & index;
    std::vector<seqan3::dna4_vector> queries{};

public:
    std::vector<std::pair<size_t, size_t>> matches{};

    explicit bi_directional_search(index_type const & bi_dir_index): index(bi_dir_index)
    {
        queries.emplace_back();
    }

    void append_5prime(seqan3::dna4 left);
    void append_3prime(seqan3::dna4 right);
    void append_stem(seqan3::dna4 left, seqan3::dna4 right);
    bool backtrack();

    size_t compute_matches();
};

} // namespace mars
