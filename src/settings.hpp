#pragma once

#include <seqan3/std/filesystem>

namespace mars
{

extern unsigned short verbose;

struct Settings
{
    std::filesystem::path alignment_file{};
    std::filesystem::path genome_file{};
    std::filesystem::path result_file{};
    std::filesystem::path structator_file{};
    unsigned char xdrop{4};
    unsigned int threads{1};
    float min_score_per_motif{5.f};

    bool parse_arguments(int argc, char ** argv);
};

} // namespace mars
