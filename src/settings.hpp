#pragma once

#include <seqan3/std/filesystem>
#include <fstream>

namespace mars
{

extern unsigned short verbose;

struct Settings
{
private:
    std::filesystem::path result_file{};
    std::ofstream file_stream{};

public:
    std::filesystem::path alignment_file{};
    std::filesystem::path genome_file{};
    unsigned char xdrop{4};
    unsigned int threads{1};
    float min_score_per_motif{5.f};

    bool parse_arguments(int argc, char ** argv, std::ostream & out);
};

} // namespace mars
