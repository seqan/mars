#pragma once

#include <seqan3/std/filesystem>


namespace mars
{

struct Settings
{
private:
    std::filesystem::path result_file{};

public:
    std::filesystem::path alignment_file{};
    std::filesystem::path genome_file{};
    unsigned char xdrop{4};

    bool parse_arguments(int argc, char ** argv, std::ostream & out);
};

} // namespace mars
