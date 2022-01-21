#pragma once

#include <seqan3/std/filesystem>
#include <memory>

#include <seqan3/core/debug_stream.hpp>

#include "ThreadPool.hpp"

#define logger(vlevel, lstr)                                 \
{                                                            \
    if (mars::settings.verbose >= (vlevel))                  \
    {                                                        \
        std::lock_guard<std::mutex> guard(mars::mutex_cerr); \
        seqan3::debug_stream << lstr;                        \
    }                                                        \
}

namespace mars
{

struct Settings
{
    std::filesystem::path alignment_file{};
    std::filesystem::path genome_file{};
    std::filesystem::path motif_file{};
    std::filesystem::path result_file{};
    std::filesystem::path structator_file{};
    unsigned char xdrop{4};
    float min_score_per_motif{5.f};
    unsigned char prune{10};
    bool compress_index{false};
    unsigned short verbose{1};

    bool parse_arguments(int argc, char ** argv);
};

extern std::unique_ptr<thread_pool::ThreadPool> pool;
extern std::mutex mutex_cerr;
extern Settings settings;

} // namespace mars
