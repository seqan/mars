#pragma once

#include <cmath>
#include <seqan3/std/filesystem>
#include <memory>

#include <seqan3/core/debug_stream.hpp>

#include "ThreadPool.hpp"

#define logger(vlevel, lstr)                                    \
{                                                               \
    if (mars::settings.verbose >= (vlevel))                     \
    {                                                           \
        std::lock_guard<std::mutex> guard(mars::mutex_console); \
        seqan3::debug_stream << lstr;                           \
    }                                                           \
}

namespace mars
{

struct Settings
{
    // input
    std::filesystem::path genome_file{};
    std::filesystem::path alignment_file{};
    // output
    std::filesystem::path result_file{};
    std::filesystem::path motif_file{};
    std::filesystem::path structator_file{};
    float min_score_per_motif{NAN};
    unsigned short verbose{1};
    // performance
    unsigned char prune{10};
    unsigned char xdrop{4};
    bool exterior{true};
    bool compress_index{false};
    unsigned int nthreads{std::thread::hardware_concurrency()};

    bool parse_arguments(int argc, char ** argv);
};

extern std::unique_ptr<thread_pool::ThreadPool> pool;
extern std::mutex mutex_console;
extern Settings settings;

} // namespace mars
