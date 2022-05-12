// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

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

//! \brief Settings for the Mars program.
struct Settings
{
    // input
    std::filesystem::path genome_file{}; //!< The filename for reading the genome.
    std::filesystem::path alignment_file{}; //!< The filename for reading the alignment.
    // output
    std::filesystem::path result_file{}; //!< The filename for writing the results (locations).
    std::filesystem::path motif_file{}; //!< The filename for writing the motifs.
    std::filesystem::path structator_file{}; //!< The filename for writing the Structator RSSPs.
    float score_filter{NAN}; //!< The minimum score per stemloop for the output, NAN = evalue criterion.
    unsigned short verbose{1}; //!< The verbosity level of the output.
    // performance
    unsigned char prune{10}; //!< Parameter for reducing the motif.
    unsigned char xdrop{4};  //!< Parameter for pruning the search.
    bool limit{false}; //!< Flag whether exterior loops are considered.
    bool compress_index{false}; //!< Flag whether the index should be compressed.
    unsigned int nthreads{std::thread::hardware_concurrency()};  //!< The number of threads in the pool.

    /*!
     * \brief Run the argument parser.
     * \param argc The arguments from the program call.
     * \param argv The number of arguments.
     * \return whether parsing was successful.
     */
    bool parse_arguments(int argc, char ** argv);
};

//! \brief The thread pool for scheduling asynchronous tasks.
extern std::unique_ptr<thread_pool::ThreadPool> pool;

//! \brief A mutex for concurrent console output.
extern std::mutex mutex_console;

//! \brief Public settings object to be accessed from anywhere in Mars.
extern Settings settings;

} // namespace mars
