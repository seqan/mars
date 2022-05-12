// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <array>
#include <mutex>
#include <ostream>
#include <string>
#include <vector>

namespace mars
{

//! \brief A MotifLocation is a region with many high-scoring StemloopHits.
struct MotifLocation
{
    double evalue; //!< Likelihood that the location is significant.
    float score; //!< Bit-score of the found location.
    uint8_t num_stemloops; //!< The number of found stemloops at this location.
    size_t position_start; //!< The start position of this location in the genome.
    size_t position_end; //!< The end position of this location in the genome.
    size_t query_length; //!< The total length of the individual stemloop hits.
    size_t sequence; //!< The sequence number within the genome.
};

/*!
 * \brief Comparison for MotifLocation.
 * \param lhs The left-hand-side for the comparison.
 * \param rhs The right-hand-side for the comparison.
 * \return whether lhs is smaller than rhs, based on evalue or bitscore.
 */
bool operator<(MotifLocation const & lhs, MotifLocation const & rhs);

//! \brief A class that stores MotifLocations and prints them in order.
class MotifLocationStore : public std::vector<MotifLocation>
{
private:
    //! \brief A reference to the sequence names for printing.
    std::vector<std::string> const & names;

    //! \brief The mutex for concurrent pushing.
    std::mutex mutex_locations;

    /*!
     * \brief Print the collected motifs, preceeded by a header line.
     * \param out The output stream.
     */
    void print(std::ostream & out);

public:
    /*!
     * \brief Constructor with a reference to the names vector.
     * \param names The sequence names.
     */
    explicit MotifLocationStore(std::vector<std::string> const & names) : names{names} {}

    //! \brief Sort all the locations and print them in order.
    void print();

    /*!
     * \brief Add a location to the collection.
     * \param loc The location to be stored.
     */
    void push(MotifLocation && loc);
};

//! \brief A StemloopHit is a genome position where a stemloop matches.
struct StemloopHit
{
    long long pos; //!< The start position within the genome sequence.
    long long length; //!< The length of the stemloop match.
    uint8_t midx; //!< The id of the matching stemloop.
    float score; //!< The score of the match.
};

/*!
 * \brief Comparison for StemloopHit.
 * \param lhs The left-hand-side for the comparison.
 * \param rhs The right-hand-side for the comparison.
 * \return whether lhs is smaller than rhs, based on position.
 */
bool operator<(StemloopHit const & lhs, StemloopHit const & rhs);

//! \brief A class that stores StemloopHits separately for each genome sequence.
class StemloopHitStore
{
private:
    std::vector<std::vector<StemloopHit>> hits; //!< The storage for hits.
    std::array<std::mutex, 256> mutexes; //!< Mutexes for concurrent pushing.

public:
    /*!
     * \brief Constructor which allocates a vector of size seq_count for storing hits.
     * \param seq_count The number of genome sequences in the index.
     */
    explicit StemloopHitStore(size_t seq_count) : hits(seq_count) {}

    /*!
     * \brief Add a hit to the collection.
     * \param hit The stemloop match.
     * \param seq The sequence where the match is located.
     */
    void push(StemloopHit && hit, size_t seq);

    /*!
     * \brief Retrieve the hits for one sequence.
     * \param seq The sequence id.
     * \return a reference to the specified hit vector.
     */
    std::vector<StemloopHit> & get(size_t seq);
};

} // namespace mars
