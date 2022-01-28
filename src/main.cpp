#include <chrono>
#include <vector>

#include "search.hpp"
#include "settings.hpp"

int main(int argc, char ** argv)
{
    std::chrono::steady_clock::time_point tm0 = std::chrono::steady_clock::now();

    // Parse arguments
    if (!mars::settings.parse_arguments(argc, argv))
        return EXIT_FAILURE;

    // Start reading the genome and creating the index asyncronously
    mars::BiDirectionalIndex index{};
    auto future_index = mars::pool->submit(&mars::BiDirectionalIndex::create, &index);

    // Generate motifs from the MSA
    std::vector<mars::StemloopMotif> motifs = mars::create_motifs();
    auto future_json = mars::pool->submit(mars::store_motifs, motifs);
    auto future_rssp = mars::pool->submit(mars::store_rssp, motifs);

    // Wait for index creation process
    future_index.wait();

    if (!motifs.empty() && !index.raw().empty())
    {
        // Search the genome for motifs
        auto locations = mars::find_motifs(index, motifs);
        mars::print_locations(locations, index);
    }
    else if (motifs.empty())
    {
        logger(1, "There are no motifs: skipping search step." << std::endl);
    }
    else
    {
        logger(1, "No genome sequence provided: skipping search step." << std::endl);
    }
    future_json.wait();
    future_rssp.wait();

    // print run time
    auto const sec = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - tm0).count();
    logger(1, argv[0] << " has finished after " << sec << " seconds." << std::endl);

    return 0;
}
