#include <chrono>
#include <fstream>
#include <future>
#include <iostream>
#include <thread>
#include <vector>

#include "index.hpp"
#include "motif.hpp"
#include "search.hpp"
#include "settings.hpp"

int main(int argc, char ** argv)
{
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    // Parse arguments
    mars::Settings settings{};
    if (!settings.parse_arguments(argc, argv))
        return EXIT_FAILURE;

    // Start reading the genome and creating the index asyncronously
    mars::BiDirectionalIndex bds{settings.xdrop};
    std::future<void> index_future = std::async(std::launch::async, &mars::BiDirectionalIndex::create, &bds,
                                                settings.genome_file);

    // Generate motifs from the MSA
    std::vector<mars::StemloopMotif> motifs = mars::create_motifs(settings.alignment_file, settings.threads);

    std::thread write_rssp([&motifs] (std::filesystem::path const & file)
    {
        if (!motifs.empty() && !file.empty())
        {
            if (mars::verbose > 0)
                std::cerr << "Exporting the stem loops in rssp format ==> " << file << std::endl;
            std::ofstream os(file);
            for (auto const & motif : motifs)
                motif.print_rssp(os);
            os.close();
        }
    }, settings.structator_file);

    // Wait for index creation process
    try
    {
        index_future.get();
    }
    catch (std::exception const & e)
    {
        std::cerr << "EXCEPTION => " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    if (!motifs.empty() && !settings.genome_file.empty())
    {
        mars::SearchGenerator search{bds, motifs.front().depth};
        search.find_motifs(motifs, settings.threads, settings.min_score_per_motif);

        auto print_results = [&bds, &search] (std::ostream & out)
        {
            if (!search.get_locations().empty())
                out << " " << std::left << std::setw(35) << "sequence name" << "\t" << "index" << "\t"
                    << "pos" << "\t" << "n" << "\t" << "score" << std::endl;
            for (mars::MotifLocation const & loc : search.get_locations())
                out << ">" << std::left << std::setw(35) << bds.get_name(loc.sequence) << "\t" << loc.sequence << "\t"
                    << loc.position << "\t" << +loc.num_stemloops << "\t" << loc.score << std::endl;
        };

        if (!settings.result_file.empty())
        {
            if (mars::verbose > 0)
                std::cerr << "Writing results ==> " << settings.result_file << std::endl;
            std::ofstream file_stream(settings.result_file);
            print_results(file_stream);
            file_stream.close();
        }
        else
        {
            print_results(std::cout);
        }
    }
    else if (motifs.empty() && mars::verbose > 0)
    {
        std::cerr << "There are no motifs: skipping search step." << std::endl;
    }
    else if (mars::verbose > 0)
    {
        std::cerr << "No genome provided: skipping search step." << std::endl;
    }
    write_rssp.join();

    if (mars::verbose > 0)
    {
        auto sec = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0).count();
        std::cerr << argv[0] << " has finished after " << sec << " seconds." << std::endl;
    }

    return 0;
}
