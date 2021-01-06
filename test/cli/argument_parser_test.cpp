#include <string>                // strings

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("mars");
    std::string expected
    {
        "mars - Motif-based aligned RNA searcher\n"
        "=======================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

//TEST_F(cli_test, fail_no_argument)
//{
//    cli_test_result result = execute_app("mars", "-o", data("out.aln"));
//    std::string expected
//    {
//        "Parsing error. Not enough positional arguments provided (Need at least 1). "
//        "See -h/--help for more information.\n"
//    };
//    EXPECT_NE(result.exit_code, 0);
//    EXPECT_EQ(result.out, std::string{});
//    EXPECT_EQ(result.err, expected);
//}

//TEST_F(cli_test, with_argument)
//{
//    cli_test_result result = execute_app("mars", data("tRNA.aln"));
//    EXPECT_EQ(result.exit_code, 0);
//    EXPECT_EQ(result.out, "> seq1\nACGTTTGATTCGCG\n> seq2\nTCGGGGGATTCGCG\n");
//    EXPECT_EQ(result.err, std::string{});
//}
//
//TEST_F(cli_test, with_argument_verbose)
//{
//    cli_test_result result = execute_app("mars", data("tRNA.aln"), "-v");
//    EXPECT_EQ(result.exit_code, 0);
//    EXPECT_EQ(result.out, "> seq1\nACGTTTGATTCGCG\n> seq2\nTCGGGGGATTCGCG\n");
//    EXPECT_EQ(result.err, "Conversion was a success. Congrats!\n");
//}
//
//TEST_F(cli_test, with_out_file)
//{
//    cli_test_result result = execute_app("mars", data("tRNA.aln"), "-o", "out.fasta");
//    seqan3::sequence_file_input fin{"out.fasta", seqan3::fields<seqan3::field::seq, seqan3::field::id>{}};
//
//    // create records to compare
//    using record_type = typename decltype(fin)::record_type;
//    using seqan3::operator""_dna5;
//    std::vector<record_type> records{};
//    records.emplace_back("ACGTTTGATTCGCG"_dna5, std::string{"seq1"});
//    records.emplace_back("TCGGGGGATTCGCG"_dna5, std::string{"seq2"});
//
//    EXPECT_RANGE_EQ(fin, records);
//    EXPECT_EQ(result.exit_code, 0);
//    EXPECT_EQ(result.out, std::string{});
//    EXPECT_EQ(result.err, std::string{});
//}
