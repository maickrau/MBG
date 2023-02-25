#include <vector>
#include <string>
#include <cxxopts.hpp>
#include "MBG.h"

int main(int argc, char** argv)
{
	std::cerr << "MBG " << VERSION << std::endl;
	cxxopts::Options options { "MBG" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("i,in", "Input reads. Multiple files can be input with -i file1.fa -i file2.fa etc (required)", cxxopts::value<std::vector<std::string>>())
		("o,out", "Output graph (required)", cxxopts::value<std::string>())
		("t", "Number of threads", cxxopts::value<size_t>()->default_value("1"))
		("k", "K-mer size. Must be odd and >=31 (required)", cxxopts::value<size_t>())
		("w", "Window size. Must be 1 <= w <= k-30 (default: k-30)", cxxopts::value<size_t>())
		("a,kmer-abundance", "Minimum k-mer abundance", cxxopts::value<size_t>()->default_value("1"))
		("u,unitig-abundance", "Minimum average unitig abundance", cxxopts::value<double>()->default_value("2"))
		("error-masking", "Error masking", cxxopts::value<std::string>()->default_value("hpc"))
		("blunt", "Output a bluntified graph without edge overlaps")
		("include-end-kmers", "Force k-mers at read ends to be included")
		("output-sequence-paths", "Output the paths of the input sequences to a file (.gaf)", cxxopts::value<std::string>())
		("r,resolve-maxk", "Maximum k-mer size for multiplex DBG resolution", cxxopts::value<size_t>())
		("R,resolve-maxk-allowgaps", "Allow multiplex resolution to add gaps up to this k-mer size", cxxopts::value<size_t>())
		("node-name-prefix", "Add a prefix to output node names", cxxopts::value<std::string>())
		("sequence-cache-file", "Use a temporary sequence cache file to speed up graph construction", cxxopts::value<std::string>())
		("keep-gaps", "Don't remove low coverage nodes if it would leave a gap in the graph")
		("hpc-variant-onecopy-coverage", "Separate k-mers based on hpc variants, using arg as single copy coverage", cxxopts::value<double>())
		("do-unsafe-guesswork-resolutions", "Use extra heuristics during multiplex resolution")
		("copycount-filter-heuristic", "Use coverage based heuristic filter during multiplex resolution")
		("only-local-resolve", "Only resolve nodes which are repetitive within a read")
		("output-homology-map", "Output a list of homologous k-mer locations", cxxopts::value<std::string>())
		("no-kmer-filter-inside-unitig", "Don't filter out k-mers which are completely contained by two other k-mers")
		("no-multiplex-cleaning", "Don't clean low coverage tips and structures during multiplex resolution")
	;
	auto params = options.parse(argc, argv);
	if (params.count("v") == 1)
	{
		std::cerr << "Version: " << VERSION << std::endl;
		exit(0);
	}
	if (params.count("h") == 1)
	{
		std::cerr << options.help();
		std::cerr << "Options for --error-masking:" << std::endl;
		std::cerr << "\tno\tNo error masking" << std::endl;
		std::cerr << "\thpc\tMask homopolymer errors (default)" << std::endl;
		std::cerr << "\tdinuc\tMask homopolymer and dinucleotide errors" << std::endl;
		std::cerr << "\tmsat\tMask homopolymer and microsatellite errors up to 6bp" << std::endl;
		std::cerr << "\tcollapse\tCollapse homopolymers" << std::endl;
		std::cerr << "\tcollapse-dinuc\tCollapse homopolymers and mask dinucleotide errors" << std::endl;
		std::cerr << "\tcollapse-msat\tCollapse homopolymers and mask microsatellite errors up to 6bp" << std::endl;
		exit(0);
	}
	bool paramError = false;
	if (params.count("i") == 0)
	{
		std::cerr << "Select input reads -i" << std::endl;
		paramError = true;
	}
	if (params.count("o") == 0)
	{
		std::cerr << "Select output filename -o" << std::endl;
		paramError = true;
	}
	if (params.count("k") == 0)
	{
		std::cerr << "Select k-mer size -k" << std::endl;
		paramError = true;
	}
	if (paramError) std::exit(1);
	std::vector<std::string> inputReads = params["i"].as<std::vector<std::string>>();
	std::string outputGraph = params["o"].as<std::string>();
	size_t kmerSize = params["k"].as<size_t>();
	size_t windowSize = 1;
	size_t maxResolveLength = 0;
	size_t maxUnconditionalResolveLength = 0;
	if (params.count("w") == 1)
	{
		windowSize = params["w"].as<size_t>();
	}
	else
	{
		windowSize = kmerSize - 30;
	}
	size_t minCoverage = params["a"].as<size_t>();
	double minUnitigCoverage = params["u"].as<double>();
	size_t numThreads = params["t"].as<size_t>();
	std::string outputSequencePaths = "";
	ErrorMasking errorMasking = ErrorMasking::Hpc;
	bool blunt = false;
	bool includeEndKmers = false;
	bool keepGaps = false;
	bool guesswork = false;
	double hpcVariantOnecopyCoverage = 0;
	bool copycountFilterHeuristic = false;
	bool onlyLocalResolve = false;
	bool filterWithinUnitig = true;
	bool doCleaning = true;
	std::string errorMaskingStr = "hpc";
	std::string nodeNamePrefix = "";
	std::string sequenceCacheFile = "";
	std::string outputHomologyMap = "";
	if (params.count("r") == 1) maxResolveLength = params["r"].as<size_t>();
	if (params.count("R") == 1) maxUnconditionalResolveLength = params["R"].as<size_t>();
	if (params.count("blunt") == 1) blunt = true;
	if (params.count("keep-gaps") == 1) keepGaps = true;
	if (params.count("do-unsafe-guesswork-resolutions") == 1) guesswork = true;
	if (params.count("error-masking") == 1)
	{
		errorMaskingStr = params["error-masking"].as<std::string>();
		if (params["error-masking"].as<std::string>() == "no")
		{
			errorMasking = ErrorMasking::No;
		}
		else if (params["error-masking"].as<std::string>() == "hpc")
		{
			errorMasking = ErrorMasking::Hpc;
		}
		else if (params["error-masking"].as<std::string>() == "dinuc")
		{
			errorMasking = ErrorMasking::Dinuc;
		}
		else if (params["error-masking"].as<std::string>() == "collapse")
		{
			errorMasking = ErrorMasking::Collapse;
		}
		else if (params["error-masking"].as<std::string>() == "msat")
		{
			errorMasking = ErrorMasking::Microsatellite;
		}
		else if (params["error-masking"].as<std::string>() == "collapse-dinuc")
		{
			errorMasking = ErrorMasking::CollapseDinuc;
		}
		else if (params["error-masking"].as<std::string>() == "collapse-msat")
		{
			errorMasking = ErrorMasking::CollapseMicrosatellite;
		}
		else
		{
			std::cerr << "unknown parameter for --error-masking: \"" << params["error-masking"].as<std::string>() << "\"" << std::endl;
			paramError = true;
		}
	}
	if (params.count("include-end-kmers") == 1) includeEndKmers = true;
	if (params.count("output-sequence-paths") == 1) outputSequencePaths = params["output-sequence-paths"].as<std::string>();
	if (params.count("node-name-prefix") == 1) nodeNamePrefix = params["node-name-prefix"].as<std::string>();
	if (params.count("sequence-cache-file") == 1) sequenceCacheFile = params["sequence-cache-file"].as<std::string>();
	if (params.count("hpc-variant-onecopy-coverage") == 1) hpcVariantOnecopyCoverage = params["hpc-variant-onecopy-coverage"].as<double>();
	if (params.count("copycount-filter-heuristic") == 1) copycountFilterHeuristic = true;
	if (params.count("only-local-resolve") == 1) onlyLocalResolve = true;
	if (params.count("output-homology-map") == 1) outputHomologyMap = params["output-homology-map"].as<std::string>();
	if (params.count("no-kmer-filter-inside-unitig") == 1) filterWithinUnitig = false;
	if (params.count("no-multiplex-cleaning")) doCleaning = false;

	if (numThreads == 0)
	{
		std::cerr << "Number of threads cannot be 0" << std::endl;
		paramError = true;
	}
	if (params.count("w") == 1 && windowSize > kmerSize)
	{
		std::cerr << "Window size cannot be greater than k-mer size" << std::endl;
		paramError = true;
	}
	if (params.count("w") == 1 && windowSize == 0)
	{
		std::cerr << "Window size must be >0" << std::endl;
		paramError = true;
	}
	if (kmerSize % 2 == 0)
	{
		std::cerr << "K-mer size must be odd" << std::endl;
		paramError = true;
	}
	if (kmerSize < 31)
	{
		std::cerr << "Minimum k-mer size is 31" << std::endl;
		paramError = true;
	}
	if (params.count("w") == 1 && windowSize > kmerSize - 30)
	{
		std::cerr << "Window size must be <= k-30. With current k (" << kmerSize << ") maximum w is " << kmerSize - 30 << std::endl;
		paramError = true;
	}
	if (outputSequencePaths != "" && blunt)
	{
		std::cerr << "--output-sequence-paths and --blunt are not supported together" << std::endl;
		paramError = true;
	}
	if (maxResolveLength > 0 && blunt)
	{
		std::cerr << "-r (--resolve-maxk) and --blunt are not supported together" << std::endl;
		paramError = true;
	}
	if (paramError) std::abort();
	
	std::cerr << "Parameters: ";
	std::cerr << "k=" << kmerSize << ",";
	std::cerr << "w=" << windowSize << ",";
	std::cerr << "a=" << minCoverage << ",";
	std::cerr << "u=" << minUnitigCoverage << ",";
	std::cerr << "t=" << numThreads << ",";
	std::cerr << "r=" << maxResolveLength << ",";
	std::cerr << "R=" << maxUnconditionalResolveLength << ",";
	std::cerr << "hpcvariantcov=" << hpcVariantOnecopyCoverage << ",";
	std::cerr << "errormasking=" << errorMaskingStr << ",";
	std::cerr << "endkmers=" << (includeEndKmers ? "yes" : "no") << ",";
	std::cerr << "blunt=" << (blunt ? "yes" : "no") << ",";
	std::cerr << "keepgaps=" << (keepGaps ? "yes" : "no") << ",";
	std::cerr << "guesswork=" << (guesswork ? "yes" : "no") << ",";
	std::cerr << "copycountfilter=" << (copycountFilterHeuristic ? "yes" : "no") << ",";
	std::cerr << "onlylocal=" << (onlyLocalResolve ? "yes" : "no") << ",";
	std::cerr << "filterwithinunitig=" << (filterWithinUnitig ? "yes" : "no") << ",";
	std::cerr << "cleaning=" << (doCleaning ? "yes" : "no") << ",";
	std::cerr << "cache=" << (sequenceCacheFile.size() > 0 ? "yes" : "no");
	std::cerr << std::endl;

	runMBG(inputReads, outputGraph, kmerSize, windowSize, minCoverage, minUnitigCoverage, errorMasking, numThreads, includeEndKmers, outputSequencePaths, maxResolveLength, blunt, maxUnconditionalResolveLength, nodeNamePrefix, sequenceCacheFile, keepGaps, hpcVariantOnecopyCoverage, guesswork, copycountFilterHeuristic, onlyLocalResolve, outputHomologyMap, filterWithinUnitig, doCleaning);
}