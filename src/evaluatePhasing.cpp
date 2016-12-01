/*
============================================================================
Phase BAM using phased SNP scaffold
============================================================================
Copyright (C) 2016 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#include <iostream>
#include <fstream>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/vcf.h>

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include "util.h"


using namespace phasebam;

struct Config {
  bool filterForPass;
  int32_t minac;
  std::string chrom;
  std::string sample;
  boost::filesystem::path outfile;
  boost::filesystem::path p1;
  boost::filesystem::path p2;
};


template<typename TConfig>
inline int
evaluatePhasing(TConfig& c) {

#ifdef PROFILE
  ProfilerStart("phasebam.prof");
#endif

  // Load variants
  typedef std::vector<Variant> TPhasedVariants;
  TPhasedVariants pV1;
  TPhasedVariants pV2;
  if (!_loadVariants(c.sample, c.chrom, c.p1.string(), c.filterForPass, c.minac, true, pV1)) return -1;
  if (!_loadVariants(c.sample, c.chrom, c.p2.string(), false, 0, true, pV2)) return -1;

  // Find common variants
  TPhasedVariants hap1;
  TPhasedVariants hap2;
  int32_t commonSites = 0;
  int32_t distinctSites1 = 0;
  int32_t distinctSites2 = 0;
  uint32_t pos1 = 0;
  uint32_t pos2 = 0;
  while ((pos1 < pV1.size()) && (pos2 < pV2.size())) {
    if (pV1[pos1].pos != pV2[pos2].pos) {
      if (pV1[pos1].pos < pV2[pos2].pos) {
	++distinctSites1;
	++pos1;
      } else {
	++distinctSites2;
	++pos2;
      }
    } else { 
      if ((pV1[pos1].ref == pV2[pos2].ref) && (pV1[pos1].alt == pV2[pos2].alt)) {
	++commonSites;
	hap1.push_back(pV1[pos1]);
	hap2.push_back(pV2[pos2]);
      } else {
	++distinctSites1;
	++distinctSites2;
      }
      ++pos1;
      ++pos2;
    }
  }

  // Hamming count
  int32_t hamcount = 0;
  for(uint32_t i=0; i<hap1.size();++i)
    if (hap1[i].hap != hap2[i].hap) ++hamcount;

  // Get differences
  typedef std::vector<bool> THap;
  int32_t diff1size = 0;
  if (!hap1.empty()) diff1size = hap1.size() - 1;
  THap diff1(diff1size, 0);
  for(uint32_t i=1; i<hap1.size();++i)
    if (hap1[i-1].hap != hap1[i].hap) diff1[i-1]=1;
  int32_t diff2size = 0;
  if (!hap2.empty()) diff2size = hap2.size() - 1;
  THap diff2(diff2size, 0);
  for(uint32_t i=1; i<hap2.size();++i)
    if (hap2[i-1].hap != hap2[i].hap) diff2[i-1]=1;

  // Output switch errors
  boost::iostreams::filtering_ostream dataOut;
  dataOut.push(boost::iostreams::gzip_compressor());
  dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
  dataOut << "distance\tacPre\tacSuc" << std::endl;

  int32_t lastPos = -1;
  int32_t switchcount = 0;
  for(uint32_t i=0; i<diff1.size(); ++i) {
    if (diff1[i] != diff2[i]) {
      int32_t distance = -1;
      if (lastPos != -1) distance = (hap1[i].pos - lastPos);
      lastPos = hap1[i+1].pos;
      if (distance != -1) dataOut << distance << "\t";
      else dataOut << "NA\t";
      if ((hap1[i].ac != -1) && (hap1[i+1].ac != -1)) dataOut << hap1[i].ac << '\t' << hap1[i+1].ac << std::endl;
      else dataOut << "NA\tNA" << std::endl;
      ++switchcount;
    }
  }


  // Switch count
  std::cout << "Sample\tChromosome\tCommonSites\tDistinctBCF1\tDistinctBCF2\tSwitchErrors\tSwitchErrorRate\tHammingDistance\tHammingErrorRate" << std::endl;
  std::cout << c.sample << "\t" << c.chrom << "\t" << commonSites << "\t" << distinctSites1 << "\t" << distinctSites2 << "\t" << switchcount << "\t" << (double) switchcount / (double) diff1.size() << "\t" << hamcount << "\t" << (double) hamcount / (double) hap1.size() << std::endl;

  
#ifdef PROFILE
  ProfilerStop();
#endif


  return 0;
}

int main(int argc, char **argv) {
  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("phased,p", boost::program_options::value<boost::filesystem::path>(&c.p1), "gold-standard phased BCF")
    ("sample,s", boost::program_options::value<std::string>(&c.sample)->default_value("HG00512"), "sample name")
    ("chrom,c", boost::program_options::value<std::string>(&c.chrom)->default_value("1"), "chromosome name")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("swerr.gz"), "switch error")
    ("minac,m", boost::program_options::value<int32_t>(&c.minac)->default_value(25), "min. allele count in gold-standard phased BCF")
    ("pass,a", "Filter sites for PASS in gold-standard phased BCF")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.p2), "input phased BCF")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("phased"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -s NHG00512 -p <phased.bcf> -c 1 <other.phased.bcf>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  } 

  // Check VCF/BCF file
  if (vm.count("input-file")) {
    if (!(boost::filesystem::exists(c.p1) && boost::filesystem::is_regular_file(c.p1) && boost::filesystem::file_size(c.p1))) {
      std::cerr << "Input SNP VCF/BCF file is missing: " << c.p1.string() << std::endl;
      return 1;
    }
    htsFile* ifile = bcf_open(c.p1.string().c_str(), "r");
    if (ifile == NULL) {
      std::cerr << "Fail to open file " << c.p1.string() << std::endl;
      return 1;
    }
    hts_idx_t* bcfidx = bcf_index_load(c.p1.string().c_str());
    if (bcfidx == NULL) {
      std::cerr << "Fail to open index file for " << c.p1.string() << std::endl;
      return 1;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    if (hdr == NULL) {
      std::cerr << "Fail to open header for " << c.p1.string() << std::endl;
      return 1;
    }
    bcf_hdr_destroy(hdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ifile);
  }

  // Check VCF/BCF file
  if (vm.count("input-file")) {
    if (!(boost::filesystem::exists(c.p2) && boost::filesystem::is_regular_file(c.p2) && boost::filesystem::file_size(c.p2))) {
      std::cerr << "Input SNP VCF/BCF file is missing: " << c.p2.string() << std::endl;
      return 1;
    }
    htsFile* ifile = bcf_open(c.p2.string().c_str(), "r");
    if (ifile == NULL) {
      std::cerr << "Fail to open file " << c.p2.string() << std::endl;
      return 1;
    }
    hts_idx_t* bcfidx = bcf_index_load(c.p2.string().c_str());
    if (bcfidx == NULL) {
      std::cerr << "Fail to open index file for " << c.p2.string() << std::endl;
      return 1;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    if (hdr == NULL) {
      std::cerr << "Fail to open header for " << c.p2.string() << std::endl;
      return 1;
    }
    bcf_hdr_destroy(hdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ifile);
  }

  // Filter for PASS
  if (vm.count("pass")) c.filterForPass = true;
  else c.filterForPass = false;
  
  return evaluatePhasing(c);
}
