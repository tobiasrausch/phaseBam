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

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

struct Config {
  std::string chrom;
  std::string sample;
  boost::filesystem::path p1;
  boost::filesystem::path p2;
};

struct Variant {
  int32_t pos;
  std::string ref;
  std::string alt;
  bool hap;

  Variant(int32_t p, std::string r, std::string a, bool h) : pos(p), ref(r), alt(a), hap(h) {}
};


template<typename TConfig, typename TPhasedVariants>
inline bool
_loadVariants(TConfig const& c, std::string const& bcffile, TPhasedVariants& pV) {
  typedef typename TPhasedVariants::value_type TVariant;

  // Load bcf file
  htsFile* ifile = hts_open(bcffile.c_str(), "r");
  hts_idx_t* bcfidx = bcf_index_load(bcffile.c_str());
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);
  int32_t sampleIndex = -1;
  for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
    if (hdr->samples[i] == c.sample) sampleIndex = i;
  if (sampleIndex < 0) {
    std::cerr << "Sample not found " << c.sample << std::endl;
    return false;
  }

  // Genotypes
  int ngt = 0;
  int32_t* gt = NULL;
  
  // Collect Snps for this chromosome
  int32_t chrid = bcf_hdr_name2id(hdr, c.chrom.c_str());
  if (chrid < 0) {
    std::cerr << "Chromosome not found " << c.chrom << std::endl;
    return false;
  }
  hts_itr_t* itervcf = bcf_itr_querys(bcfidx, hdr, c.chrom.c_str());
  if (itervcf != NULL) {
    bcf1_t* rec = bcf_init1();
    while (bcf_itr_next(ifile, itervcf, rec) >= 0) {
      bcf_unpack(rec, BCF_UN_ALL);
      bcf_get_genotypes(hdr, rec, &gt, &ngt);
      if ((bcf_gt_allele(gt[sampleIndex*2]) != -1) && (bcf_gt_allele(gt[sampleIndex*2 + 1]) != -1) && (!bcf_gt_is_missing(gt[sampleIndex*2])) && (!bcf_gt_is_missing(gt[sampleIndex*2 + 1])) && (bcf_gt_is_phased(gt[sampleIndex*2 + 1]))) {
	int gt_type = bcf_gt_allele(gt[sampleIndex*2]) + bcf_gt_allele(gt[sampleIndex*2 + 1]);
	if (gt_type == 1) {
	  std::vector<std::string> alleles;
	  for(std::size_t i = 0; i<rec->n_allele; ++i) alleles.push_back(std::string(rec->d.allele[i]));
	  // Only bi-allelic variants
	  if (alleles.size() == 2) pV.push_back(TVariant(rec->pos, std::string(alleles[0]), std::string(alleles[1]), bcf_gt_allele(gt[sampleIndex*2])));
	}
      }
    }
    bcf_destroy(rec);
    hts_itr_destroy(itervcf);
  }
  if (gt != NULL) free(gt);
  
  // Close VCF
  bcf_hdr_destroy(hdr);
  hts_idx_destroy(bcfidx);
  bcf_close(ifile);

  return true;
}


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
  if (!_loadVariants(c, c.p1.string(), pV1)) return -1;
  if (!_loadVariants(c, c.p2.string(), pV2)) return -1;

  // Find common variants
  typedef std::vector<bool> THap;
  THap hap1;
  THap hap2;
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
	hap1.push_back(pV1[pos1].hap);
	hap2.push_back(pV2[pos2].hap);
      } else {
	++distinctSites1;
	++distinctSites2;
      }
      ++pos1;
      ++pos2;
    }
  }

  // Get differences
  THap diff1(hap1.size() - 1, 0);
  for(uint32_t i=1; i<hap1.size();++i)
    if (hap1[i-1]!=hap1[i]) diff1[i-1]=1;
  THap diff2(hap2.size() - 1, 0);
  for(uint32_t i=1; i<hap2.size();++i)
    if (hap2[i-1]!=hap2[i]) diff2[i-1]=1;

  // Switch errors
  int32_t switchcount = 0;
  for(uint32_t i=0; i<diff1.size(); ++i)
    if (diff1[i] != diff2[i]) ++switchcount;

  // Switch count
  std::cout << "Sample\tChromosome\tCommonSites\tDistinctBCF1\tDistinctBCF2\tSwitchErrors\tSwitchErrorRate" << std::endl;
  std::cout << c.sample << "\t" << c.chrom << "\t" << commonSites << "\t" << distinctSites1 << "\t" << distinctSites2 << "\t" << switchcount << "\t" << (double) switchcount / (double) diff1.size() << std::endl;

  
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

  return evaluatePhasing(c);
}
