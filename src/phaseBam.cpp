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
  int32_t blockcounter;
  std::string sample;
  boost::filesystem::path h1bam;
  boost::filesystem::path h2bam;
  boost::filesystem::path bamfile;
  boost::filesystem::path vcffile;
};

struct Snp {
  uint32_t blockid;
  uint32_t pos;
  char h1;
  char h2;
  
  Snp() : blockid(0), pos(0), h1('N'), h2('N') {}
  Snp(uint32_t b, uint32_t p, char r, char a) : blockid(b), pos(p), h1(r), h2(a) {}
};

template<typename TRecord>
struct SortSnps : public std::binary_function<TRecord, TRecord, bool> {
  inline bool operator()(TRecord const& s1, TRecord const& s2) const {
    return s1.pos < s2.pos;
  }
};

inline uint32_t
lastAlignedPosition(bam1_t const* rec) {
  uint32_t* cigar = bam_get_cigar(rec);
  uint32_t alen = 0;
  for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CDEL)) alen += bam_cigar_oplen(cigar[i]);
  return rec->core.pos + alen;
}

template<typename TConfig, typename TSnpVector>
inline bool
_loadMarkers(TConfig& c, std::string const& chrName, int32_t const chrLength, TSnpVector& snps) {
  typedef typename TSnpVector::value_type TSnp;
  
  // Load bcf file
  htsFile* ifile = hts_open(c.vcffile.string().c_str(), "r");
  hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);
  int32_t sampleIndex = -1;
  for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
    if (hdr->samples[i] == c.sample) sampleIndex = i;
  if (sampleIndex < 0) return false;

  // Genotypes
  int ngt = 0;
  int32_t* gt = NULL;

  // Collect Snps for this chromosome
  int32_t chrid = bcf_hdr_name2id(hdr, chrName.c_str());
  if (chrid < 0) return false;
  hts_itr_t* itervcf = bcf_itr_queryi(bcfidx, chrid, 0, chrLength);
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
	  // Only bi-allelic SNPs
	  if ((alleles.size() == 2) && (alleles[0].size() == 1) && (alleles[1].size() == 1)) {
	    if (bcf_gt_allele(gt[sampleIndex*2]) == 1) snps.push_back(TSnp(c.blockcounter, rec->pos, alleles[1][0], alleles[0][0]));
	    else snps.push_back(TSnp(c.blockcounter, rec->pos, alleles[0][0], alleles[1][0]));
	  }
	}
      } else {
	// New phased block
	++c.blockcounter;
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
    
  // Sort Snps by position
  std::sort(snps.begin(), snps.end(), SortSnps<TSnp>());
  
  return true;
}


template<typename TConfig>
inline int
phaseBamRun(TConfig& c) {

#ifdef PROFILE
  ProfilerStart("phasebam.prof");
#endif

  // Load bam files
  samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
  hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
  bam_hdr_t* hdr = sam_hdr_read(samfile);

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Assign reads to haplotypes" << std::endl;
  boost::progress_display show_progress(hdr->n_targets);

  // Open output file
  samFile* h1bam = sam_open(c.h1bam.string().c_str(), "wb");
  sam_hdr_write(h1bam, hdr);

  samFile* h2bam = sam_open(c.h2bam.string().c_str(), "wb");
  sam_hdr_write(h2bam, hdr);

  // Assign reads to SNPs
  uint32_t assignedReads = 0;
  uint32_t unassignedReads = 0;
  uint64_t assignedBases = 0;
  uint64_t unassignedBases = 0;
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    ++show_progress;
    // New chromosome -> new phased block
    ++c.blockcounter;
    
    // Load variation data
    typedef std::vector<Snp> TSnpVector;
    TSnpVector snps;
    std::string chrName(hdr->target_name[refIndex]);
    if (_loadMarkers(c, chrName, hdr->target_len[refIndex], snps)) {
      // Assign reads to haplotypes
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	uint32_t hp1votes = 0;
	uint32_t hp2votes = 0;
	TSnpVector::const_iterator iSnp = std::lower_bound(snps.begin(), snps.end(), Snp(0, rec->core.pos, 'A', 'A'), SortSnps<Snp>());
	TSnpVector::const_iterator iSnpEnd = std::upper_bound(snps.begin(), snps.end(), Snp(0, lastAlignedPosition(rec), 'A', 'A'), SortSnps<Snp>());
	if (iSnp != iSnpEnd) {
	  // Get read sequence
	  std::string sequence;
	  sequence.resize(rec->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(rec);
	  for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	  
	  // Parse CIGAR
	  uint32_t* cigar = bam_get_cigar(rec);
	  for(;iSnp != iSnpEnd; ++iSnp) {
	    uint32_t gp = rec->core.pos; // Genomic position
	    uint32_t sp = 0; // Sequence position
	    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	      if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
		if (gp + bam_cigar_oplen(cigar[i]) < iSnp->pos) {
		  gp += bam_cigar_oplen(cigar[i]);
		  sp += bam_cigar_oplen(cigar[i]);
		} else {
		  for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
		    if (gp == iSnp->pos) {
		      if (sequence[sp] == iSnp->h1) ++hp1votes;
		      else if (sequence[sp] == iSnp->h2) ++hp2votes;
		    }
		  }
		  // SNP has been found -> break
		  break;
		}
	      }
	    }
	  }
	}
	int32_t hp = 0;
	if (hp1votes > 2*hp2votes) hp = 1;
	else if (hp2votes > 2*hp1votes) hp = 2;
	if (hp) {
	  ++assignedReads;
	  assignedBases += rec->core.l_qseq;
	  int32_t ps = c.blockcounter;
	  bam_aux_append(rec, "PS", 'i', 4, (uint8_t*)&ps);
	  if (hp == 1) {
	    bam_aux_append(rec, "HP", 'i', 4, (uint8_t*)&hp);
	    sam_write1(h1bam, hdr, rec);
	  } else {
	    bam_aux_append(rec, "HP", 'i', 4, (uint8_t*)&hp);
	    sam_write1(h2bam, hdr, rec);
	  }
	} else {
	  ++unassignedReads;
	  unassignedBases += rec->core.l_qseq;
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
  }

  // Close bam
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);

  // Close output BAMs
  sam_close(h1bam);
  sam_close(h2bam);

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

  // Statistics
  uint64_t sumReads = assignedReads + unassignedReads;
  uint64_t sumBases = assignedBases + unassignedBases;
  std::cout << "Assigned Reads=" << assignedReads << ", Unassigned Reads=" << unassignedReads << ", Fraction assigned=" << (float) assignedReads / (float) sumReads << std::endl;
  std::cout << "Assigned Bases=" << assignedBases << ", Unassigned Bases=" << unassignedBases << ", Fraction assigned=" << (float) assignedBases / (float) sumBases << std::endl;


#ifdef PROFILE
  ProfilerStop();
#endif


  return 0;
}

int main(int argc, char **argv) {
  Config c;
  c.blockcounter = 0;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("hap1,p", boost::program_options::value<boost::filesystem::path>(&c.h1bam)->default_value("h1.bam"), "haplotype 1 BAM file")
    ("hap2,q", boost::program_options::value<boost::filesystem::path>(&c.h2bam)->default_value("h2.bam"), "haplotype 2 BAM file")
    ("sample,s", boost::program_options::value<std::string>(&c.sample)->default_value("NA12878"), "sample name")
    ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "phased SNP BCF file")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamfile), "input bam file")
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
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("vcffile"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -s NA12878 -v <snps.bcf> --hap1 <h1.bam> --hap2 <h2.bam> <unphased.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  } 

  // Check input BAM file
  if (vm.count("input-file")) {
    if (!(boost::filesystem::exists(c.bamfile) && boost::filesystem::is_regular_file(c.bamfile) && boost::filesystem::file_size(c.bamfile))) {
      std::cerr << "Input BAM file is missing: " << c.bamfile.string() << std::endl;
      return 1;
    }
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    if (samfile == NULL) {
      std::cerr << "Fail to open file " << c.bamfile.string() << std::endl;
      return 1;
    }
    hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
    if (idx == NULL) {
      std::cerr << "Fail to open index for " << c.bamfile.string() << std::endl;
      return 1;
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    if (hdr == NULL) {
      std::cerr << "Fail to open header for " << c.bamfile.string() << std::endl;
      return 1;
    }
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }
  
  // Check VCF/BCF file
  if (vm.count("vcffile")) {
    if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
      std::cerr << "Input SNP VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
      return 1;
    }
    htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
    if (ifile == NULL) {
      std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
      return 1;
    }
    hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
    if (bcfidx == NULL) {
      std::cerr << "Fail to open index file for " << c.vcffile.string() << std::endl;
      return 1;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    if (hdr == NULL) {
      std::cerr << "Fail to open header for " << c.vcffile.string() << std::endl;
      return 1;
    }
    bcf_hdr_destroy(hdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ifile);
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  return phaseBamRun(c);
}
