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
#include <htslib/faidx.h>
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
  bool attemptPhasing;
  unsigned short minMapQual;
  int32_t mincov;
  std::string sample;
  boost::filesystem::path genome;
  boost::filesystem::path bamfile;
  boost::filesystem::path vcffile;
  boost::filesystem::path outfile;
};

template<typename TConfig>
inline int
annotateRefAlt(TConfig& c) {

#ifdef PROFILE
  ProfilerStart("phasebam.prof");
#endif
  // Load bam files
  samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
  hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
  bam_hdr_t* hdr = sam_hdr_read(samfile);

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Annotate REF & ALT" << std::endl;
  boost::progress_display show_progress(hdr->n_targets);

  // Output BCF
  htsFile *fp = hts_open(c.outfile.string().c_str(), "wb");
  bcf_hdr_t *hdr_out = bcf_hdr_init("w");
  now = boost::posix_time::second_clock::local_time();
  boost::gregorian::date today = now.date();
  std::string datestr("##fileDate=");
  datestr += boost::gregorian::to_iso_string(today);
  bcf_hdr_append(hdr_out, datestr.c_str());
  bcf_hdr_append(hdr_out, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">");
  bcf_hdr_append(hdr_out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  bcf_hdr_append(hdr_out, "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# REF support\">");
  bcf_hdr_append(hdr_out, "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# ALT support\">");
  std::string refloc("##reference=");
  refloc += c.genome.string();
  bcf_hdr_append(hdr_out, refloc.c_str());
  for (int i = 0; i<hdr->n_targets; ++i) {
    std::string refname("##contig=<ID=");
    refname += std::string(hdr->target_name[i]) + ",length=" + boost::lexical_cast<std::string>(hdr->target_len[i]) + ">";
    bcf_hdr_append(hdr_out, refname.c_str());
  }
  bcf_hdr_add_sample(hdr_out, c.sample.c_str());
  bcf_hdr_add_sample(hdr_out, NULL);
  bcf_hdr_write(fp, hdr_out);
  
  // Iterate chromosomes
  faidx_t* fai = fai_load(c.genome.string().c_str());
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    std::string chrName(hdr->target_name[refIndex]);
    ++show_progress;

    // Load het. markers
    typedef std::vector<Variant> TVariants;
    TVariants var;
    if (!_loadVariants(c.sample, chrName, c.vcffile.string(), false, 0, c.attemptPhasing, var)) continue;
    if (var.empty()) continue;

    // Load reference
    int32_t seqlen = -1;
    char* seq = NULL;
    seq = faidx_fetch_seq(fai, chrName.c_str(), 0, hdr->target_len[refIndex], &seqlen);

    // Iter bam file
    typedef std::vector<int32_t> TAlleleSupport;
    TAlleleSupport ref(var.size(), 0);
    TAlleleSupport alt(var.size(), 0);
    hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
    bam1_t* rec = bam_init1();
    while (sam_itr_next(samfile, iter, rec) >= 0) {
      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
      if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
      typename TVariants::const_iterator vIt = std::lower_bound(var.begin(), var.end(), Variant(rec->core.pos), SortVariants<Variant>());
      typename TVariants::const_iterator vItEnd = std::upper_bound(var.begin(), var.end(), Variant(lastAlignedPosition(rec)), SortVariants<Variant>());
      if (vIt != vItEnd) {
	// Get read sequence
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	
	// Parse CIGAR
	uint32_t* cigar = bam_get_cigar(rec);
	for(;vIt != vItEnd; ++vIt) {
	  int32_t gp = rec->core.pos; // Genomic position
	  int32_t sp = 0; // Sequence position
	  bool varFound = false;
	  for (std::size_t i = 0; ((i < rec->core.n_cigar) && (!varFound)); ++i) {
	    if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	      if (gp + (int32_t) bam_cigar_oplen(cigar[i]) < vIt->pos) {
		gp += bam_cigar_oplen(cigar[i]);
		sp += bam_cigar_oplen(cigar[i]);
	      } else {
		for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
		  if (gp == vIt->pos) {
		    varFound = true;
		    // Check REF allele
		    if (vIt->ref == std::string(seq + gp, seq + gp + vIt->ref.size())) {
		      // Check ALT allele
		      if ((sp + vIt->alt.size() < sequence.size()) && (sp + vIt->ref.size() < sequence.size())) {
			if (vIt->ref.size() == vIt->alt.size()) {
			  // SNP
			  if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) ++alt[vIt-var.begin()];
			  else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) ++ref[vIt-var.begin()];
			} else if (vIt->ref.size() < vIt->alt.size()) {
			  // Insertion
			  int32_t diff = vIt->alt.size() - vIt->ref.size();
			  std::string refProbe = vIt->ref + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			  if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->alt.size()) != refProbe)) ++alt[vIt-var.begin()];
			  else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->alt.size()) == refProbe)) ++ref[vIt-var.begin()];
			} else {
			  // Deletion
			  int32_t diff = vIt->ref.size() - vIt->alt.size();
			  std::string altProbe = vIt->alt + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			  if ((sequence.substr(sp, vIt->ref.size()) == altProbe) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) ++alt[vIt-var.begin()];
			  else if ((sequence.substr(sp, vIt->ref.size()) != altProbe) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) ++ref[vIt-var.begin()];
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    bam_destroy1(rec);
    hts_itr_destroy(iter);
    if (seqlen) free(seq);

    int32_t *gts = (int*) malloc(2 * sizeof(int));
    int32_t *rrcount = (int*) malloc(sizeof(int));
    int32_t *rvcount = (int*) malloc(sizeof(int));

    bcf1_t *rec_out = bcf_init();
    for(typename TVariants::const_iterator ip = var.begin(); ip != var.end(); ++ip) {
      rec_out->rid = bcf_hdr_name2id(hdr_out, hdr->target_name[refIndex]);
      rec_out->pos = ip->pos;
      std::string id(".");
      bcf_update_id(hdr_out, rec_out, id.c_str());
      std::string alleles;
      alleles += ip->ref + "," + ip->alt;
      bcf_update_alleles_str(hdr_out, rec_out, alleles.c_str());
      int32_t tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "PASS");
      bcf_update_filter(hdr_out, rec_out, &tmpi, 1);
      if (ip->ac >= 0) {
	tmpi = ip->ac;
	bcf_update_info_int32(hdr_out, rec_out, "AC", &tmpi, 1);
      }
      rrcount[0] = ref[ip - var.begin()];
      bcf_update_format_int32(hdr_out, rec_out, "RR", rrcount, 1);
      rvcount[0] = alt[ip - var.begin()];
      bcf_update_format_int32(hdr_out, rec_out, "RV", rvcount, 1);
      bool phased = false;
      if ((c.attemptPhasing) && (rrcount[0] + rvcount[0] >= c.mincov)) {
	if (rrcount[0] > 2 * rvcount[0]) {
	  phased = true;
	  gts[0] = bcf_gt_phased(0);
	  gts[1] = bcf_gt_phased(1);
	} else if (rvcount[0] > 2 * rrcount[0]) {
	  phased = true;
	  gts[0] = bcf_gt_phased(1);
	  gts[1] = bcf_gt_phased(0);
	}
      }
      if (!phased) {
	gts[0] = bcf_gt_missing;
	gts[1] = bcf_gt_missing;
      }
      bcf_update_genotypes(hdr_out, rec_out, gts, 2);
      bcf_write1(fp, hdr_out, rec_out);
      bcf_clear1(rec_out);
    }
    
    // Clean-up
    free(gts);
    free(rrcount);
    free(rvcount);
    bcf_destroy1(rec_out);
  }
  fai_destroy(fai);
  
  // Close BCF file
  bcf_hdr_destroy(hdr_out);
  hts_close(fp);

  // Index
  bcf_index_build(c.outfile.string().c_str(), 14);
  
  // Close bam
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);
  
  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

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
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "reference fasta file")
    ("mincov,m", boost::program_options::value<int32_t>(&c.mincov)->default_value(10), "min. coverage for phasing")
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. mapping quality")
    ("sample,s", boost::program_options::value<std::string>(&c.sample)->default_value("NA12878"), "sample name")
    ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "phased BCF file")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("phased.bcf"), "phased BCF output file")
    ("phase,p", "emit phased GTs for LOH chromosomes")
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
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("vcffile"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -g <ref.fa> -s NA12878 -v <snps.bcf> <input.bam>" << std::endl;
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

  // Check reference
  if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
    std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
    return 1;
  } else {
    faidx_t* fai = fai_load(c.genome.string().c_str());
    if (fai == NULL) {
      if (fai_build(c.genome.string().c_str()) == -1) {
	std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	return 1;
      } else fai = fai_load(c.genome.string().c_str());
    }
    fai_destroy(fai);
  }

  // Attempt Phasing
  if (vm.count("phase")) c.attemptPhasing = true;
  else c.attemptPhasing = false;

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  return annotateRefAlt(c);
}
