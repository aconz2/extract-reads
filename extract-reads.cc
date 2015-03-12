#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#include <string>
#include <fstream>
#include <iostream>
#include <iterator>
#include <utility>
#include <algorithm>
#include <functional>
#include <memory>

#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/large_hash_iterator.hpp>

#include "tbb/pipeline.h"
#include "tbbff.h"

using jellyfish::mer_dna;
typedef std::vector<std::string>::iterator string_iter;
typedef jellyfish::stream_manager<string_iter> stream_manager;
typedef jellyfish::mer_overlap_sequence_parser<stream_manager> sequence_parser;
typedef jellyfish::mer_iterator<sequence_parser, mer_dna> mer_iterator;
typedef jellyfish::whole_sequence_parser<std::string> read_parser;

/* ========== mer_counter ========== */
class mer_counter : public jellyfish::thread_exec {
  bool            canonical_;
  int             nb_threads_;
  mer_hash&       ary_;
  stream_manager  streams_;
  sequence_parser parser_;

public:
  mer_counter(int nb_threads, mer_hash& ary, string_iter file_begin, string_iter file_end, bool canonical) :
    canonical_(canonical),
    nb_threads_(nb_threads),
    ary_(ary),
    streams_(file_begin, file_end),
    parser_(mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
  { }

  virtual void start(int thid) {
    for(mer_iterator mers(parser_, canonical_) ; mers; ++mers) {
      ary_.add(*mers, 1);
    }
    ary_.done();
  }

};

/* ========== fastq struct ========== */
struct fastq {
  std::string id;
  std::string sequence;
  std::string plus;
  std::string quality;

  friend std::istream & operator>>(std::istream & is, fastq & record);
  friend std::ostream & operator<<(std::ostream & os, fastq & record);
  friend std::ostream & operator<<(std::ostream & os, const fastq & record);
};

std::istream & operator>>(std::istream & is, fastq & record) {
  std::getline(is, record.id);
  std::getline(is, record.sequence);
  std::getline(is, record.plus);
  std::getline(is, record.quality);
  return is;
};

std::ostream & operator<<(std::ostream & os, fastq & record) {
  return os << record.id << '\n'
            << record.sequence << '\n'
            << record.plus << '\n'
            << record.quality << '\n';
};

std::ostream & operator<<(std::ostream & os, const fastq & record) {
  return os << record.id << '\n'
            << record.sequence << '\n'
            << record.plus << '\n'
            << record.quality << '\n';
};

bool is_good_read(std::string &seq, mer_array *ary, unsigned int kmer_length, unsigned int cutoff) {
  if(seq.size() < kmer_length) {
    return false;
  }
  unsigned int count = 0;
  uint64_t val;
  // should eventually make an iterator for this!
  mer_dna mer;
  mer = seq.substr(0, kmer_length);
  size_t i = kmer_length;
  do {
    if(ary->get_val_for_key(mer.get_canonical(), &val)) {
      count++;
    }
    if(count >= cutoff) {
      return true;
    }
    mer.shift_left(seq[i++]);
  } while(i < seq.size());
  return false;
}

/* ========== main ========== */
namespace po = boost::program_options;

int main(int argc, char *argv[]) {

  const int max_reprobe = 126;
  const bool canonical = true;

  /* parallel parameters */
  // Experiments show these are pretty good numbers
  int buffer_size = 500;
  int token_multiplier = 4;

  /* default values */
  int counter_length = 7;
  int num_threads = 1;
  std::string out_dir = "./";

  /* manadatory arguments */
  unsigned int kmer_length;
  unsigned int cutoff;
  unsigned long hash_size;
  std::string single_reads;
  std::vector<std::string> paired_reads;
  std::string genes;

  po::options_description desc("Extract reads which likely appear in reference. Give one of --single or --paired");
  desc.add_options()
    ("help,h", "Show help message")
    ("kmer-length,k", po::value<unsigned int>(&kmer_length)->required(), "Kmer length to use")
    ("hash-size,s", po::value<unsigned long>(&hash_size)->required(), "Initial hash table size for Jellyfish (will grow if full)")
    ("cutoff,c", po::value<unsigned int>(&cutoff)->required(), "A read with a k-mer which appears >= this is kept")
    ("single", po::value<std::string>(&single_reads), "Single fastq file")
    ("paired", po::value<std::vector<std::string> >(&paired_reads)->multitoken(), "Paired reads fastq files. give in order <left> <right>")
    ("reference,r", po::value<std::string>(&genes)->required(), "Reference file, fasta")
    ("out,o", po::value<std::string>(&out_dir), "Output directory, MUST exist already")
    //("buffer,b", po::value<int>(&buffer_size)->default_value(buffer_size), "Internal! buffer size")
    //("multiplier,m", po::value<int>(&token_multiplier)->default_value(token_multiplier), "Internal! token multiplier")
    ("threads,t", po::value<int>(&num_threads)->default_value(num_threads), "Number of threads to use");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if(vm.count("help")) {
      std::cout << desc;
      exit(0);
    }
    po::notify(vm);
  } catch(po::error &e) {
    std::cerr << "ERROR:" << e.what() << std::endl << desc;
    exit(1);
  }

  if(out_dir[out_dir.length() - 1] != '/') {
    out_dir += "/";
  }

  mer_dna::k(kmer_length);
  mer_hash hash(hash_size, mer_dna::k() * 2, counter_length, num_threads, max_reprobe);
  /* ---------- Fill Jellyfish Hashtable ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Filling jellyfish hash ===" << std::endl;
    std::vector<std::string> in_files;
    in_files.push_back(genes); 

    std::vector<std::string>::iterator it = in_files.begin();

    mer_counter counter(num_threads, hash, in_files.begin(), in_files.end(), canonical);
    counter.exec_join(num_threads);
  }

  mer_array* ary = hash.ary();

  bool paired_end = paired_reads.size() == 2;

  /* ---------- Count kmers in genes ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Extracting reads ===" << std::endl;

    if(paired_end) {
      std::ifstream ifs_left;
      std::ifstream ifs_right;
      ifs_left.open(paired_reads[0]);
      ifs_right.open(paired_reads[1]);
      if(!ifs_left.good() || !ifs_right.good()) {
        std::cerr << "Problem opening paired reads file" << std::endl;
        exit(1);
      }
      std::string out_file_left = out_dir + paired_reads[0] + ".extracted";
      std::string out_file_right = out_dir + paired_reads[1] + ".extracted";
      std::ofstream ofs_left(out_file_left);
      std::ofstream ofs_right(out_file_right);
      if(!ofs_left.good() || !ofs_right.good()) {
        std::cerr << "Problem opening paired reads output file" << std::endl;
        exit(1);
      }
      std::istream_iterator<fastq> EOS;
      std::istream_iterator<fastq> iter_left(ifs_left);
      std::istream_iterator<fastq> iter_right(ifs_right);

      // this is pretty nasty, but kinda necessary to play nice with templates
      typedef boost::tuple<fastq, fastq> fastq_pair;
      auto pair_iterator_begin = boost::make_zip_iterator(boost::make_tuple(iter_left, iter_right));
      auto pair_iterator_end = boost::make_zip_iterator(boost::make_tuple(EOS, EOS));
      typedef decltype(pair_iterator_begin) fastq_pair_iterator;

      tbb::parallel_pipeline(num_threads * token_multiplier,
        tbbff::input_buffer<fastq_pair, fastq_pair_iterator>(pair_iterator_begin, pair_iterator_end, buffer_size) &
        tbbff::keep_if<fastq_pair>([ary, kmer_length, cutoff](fastq_pair &pair) -> bool {
          return is_good_read(pair.get<0>().sequence, ary, kmer_length, cutoff) ||
                 is_good_read(pair.get<1>().sequence, ary, kmer_length, cutoff);
        }) &
        tbb::make_filter<std::vector<fastq_pair>*, void>(
          tbb::filter::serial_out_of_order,
          [&ofs_left, &ofs_right](std::vector<fastq_pair> *buffer) -> void {
            for(fastq_pair pair: *buffer) {
              ofs_left << pair.get<0>();
              ofs_right << pair.get<1>();
            }
            buffer->clear();
            delete buffer;
          }
        )
      );
    } else {
      std::ifstream ifs;
      ifs.open(single_reads);
      if(!ifs.good()) {
        std::cerr << "Problem opening single reads file" << std::endl;
        exit(1);
      }
      std::string out_file = out_dir + single_reads + ".extracted";
      std::ofstream ofs(out_file);
      if(!ofs.good()) {
        std::cerr << "Problem opening single reads output file" << std::endl;
        exit(1);
      }
      std::istream_iterator<fastq> EOS;
      std::istream_iterator<fastq> ifs_iter(ifs);
      std::ostream_iterator<fastq> ofs_iter(ofs);

      tbb::parallel_pipeline(num_threads * token_multiplier,
        tbbff::input_buffer<fastq, std::istream_iterator<fastq> >(ifs_iter, EOS, buffer_size) &
        tbbff::keep_if<fastq>([ary, kmer_length, cutoff](fastq &read) -> bool {
          return is_good_read(read.sequence, ary, kmer_length, cutoff);
        }) &
        tbbff::output_buffer<fastq, std::ostream_iterator<fastq> >(ofs_iter)
      );
    } // else
  }// block scope

  return 0;
}
