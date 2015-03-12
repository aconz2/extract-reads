## extract-reads

This is a tool to filter reads which likely appear in a set of reference contigs.
It processes reads in parallel and is quick.
It uses jellyfish internally to do kmer counting of the references.

## Usage

```
Extract reads which likely appear in reference:
  -h [ --help ]             Show help message
  -k [ --kmer-length ] arg  *Kmer length to use
  -s [ --hash-size ] arg    *Initial hash table size for Jellyfish (will grow 
                            if full)
  -c [ --cutoff ] arg       *A read with a k-mer which appears >= this is kept
  --single arg              Single fastq file
  --paired arg              Paired reads fastq files. give in order <left> 
                            <right>
  -r [ --reference ] arg    Reference file, fasta
  -t [ --threads ] arg (=1) Number of threads to use
```

## Install
```
make
make install
```

Tested with:
  - `gcc >= 4.8`
  - `boost >= 1.5.4`
  - `intel tbb >= 2015.0`
  - `jellyfish >= 2.2`
Requires `BOOST_ROOT` and `TBBROOT` to be set in order to compile. And 
`pkg-config` must be able to find the relevant jellyfish files.
You may change `PREFIX` to change the install location, defaulting to `/usr/local`

