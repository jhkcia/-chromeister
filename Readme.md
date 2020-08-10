For running the code you can use:
python Main.py -db [db-sequece-path] -query [query-sequence-path]
remember to create output dir before running the code.

for example:
python Main.py -db /home/jalal/Uni/term2/bioAlg/project/Prij/canfam3_1_0009.fsa -query /home/jalal/Uni/term2/bioAlg/project/Prij/sequence.fasta

extra-args can be:

CHROMEISTER ALGORITHME

optional arguments:
  -h, --help            show this help message and exit
  -db DB                database sequence fasta file path
  -query QUERY          database sequence fasta file path
  -k K                  k-mer length
  -hash-key-len HASH_KEY_LEN
                        hash function key length
  -z Z                  hash function Z parameter
  -dim DIM              hit matrix dimension
  -o O                  output directory
  -diag-len DIAG_LEN    diagonal length
  -event-sampling-size EVENT_SAMPLING_SIZE
                        event sampling size
  -reward REWARD        reward for growing regions
  -penalty PENALTY      penalty for growing regions
  -sidePenalty SIDEPENALTY
                        sidePenalty for growing regions
  -max-hsps MAX_HSPS    max HSPS size for growing regions
  -gr-th GR_TH          threshold for growing regions
  -w-size W_SIZE        window size for growing regions