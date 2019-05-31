
# About the data

# Datasets

## 1\. Krause/Niazi et al. DNA data

The datasets below contains data from two replicates: first replicate is
obtained using SQK\_LSK108, while the second one is obtained using
SQK-LSK109.

##### 1\. dna-krause-lsk108\_109-standard\_basecalling-tailfindr\_estimates-two\_replicates-with\_filepaths.csv

tailfindr estimations using data that has been basecalled with standard
model.

##### 2\. dna-krause-lsk108\_109-flipflop\_basecalling-tailfindr\_estimates-two\_replicates-with\_filepaths.csv

tailfindr estimations using data that has been basecalled with flipflop
model.

##### 3\. dna-krause-lsk108\_109-flipflop\_basecalling-decoded\_barcodes-two\_replicates-with\_filepaths.csv

Barcode assignment output using data that has been basecalled with
flipflop
model.

##### 4\. dna-krause-lsk108\_109-standard\_basecalling-transcript\_alignment\_start-two\_replicates-with\_filepaths.csv

Transcript alignment start information using data that has been
basecalled with standard
model.

##### 5\. dna-krause-lsk108\_109-flipflop\_basecalling-transcript\_alignment\_start-two\_replicates-with\_filepaths.csv

Transcript alignment start information using data that has been
basecalled with flipflop
model.

##### 6\. dna-krause-lsk108\_109-standard\_basecalling-moves\_in\_tail-two\_replicates-with\_filepaths.csv

Total moves within the poly(A)/(T) tail boundaries in the data that has
been basecalled with standard
model.

##### 7\. dna-krause-lsk108\_109-flipflop\_basecalling-moves\_in\_tail-two\_replicates-with\_filepaths.csv

Total moves within the poly(A)/(T) tail boundaries in the data that has
been basecalled with flipflop model.

## 2\. Krause/Niazi et al. RNA data

For RNA data, we obtained two replicates using SQK-RNA001 sequencing
kit, and a third replicate (after receiving the reviews) using
SQK-RNA002 kit. Replicates using SQK-RNA001 had the reverse
transcription step, whereas the SQK-RNA002 replicate was obtained by
omitting this
step.

##### rna-krause-rna001-standard\_basecalling-tailfindr\_estimates-two\_replicates-with\_filepaths.csv

tailfindr estimations using data that has been basecalled with standard
model. The data contains both of the SQK-RNA001
replicates.

##### rna-krause-rna001-standard\_basecalling-nanopolish\_estimates-two\_replicates-with\_filepaths.csv

Nanopolish estimations using data that has been basecalled with standard
model. The data contains both of the SQK-RNA001
replicates.

##### rna-krause-rna001-standard\_basecalling-decoded\_barcodes-two\_replicates-with\_filepaths.csv

Barcode assignment output using data that has been basecalled with
standard model. The data contains both of the SQK-RNA001
replicates.

##### rna-krause-rna001-standard\_basecalling-transcript\_alignment\_start-two\_replicates-with\_filepaths.csv

Transcript alignment start information using data that has been
basecalled with standard model. The data contains both of the SQK-RNA001
replicates.

##### rna-krause-rna002-standard\_basecalling-tailfindr\_estimates-with\_filepaths.csv

tailfindr estimations using data that has been basecalled with standard
model. The data contains a single replicate from
SQK-RNA002.

##### rna-krause-rna002-standard\_basecalling-nanopolish\_estimates-with\_filepaths.csv

Nanopolish estimations using data that has been basecalled with standard
model. The data contains a single replicate from
SQK-RNA002.

##### rna-krause-rna002-standard\_basecalling-decoded\_barcodes-with\_filepaths.csv

Barcode assignment output using data that has been basecalled with
standard model. The data contains a single replicate from
SQK-RNA002.

##### rna-krause-rna002-standard\_basecalling-transcript\_alignment\_start-with\_filepaths.csv

Transcript alignment start information using data that has been
basecalled with standard model. The data contains a single replicate
from
SQK-RNA002.

## 3\. Workman et al. DNA data

##### rna-workman-rna001-standard\_basecalling-tailfindr\_estimates-all\_datasets-with\_filepaths.csv

tailfindr estimations for Workman et al. data re-basecalled with
standard
model.

##### rna-workman-rna001-standard\_basecalling-nanopolish\_estimates-all\_datasets-with\_filepaths.csv

Nanopolish estimations for Workman et al. data re-basecalled with
standard model.

## 4\. ONT mer model files

##### r9.4\_180mv\_450bps\_6mer\_template\_median68pA.model

DNA model. This is used for calculation of thresholds used in our DNA
tail-finding algorithm.

##### r9.4\_180mv\_70bps\_5mer\_RNA\_template\_median69pA.model

RNA model. This is used for calculation of thresholds used in our RNA
tail-finding algorithm.

# Column descriptions

#### tailfindr CSV files

1.  read\_id: Read ID
2.  read\_type: Whether the read is poly(A)/poly(T) or invalid. Only
    reported for DNA datasets.
3.  tail\_is\_valid: Whether a poly(A)-tailed read is a full-length read
    or not. This is important because a poly(A) tail is at the end of
    the read, and premature termination of reads is prevelant in cDNA.
    Only reported for DNA datasets.
4.  tail\_start: Sample index of start site of the tail in raw data
5.  tail\_end: Sample index of end site of the tail in raw data
6.  samples\_per\_nt: Read rate in terms of samples per nucleotide
7.  tail\_length: Tail length in nucleotides. It is the difference
    between tail\_end and tail\_start divided by samples\_per\_nt
8.  file\_path: Full read path. Only relevant for internal use within
    Valen lab.
9.  replicate: Replicate number

#### Barcode assignment/decoding CSV files

1.  read\_id: Read ID
2.  file\_path: Full read path. Only relevant for internal use within
    Valen lab.
3.  read\_type: In case of DNA whether a read is GFP-containing Poly(A)
    or poly(T) read, or an invalid read. Incase of RNA, whether a read
    is GFP-containing Poly(A) read, or an invalid read.
4.  read\_length: Length of the read in terms of bases reported in that
    read  
5.  read\_too\_long: Is the read greater that 900 nt
6.  read\_too\_short: Is the read shorter that 900 nt
7.  nas\_gfp: Normalized alignment score for GFP alignment
8.  nas\_rc\_gfp: Normalized alignment score for reverse complement of
    GFP alignment  
9.  nas\_fp: Normalized alignment score for front primer alignment  
10. nas\_rc\_fp: Normalized alignment score for alignment of reverse
    complement of front primer
11. nas\_10bp: Normalized alignment score for alignment of barcode 10
12. nas\_30bp: Normalized alignment score for alignment of barcode 30  
13. nas\_40bp: Normalized alignment score for alignment of barcode 40
14. nas\_60bp: Normalized alignment score for alignment of barcode 60
15. nas\_100bp: Normalized alignment score for alignment of barcode 100
16. nas\_150bp: Normalized alignment score for alignment of barcode 150
17. nas\_slip: Normalized alignment score for alignment of slip barcode.
    Ignore it. For internal use only.
18. barcode: Barcode with the highest normalized alignment score
19. barcode\_tie: Is there any other barcode with same highest alignment
    score.
20. barcode\_2: Which barcode is it that has the same highest alignment
    score as the first barcode.  
21. barcode\_passed\_threshold: Does the normalized alignment score of
    the barcode with the highest alignment score pass the minimum
    threshold of 0.6.
22. replicate: Replicate number

#### Transcript alignment start CSV files

1.  transcript\_alignment\_start: Sample index of the junction point of
    the eGFP transcript and the poly(A)/(T) tail
2.  read\_id: Read ID
3.  file\_path: Full file path. Only relevant for internal use within
    Valen lab.

#### Moves in the tail CSV files

1.  read\_id: Read ID
2.  moves\_in\_tail\_st: Moves in the tail between tail start and end
    boundaries for data that has been basecalled with standard model  
3.  moves\_in\_tail\_ff: Moves in the tail between tail start and end
    boundaries for data that has been basecalled with flipflop model

#### Nanopolish estimation CSV files

1.  readname: Read ID
2.  contig
3.  position
4.  leader\_start
5.  adapter\_start  
6.  polya\_start: Sample index of the start of the tail
7.  transcript\_start: Sample index of the end of the tail
8.  read\_rate
9.  polya\_length
10. qc\_tag  
11. file\_path: Full read path. Only relevant for internal use within
    Valen lab.
