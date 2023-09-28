# TCGA methylation curation

This repository contains script to post-process methylation matrices originating from the GDC, specifically the TCGA project. The expected input is a probe-level matrix with methylation "beta" values and the output is a gene-level matrix aggregated according to user-defined parameters.

As written, we use a probe-mapping file for the Illumina Infinium 450k platform, obtained from https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv 


## Details

Aggregating probe-level results to the gene level can be a bit murky/ambiguous, and this section provides details which are reflected in the code. Of course, the code itself is also commented to highlight these features/decisions.

Probe-level annotation files typically have multiple mappings. As an example, consider the following row for probe ID "cg27416437":

```
IlmnID                                                                cg27416437
Name                                                                  cg27416437
AddressA_ID                                                             17807441
AlleleA_ProbeSeq               ACCRAAAAAACTCTCAACACRAAACTAAACTATTTAAAATAATTCC...
AddressB_ID                                                                  NaN
AlleleB_ProbeSeq                                                             NaN
Infinium_Design_Type                                                          II
Next_Base                                                                    NaN
Color_Channel                                                                NaN
Forward_Sequence               CTCCCGCCGATCTCCGTCTCCCTGCAGGGCCGACTCTTCAGCGACC...
Genome_Build                                                                37.0
CHR                                                                         22.0
MAPINFO                                                               29664006.0
SourceSeq                      CCGGAGAGGCTCTCAGCACGAGACTAGGCTGTTTGGAATGGTTCCG...
Chromosome_36                                                                 22
Coordinate_36                                                           27994006
Strand                                                                         F
Probe_SNPs                                                                   NaN
Probe_SNPs_10                                                                NaN
Random_Loci                                                                  NaN
Methyl27_Loci                                                               True
UCSC_RefGene_Name              EWSR1;EWSR1;EWSR1;EWSR1;RHBDD3;EWSR1;EWSR1;EWS...
UCSC_RefGene_Accession         NM_001163286;NM_013986;NM_001163286;NM_0011632...
UCSC_RefGene_Group             5'UTR;5'UTR;1stExon;1stExon;TSS200;1stExon;5'U...
UCSC_CpG_Islands_Name                                    chr22:29663622-29664903
Relation_to_UCSC_CpG_Island                                               Island
Phantom                                                                      NaN
DMR                                                                          NaN
Enhancer                                                                     NaN
HMM_Island                                                  22:27993619-27995436
Regulatory_Feature_Name                                     22:29663777-29664298
Regulatory_Feature_Group                                     Promoter_Associated
DHS                                                                          NaN
```

Looking at the UCSC build (37), we see that EWSR1 has coordinates of chr22:29,663,998 - 29,696,515 (+ strand) and RHBDD3 has coordinates chr22:29,655,844 - 29,663,914 (- strand). Thus, the probe location (anchored at chr22:29664006) is just barely "inside" EWSR1 and just upstream of RHBDD3. If we parse out the genes/groups (UCSC_RefGene_Name/UCSC_RefGene_Group) from the semicolon-delimited strings, we get (removing duplicates):

```
RHBDD3: TSS200
EWSR1: 1stExon
EWSR1: 5'UTR
```
Given the probe length of 128bp, it apparently covers the 5'UTR *and* the first exon.

Thus, when creating gene-level summaries, we might need to include the "beta" values for both genes, depending on the user-defined filters. For instance, if the user requests "TSS200", then we would only report RHBDD3. If they ask for "TSS200" and "5'UTR", then we include both genes.

In total, we have the following annotations:
```
TSS1500
TSS200
1stExon
Body
5'UTR
3'UTR
NaN
```
There is also a boolean(ish) field named "Enhancer" (values are "True" and `nan`), which we can also include. This would further limit those features to ones required to be enhancers.

To save time when preparing user files, we have a pre-processing script for the "raw" file we downloaded from Illumina. This takes the records as formatted above and creates the following (here for only `cg00018261` and `cg27416437`):
```
probe          gene       feature        enhancer
cg27416437     EWSR1      1stExon        NaN
cg27416437     EWSR1      5'UTR          NaN
cg27416437     RHBDD3     TSS200         NaN
cg00018261     ZCCHC12    TSS200         True
```

Thus, if a user requests "TSS200", we would perform an inner join between the full probe-level matrix and the following two rows:
```
probe          gene       feature        enhancer
cg27416437     RHBDD3     TSS200         NaN
cg00018261     ZCCHC12    TSS200         True
```
If "enhancer" flag is selected, we would also drop probe `cg27416437`.

If a user selects the union of 5'UTR and 1stExon, we would get
```
probe          gene       feature        enhancer
cg27416437     EWSR1      1stExon        NaN
cg27416437     EWSR1      5'UTR          NaN
```
However, we then need to recognize that we don't want to double count, so we need to only retain unique combinations of the `(probe, gene)` tuple.

In the end, we end up with a joined matrix something like (e.g. user selected union of TSS200, 5'UTR, and 1stExon):
```
probe          gene       sample1    sample2     ...    sampleN
cg27416437     EWSR1      0.91       0.76        ...      0.19
cg27416437     RHBDD3     0.91       0.76        ...      0.19
cg00018261     ZCCHC12    0.22       0.98        ...      0.79
```
Thus, we saw the repeat of the (cg27416437, EWSR1) key since the cg27416437 probe covered both 5'UTR *and* 1stExon, dropped it, and joined with the raw data matrix. 

As a more general case, consider the situation where we have multiple probes assigned to the gene. For instance, consider cg1 and cg3 which are both assigned to EWSR1; cg1 covers the first exon and cg3 covers the body.
```
probe          gene       sample1    sample2     ...    sampleN
cg1            EWSR1      0.91       0.76        ...      0.19
cg1            RHBDD3     0.91       0.76        ...      0.19
cg2            ZCCHC12    0.22       0.98        ...      0.79
cg3            EWSR1      0.71       0.26        ...      0.08
```

**Important notes/caveats:**

- We drop probes that are not annotated to any gene. For instance, some probes do not have a corresponding gene and look like:
```
IlmnID                                                                cg00213748
Name                                                                  cg00213748
AddressA_ID                                                             30703409
AlleleA_ProbeSeq               TTTTAACACCTAACACCATTTTAACAATAAAAATTCTACAAAAAAA...
AddressB_ID                                                           36767301.0
AlleleB_ProbeSeq               TTTTAACGCCTAACACCGTTTTAACGATAAAAATTCTACAAAAAAA...
Infinium_Design_Type                                                           I
Next_Base                                                                      A
Color_Channel                                                                Red
Forward_Sequence               TCTGTGGGACCATTTTAACGCCTGGCACCGTTTTAACGATGGAGGT...
Genome_Build                                                                37.0
CHR                                                                            Y
MAPINFO                                                                8148233.0
SourceSeq                      CGCCCCCTCCTGCAGAACCTCCATCGTTAAAACGGTGCCAGGCGTT...
Chromosome_36                                                                  Y
Coordinate_36                                                            8208233
Strand                                                                         R
Probe_SNPs                                                                   NaN
Probe_SNPs_10                                                                NaN
Random_Loci                                                                  NaN
Methyl27_Loci                                                                NaN
UCSC_RefGene_Name                                                            NaN
UCSC_RefGene_Accession                                                       NaN
UCSC_RefGene_Group                                                           NaN
UCSC_CpG_Islands_Name                                       chrY:8147877-8148210
Relation_to_UCSC_CpG_Island                                              S_Shore
Phantom                                                                      NaN
DMR                                                                          NaN
Enhancer                                                                     NaN
HMM_Island                                                     Y:8207555-8208234
Regulatory_Feature_Name                                                      NaN
Regulatory_Feature_Group                                                     NaN
DHS                                                                          NaN
```
For our purposes, we drop these rows from further consideration.  


## Aggregation

Given the example matrix immediately above, we have a choice to be made when considering how we combine multiple probe values targeting the same gene. Recall the goal of this module is to aggregate methylation "beta" values at the gene level for use with a tool like NetZoo's DRAGON (https://netzoo.github.io/zooanimals/dragon/).

The options can include
- sum (as implemented in the DRAGON publication)
- arithmetic mean
- geometric mean

Therefore, the resulting matrix would be a proxy for the "methylation load" of a gene. The "sum" aggregation possibly biases genes with high probe coverage (or longer genes), but we include it since this is how the data was prepared for DRAGON.