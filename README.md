# `BCB420.2019.hu.MAP`



###### [Yin Yin](https://orcid.org/0000-0001-9168-488X), University of Toronto, Canada. &lt;yin.yin@mail.utoronto.ca&gt;
## 1 About this package:
This package describes the workflow to utilize human protein complex data from the hu.MAP. The workflow demonstrates how to check if HGNC symbols inside database is latest, how to annotate the example gene set, and also provides examples of computing database statistics.

The package serves dual duty, as an RStudio project, as well as an R package that can be installed. Package checks pass without errors, warnings, or notes.
#### In this project...

## 2 hu.MAP Data
hu.MAP refers to Human Protein Complex Map which attempts to address the lack of understanding of protein complex in human. This database combines serveral large scale protein interaction datasets to get most comprehensive understanding of human protein complexes. Protein inteaction network shows co-compelx potein pairs that are observed in protein complex map with corresponding probability score.The enrichments data and edge data are also included in the additional files of hu.MAP. 

This main document describes work with [Protein Interaction Network with probability scores (genenames)](http://proteincomplexes.org/download) [Drew, Kevin, et al. ](http://msb.embopress.org/content/13/6/932)

You may also need to download additional file: Enrichment table and Edge table.
#### 2.1 Data semantics
![alt text](http://msb.embopress.org/sites/default/files/highwire/msb/13/6/932/embed/graphic-1.gif)
<br /> This data integrated three protein interaction networks, BioPlex, Hein et al, and Wan et al into an combined protein complex network and clustered to identify protein complexes through the synthesis of over 9,000 published mass spectrometry experiment. Parameters for the SVM and clustering algorithms were optimized on a training set of literature‐curated complexes and validated on a test set of complexes.
## 3 Data download and cleanup
To download this package related hu.MAP data:
1. Navigate to the [hu.MAP database download page](http://proteincomplexes.org/download)
2. Click following files to download:
Protein Interaction Network with probability scores (genenames)
Enrichment table
Edge table
3. Place them in the same directory.
## 4 Mapping Entrez IDs to HGNC symbols
2. Click following files to download:
<br />Protein Complex Map
```txt
# Example
# Each line represent a complex and the listed genes() are composed 
153129 10670 64121
441502 3024
2648 26009 10474 57325 8850 55689 6871
3189 144983

```
<br />Protein Complex Map (genenames)
```txt
# Example
# Each line represent a complex and the listed genes(HGNC symbols) are composed 
SLC38A9	RRAGA	RRAGC
441502	HIST1H1A
KAT2A	ZZZ3	TADA3	KAT14	KAT2B	YEATS2	TADA2A
HNRNPH3	HNRNPA1L2

```
<br />Protein Interaction Network with probability scores
```txt
# Example
# Each line represent a co-protein (Entrez IDs) complex with the corresponding svm probability score 
996	64682	1.0
9861	5706	1.0
9861	5700	1.0
9774	9967	1.0
```
<br />Protein Interaction Network with probability scores (genenames)
```txt
# Example
# Each line represent a co-protein (HGNC symbols) complex with the corresponding svm probability score 
9YEATS4	VPS72	1.0
STON2	AP2M1	1.0
SRSF9	SRSF1	1.0
SNRPA1	SNRPD1	1.0
```
<br />Enrichment table
```txt
# Example
# Output from gprofiler for each complex, FDR-corrected hypergeometric p <= 0.05
complex_id	corr_pval	t_count	q_count	qandt_count	qandt_by_q	qandt_by_t	term_id	t_type	t_group	t_name	depth_in_group	qandt_list
0	1.60e-08	24	3	3	1.000	0.125	GO:0071230	BP	1	cellular response to amino acid stimulus	1	Q9HB90,Q7L523,Q8NBW4
0	1.82e-08	25	3	3	1.000	0.120	GO:0032008	BP	1	positive regulation of TOR signaling	1	Q9HB90,Q7L523,Q8NBW4
0	2.31e-08	27	3	3	1.000	0.111	GO:0043200	BP	1	response to amino acid	1	Q9HB90,Q7L523,Q8NBW4
```
<br />Edge table
```txt
# Example
# List of edges in the complex map with svm probability score and boolean values for each evidence type determining support for the edge
# Complex_id is made of complex number and gene ids
complex_id	corr_pval	t_count	q_count	qandt_count	qandt_by_q	qandt_by_t	term_id	t_type	t_group	t_name	depth_in_group	qandt_list
id1	score	fractions	bioplex	hein	bioplex_prey	hein_prey
0_153129 (pp) 0_10670	0.6968840000000001	False	True	False	True	False
0_153129 (pp) 0_64121	0.5299699999999999	False	False	False	True	False
0_10670 (pp) 0_64121	0.9627049999999999	False	True	False	True	False
```
3. Place them in the same directory.
## 4 Update HGNC symbols
hu.Map's data has two versions. One version uses Entrez ID and another version uses gene symbol. Although this database came out only at 2017 and the possibility that HGNC symbols update is small, the version of using gene symbol maybe out-dated in the future. Therefore, we provide a way to check and update HGNC symbols.

Preparation:
BioMart

### 4.1 Import all data into R
   ```R
   #for all data file
   tmp <- read.table(filepath,
                          sep  = " ",
                         fill=TRUE)
   tmp_genename <- read.table(filepath,
                          sep  ='\t',
                         fill=TRUE,stringsAsFactors = FALSE)
   ```

### 4.2 Update HGNC symbol 
&nbsp;

```R 
  # Fetch HGNC from github
  myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
  load(url(myURL))
  # Use bioMart to map Entrez to HGNC symbols
  
  myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

  dbIDs <- biomaRt::getBM(attributes=c('entrezgene','hgnc_symbol'),mart = ensembl)

  # Check if genenames in the file is the same as in dbIDS
 sej <- matrix( ncol = 3)
 for (i in 1:20){
           for (n in 1:ncol(normal_genename)){
                   if (!((normal_genename[i,n]) %in% HGNC$sym) & nchar(sub('\\.[0-9]+', '', normal_genename[i,n])) != 0){ sej <- rbind(sej,c(toString(normal_genename[i,n]),  i,n))}
                   
          }}
  #Delete first line that only contain NA
  sej <- sej[2:nrow(sej),]
  ```
  Now, we face serveral conditions for changing unsame genenames to latest genenames
  1) The unsame symbols are in either HGNC$prev or HGNC$synonym.
  We update in file
  2) The unsame symbols are gene IDs
   We update in file
  3) The unsame symbols are Ensmbol IDs
  ``` 
   
 for (i in 1:nrow(sej)) {
    #The unsame symbols are in either HGNC$prev or HGNC$synonym.
    iPrev <- grep(sej[i, 1], HGNC$prev)[1] # take No. 1 if there are several
    if (!is.na(iPrev)) {
        row <- as.numeric(sej[i, 2])
        col <- as.numeric(sej[i,3])
        tmp_genename[row,col] <- HGNC$sym[iPrev]
        sej[i,1] <- "changed"
    } else {
        iSynonym <- grep(sej[i, 1], HGNC$synonym)[1]
        if (!is.na(iSynonym)) {
            row <- as.numeric(sej[i,2])
            col <- as.numeric(sej[i,3])
            tmp_genename[row, col] <- HGNC$sym[iSynonym]
            sej[i,1] <- "changed"
        }else {
            #The unsame symbols are gene IDS
            iGeneids <- grep(sej[i, 1], HGNC$GeneID)[1]
            if (!is.na(iGeneids)) {
                row <- as.numeric(sej[i, 2])
                col <- as.numeric(sej[i, 3])
                tmp_genename[row, col] <- HGNC$sym[iGeneids]
                sej[i, 1] <- "changed"
            } else {
                #The unsame symbols are ensembl ID
                iEnsids <- grep(sej[i, 1], HGNC$EnsID)[1]
                if (!is.na(iEnsids)) {
                    row <- as.numeric(sej[i, 2])
                    col <- as.numeric(sej[i, 3])
                    tmp_genename[row, col] <- HGNC$sym[iEnsids]
                    sej[i, 1] <- "changed"
                }            
            }
        }}}
```
Now, in ```sej``` we have left ids that cannot mapping to other known symbols. So, we will change them into NA.
```
for (i in 1:nrow(sej)) {
         if (!identical(sej[i,1], "changed")){row <- as.numeric(sej[i, 2])
        col <- as.numeric(sej[i, 3]) 
    dup_geno[row, col] <- NA}}
```

### 4.3 Final validation
## 5 Annotation  gene set
First, we need to analyze our hu.MAP.

### 5.1 Complex statistics
We need to know how number of subunit of complex distribution

```
df <- matrix(0,nrow = ncol(dup_geno))

for (i in 1:nrow(dup_geno)){
          nSubunit <- 0
          n <- 1
          for (n in 1:ncol(dup_geno)){
                  if (is.na(dup_geno[i,n])){
                        nSubunit <- nSubunit + 1
             } else if (!nchar(dup_geno[i,n]) == 0){
                  nSubunit <- nSubunit + 1}
              df[i,1] <- nSubunit}
     }

hist(newdf,
     
     ylim=c(0,4000),col = "#3fafb388",
     main = "Number of subunit of complex distribution",
     xlab = "number of subunit in a complex ",
     ylab = "Counts")
```

### 5.2 Co-complex score statistics
We want to know the distribution of co-complex

### 5.3 Network statistics
we want
###
## 6 References
## 7 Acknowledgements
