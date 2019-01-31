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
## 4 Mapping ENSEMBL IDs to HGNC symbols
## 5 Network statistics
## 6 Biological validation: network properties
## 7 Annotation of the example gene set
## 8 References
## 9 Acknowledgements
