# Team Project, PUBH 6885: Computational Biology

### Project info

-   Team members: Chiraag, Varun
-   Link to paper: [Cell](https://www.cell.com/cell/fulltext/S0092-8674(21)00882-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421008825%3Fshowall%3Dtrue)
-   Link to data: [Broad Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1289/impaired-local-intrinsic-immunity-to-sars-cov-2-infection-in-severe-covid-19#/)
-   Our project involves the usage of various tools from `omicsEye` to analyze data from a scRNA-seq study analyzing differences in gene expression across Covid-19 severity. 

### Reproducing our results

- Scripts used for our final report can be found in `R/`.
  - `loading-data.R` performs data scraping and centralizing of large data files to a shared google drive. This contains functions that load in large amounts of data (that are now already processed and saved), and can be skipped when reproducing results
  - `aims.pdf` describes usage of the other scripts
  
