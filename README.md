# Team Project, PUBH 6885: Computational Biology

### Project Information

-   Team members: Chiraag Gohel, Varun Subramaniam
-   Link to paper: [Cell](https://www.cell.com/cell/fulltext/S0092-8674(21)00882-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421008825%3Fshowall%3Dtrue)
-   Link to data: [Broad Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1289/impaired-local-intrinsic-immunity-to-sars-cov-2-infection-in-severe-covid-19#/)
-   Our project involves the usage of various tools from `omicsEye` to analyze data from a scRNA-seq study analyzing differences in gene expression across COVID-19 severity. 

### Reproducing our Results

- Scripts used for our final report can be found in `R/`.
  - `loading-data.R` performs data scraping and centralizing of large data files to a shared Google Drive. This contains functions that load in large amounts of data (that are now already processed and saved), and can be skipped when reproducing results

  - The shared Google Drive, containing all raw data files for use in `loading-data.R`, can be found [here](https://drive.google.com/drive/u/3/folders/0AJsoMBJBwr6MUk9PVA)

  - Our final team project report, `aims.pdf`, details our process, aims, and outputs from this investigation. It also contains the exact command line codes written to implement omeClust for cluster visualizations.
  
