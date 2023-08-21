# omega
dN/dS analysis

## Run preprocessing:
$ python src/preprocessing/preprocessing.py 
Use this file to run Ensembl VEP:
./test/input/full_case/KidneyPanel.sites4VEP.tsv
Introduce the path to the EnsemblVEP annotated file:
## Input
- BED file defining the regions sequenced:
    - Format: BED6-like (not BED12)
    - The start and the end of the intervals are included as part of the analysis.

- File with somatic mutations:
    - Format:
        - The file should contain at least these 5 columns: CHROM, POS, REF, ALT, SAMPLE_ID

- Depths:
    - Format:
        - The first two columns of the file have to be informing about the chromosome and the position.
        - The next columns must have the name of the sample in the header and for each row, the corresponding value of depth in that chromosome and position.
        - Positions without coverage have to be filled with 0s.

Make sure that the chromosome names are defined in the same way across the three files.

Check the presence or absence of the "chr" prefix.

### Running Ensembl VEP



## 