# omega
dN/dS analysis

## Run preprocessing
`$ python src/preprocessing/preprocessing.py test/input.json`

## Input
- File with somatic mutations:
    - Format:
        - The file should contain at least these 5 columns: CHROM, POS, REF, ALT, SAMPLE_ID

- Depths:
    - Format:
        - The first two columns of the file have to be informing about the chromosome and the position.
        - The next columns must have the name of the sample in the header and for each row, the corresponding value of depth in that chromosome and position.
        - Positions without coverage have to be filled with 0s.

- File defining the regions/sites where you want to focus the analysis:
    - BED file defining the regions sequenced: **(always do this in the first run)**
    This will make the execution start from the beginning and generate a file with all possible sites that will need to be annotated by VEP.
        - Format: BED6-like (not BED12)
        - The start and the end of the intervals are included as part of the analysis.

    - VEP output file:
        - This must be the file downloaded from Ensembl VEP execution using as input the file of all possible sites created by the `omega` preprocessing module, otherwise the downstream steps might not work.

    - VEP annotated sites postprocessed:
        - If you have already run the full preprocessing, this is the file that you should keep and then as long as you do not change the regions where you focus the analysis you can keep using this same file.
        - You could try to generate this file (see test/output/preprocessing/*.annotation_summary.tsv) or a file with the same information but the easiest way to proceed is to run the `omega` preprocessing module.


Make sure that the chromosome names are defined in the same way across the three files. For safety reasons, the "chr" prefix is enforced in all the three files.


### Running Ensembl VEP

1. Take the file provided by omega preprocessing module.

2. Upload it to the server see here.

![Ensembl VEP upload file](./docs/img/ChooseFile_button.png)

3. Select EnsemblVEP parameters:
    - Leave everything as default except for the options to choose which variants to output. (see picture below)

![Ensembl VEP select output](./docs/img/EnsemblVEP_select_output.png)

4. Download the results, back to your computer/cluster

![Download Ensembl VEP results](./docs/img/DownloadEnsemblVEPoutput.png)

**NOTE:** If you are using a HPC clusterAlternatively you can right click on this button and copy the link. Then use `wget` to download the file directly to your desired location in the cluster, like this:

    wget "https://www.ensembl.org/Homo_sapiens/Download/Tools/VEP?format=txt;tl=....." -O regions_sites.VEP_annotated.tsv

5. Provide the full path to this file back to the omega preprocessing execution that will be waiting for this.

## 