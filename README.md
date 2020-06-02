# Merino - Pre-publication Release

Pre-publication version of data Flow for converting raw PR500 data into processed results

## Setting up the environment for development:

1.  Create a local conda environment

    ```<path to conda bin>/conda create --name merino -c bioconda python=2.7.11 numpy=1.11.2 pandas=0.20 h5py=2.6.0 requests=2.13.0 cmapPy scipy setuptools pathlib patsy matplotlib seaborn statsmodels yaml pyyaml jinja2```

2.  Activate the environment:
    ```source <path to conda bin>/activate merino```
3.  Add the merino code directory to the python path: 
    
    `cd <directory above path to merino code>`
    
    `python <path to merino code>/setup_merino.py develop`
    
4.  Run tests to verify setup
    
    `cd <path to merino code>`
    
    `for a in test_*.py; do python $a; echo $a; read; done`
    
## Running Merino

### Assemble
Assemble combines a data CSV with the appropriate annotations to create the GCT file used in downstream processing. 

For backwards compatibility with DP78 assay, command:
`python assemble.py -assay_type DP78 -davepool_id_csv_filepairs DP7 {DP7_csv_path} DP8 {DP8_csv_path} -plate_map_path {plate_map_path} -out {out_location}`

For all other assay types : 
`python assemble.py -csv_filepath {csv_filepath} -plate_map_path {plate_map_path} -out {out_location}`

Additional arguments can be used to override defaults. By default, assemble determine assay type by querying the CLUE API with the bead_batch (the last token in the plate name).
The pipeline config, cell_set_definition_files, and analyte_mapping_files definitive resources are stored in s3. These too can be overriden with optional arguments.
