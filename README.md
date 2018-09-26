# Merino - !! Pre-publication !!

Pre-publication version of data Flow for converting raw PR500 data into processed results

## Setting up the environment for development:

1.  Create a local conda environment
    
    ```<path to conda bin>/conda create --name merino -c bioconda python=2.7.11 numpy=1.11.2 pandas=0.18 h5py=2.6.0 requests=2.13.0 cmapPy scipy setuptools```
2.  Activate the environment:
    ```source <path to conda bin>/activate merino```
3.  Add the merino code directory to the python path:
    
    ```cd <directory above path to merino code>
    python <path to merino code>/setup_merino.py develop```
    
4.  Run tests to verify setup
    ```cd <path to merino code>
    for a in test_*.py; do python $a; echo $a; read; done```
    
## TODO's
1. fix broken test:  test_assemble_no_davepool.py test_main
1. get environment working with mysql / caldaia / lims to register plates
