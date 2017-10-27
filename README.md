# Merino

Data Flow for converting raw PR500 data into processed results

## Setting up the environment for development:

1.  Create a local conda environment
    1.  \<path to conda bin>/conda create --name merino python=2.7.11 numpy=1.11.2 pandas=0.18 h5py=2.6.0 requests==2.13.0 cmappy scipy setuptools
1.  Activate the environment:
    1.  source \<path to conda bin>/activate merino
1.  Add the merino code directory to the python path:
    1.  cd \<directory above path to merino code>
    1.  python \<path to merino code>/setup_merino.py develop
1.  Run tests to verify setup
    1.  cd \<path to merino code>
    1.  for a in test_*.py; do python $a; echo $a; read; done
    
## TODO's
1. fix broken test:  test_assemble_no_davepool.py test_main
1. get environment working with mysql / caldaia / lims to register plates

