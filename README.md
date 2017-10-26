# Merino

Data Flow for converting raw PR500 data into processed results

## Setting up the environment for development:

0.  Create a local conda environment
    0.  <path to conda bin>/conda create --name merino python=2.7.11 numpy=1.11.2 pandas=0.18 h5py=2.6.0 requests==2.13.0 cmappy scipy setuptools
0.  Activate the environment:
    0.  source \<path to conda bin>/activate merino
0.  Add the merino code directory to the python path:
    0.  cd \<directory above path to merino code>
    0.  python \<path to merino code>/setup_merino.py develop
0.  Add GCToo code to path TODO fix this - needs to use cmapPy
    0.  cd \<directory above path to GCToo e.g. pestle/cmap/io>
    0.  python GCToo/setup_GCToo.py deveop
0.  Run tests to verify setup
    0.  cd \<path to merino code>
    0.  for a in test_*.py; do python $a; echo $a; read; done
    
    
TODO's
0. fix broken test:  test_assemble_no_davepool.py test_main
0. get environment working with mysql / caldaia / lims to register plates

