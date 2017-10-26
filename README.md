# Merino

Data Flow for converting raw PR500 data into processed results

## Setting up the environment for development:

0.  Create a local conda environment
    0.  conda env create -f \<path to merino code>/merino.conda
0.  Activate the environment:
    0.  source \<path to conda bin>/activate merino
0.  Add the merino code directory to the python path:
    0.  cd \<directory above path to merino code>
    0.  python \<path to merino code>/setup_merino.py develop
0.  Add GCToo code to path TODO fix this - needs to use cmapPy
    0.  cd \<directory above path to GCToo e.g. pestle/cmap/io>
    0.  python GCToo/setup_GCToo.py deveop
0.  Run tests to verify setup
    0.  cd \<path to GCToo e.g. pestle/cmap/io/GCToo>
    0.  for a in test_*.py; do python $a; echo $a; read; done
    0.  cd \<path to merino code>
    0.  for a in test_*.py; do python $a; echo $a; read; done
