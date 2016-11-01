# PRISM Data Flow

Data Flow for converting raw PRISM data into processed results

## Setting up the environment for development:

0.  Create a local conda environment
    0.  conda env create -f \<path to prism_pipeline code>/prism_pipeline.conda
     Note: This will only work with conda version 3.19
0.  Activate the environment:
    0.  source \<path to conda bin>/activate prism_pipeline
0.  Add the prism pipeline code directory to the python path:
    0.  cd \<directory above path to prism_pipeline code>
    0.  python \<path to prism_pipeline code>/setup_prism_pipeline.py develop
0.  Add GCToo code to path
    0.  cd \<directory above path to GCToo e.g. l1ktools/python/broadinstitute_cmap/io>
    0.  python GCToo/setup_GCToo.py develop
0.  Run tests to verify setup
    0.  cd \<path to GCToo e.g. pestle/cmap/io/GCToo>
    0.  for a in test_*.py; do python $a; echo $a; read; done
    0.  cd \<path to prism_pipeline code>
    0.  for a in test_*.py; do python $a; echo $a; read; done
