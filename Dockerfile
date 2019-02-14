FROM continuumio/miniconda
MAINTAINER Desiree Davison <cmap-soft@broadinstitute.org>
LABEL merino.pipeline.clue.io.version="0.0.1"
LABEL merino.pipeline.clue.io.vendor="Connectivity Map"
RUN mkdir -p /cmap/bin && \
mkdir -p /cmap/merino/ && \
cd /cmap/merino/ && \
conda create -y --name merino -c bioconda python=2.7.11 numpy=1.11.2 pandas=0.20 h5py=2.6.0 requests=2.13.0 cmapPy=3.2.0 scipy setuptools argparse=1.4.0 boto3 pathlib patsy matplotlib seaborn statsmodels yaml pyyaml && \
cd /cmap && \
git clone -b automation https://github.com/cmap/merino.git && \
cp /cmap/merino/assemble/batch_assemble.sh /cmap/bin/assemble && \
cp /cmap/merino/card/batch_card.sh /cmap/bin/card && \
cp /cmap/merino/normalization/batch_weave.sh /cmap/bin/weave && \
cp /cmap/merino/normalization/batch_mk_build.sh /cmap/bin/mk_build  && \
cp /cmap/merino/build_summary/qc.sh /cmap/bin/qc

WORKDIR /cmap/bin
env PATH /cmap/bin:$PATH
RUN ["chmod","-R", "+x", "/cmap/merino/"]
RUN ["chmod","-R", "+x", "."]
ENTRYPOINT ["/bin/bash"]
CMD ["-help"]