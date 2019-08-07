FROM cmap/base-merino:latest
MAINTAINER Desiree Davison <cmap-soft@broadinstitute.org>
LABEL merino.pipeline.clue.io.version="0.0.1"
LABEL merino.pipeline.clue.io.vendor="Connectivity Map"

RUN cd /cmap && \
git clone https://github.com/cmap/merino.git && \
cp /cmap/merino/merino/assemble/batch_assemble.sh /cmap/bin/assemble && \
cp /cmap/merino/merino/card/batch_card.sh /cmap/bin/card && \
cp /cmap/merino/merino/normalization/batch_weave.sh /cmap/bin/weave && \
cp /cmap/merino/merino/normalization/batch_mk_build.sh /cmap/bin/mk_build  && \
cp /cmap/merino/merino/build_summary/batch_plate_qc.sh /cmap/bin/plate_qc && \
cp /cmap/merino/merino/build_summary/build_qc.sh /cmap/bin/qc

WORKDIR /cmap/bin
env PATH /cmap/bin:$PATH
RUN ["chmod","-R", "+x", "/cmap/merino/"]
RUN ["chmod","-R", "+x", "."]
ENTRYPOINT ["/bin/bash"]

CMD ["-help"]

