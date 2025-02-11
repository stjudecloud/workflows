FROM aidenlab/juicer:v2.0.1 AS builder

FROM alpine/java:22

RUN mkdir -p /HiC-Pro_3.0.0/bin/utils/ \
    && cd /HiC-Pro_3.0.0/bin/utils/ \
    && wget https://github.com/nservant/HiC-Pro/raw/refs/tags/v3.0.0/bin/utils/hicpro2juicebox.sh \
    && chmod a+x /HiC-Pro_3.0.0/bin/utils/hicpro2juicebox.sh

COPY --from=builder /opt/juicer /opt/juicer
RUN chmod a+rwx /opt/juicer/CPU/common/juicer_tools.jar

RUN apk add --no-cache bash

ENTRYPOINT [ "bash" ]