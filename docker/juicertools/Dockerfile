FROM nservant/hicpro:3.0.0 AS hicpro
FROM aidenlab/juicer:v2.0.1

COPY --from=hicpro /HiC-Pro_3.0.0/bin/utils/hicpro2juicebox.sh /HiC-Pro_3.0.0/bin/utils/hicpro2juicebox.sh
RUN chmod a+rwx /opt/juicer/CPU/common/juicer_tools.jar

ENTRYPOINT [ "bash" ]