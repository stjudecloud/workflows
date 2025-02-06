FROM nservant/hicpro:3.0.0 AS hicpro
FROM aidenlab/juicer:1.0.13

COPY --from=hicpro /HiC-Pro_3.0.0/bin/utils/hicpro2juicebox.sh /HiC-Pro_3.0.0/bin/utils/hicpro2juicebox.sh
RUN chmod a+rwx /opt/juicer-1.6.2/CPU/common/juicer_tools.1.7.6_jcuda.0.8.jar

ENTRYPOINT [ "bash" ]