FROM ubuntu:noble-20260113

RUN apt-get update \
    && apt-get upgrade -y \
    &&  apt-get install wget zip gcc zlib1g-dev libbz2-dev \
      liblzma-dev jq csvkit python3 python3-pip parallel \
      python3-csvkit python3-pysam python3-gtfparse -y \
    && rm -r /var/lib/apt/lists/*

COPY --from=scripts --chmod=777 util/calc_global_phred_scores.py /scripts/util/calc_global_phred_scores.py
COPY --from=scripts --chmod=777 util/calc_feature_lengths.py /scripts/util/calc_feature_lengths.py
COPY --from=scripts --chmod=777 util/check_FQs_and_RGs.py /scripts/util/check_FQs_and_RGs.py
COPY --from=scripts --chmod=777 htseq/calc_tpm.py /scripts/htseq/calc_tpm.py
COPY --from=scripts --chmod=777 data_structures/get_read_groups.py /scripts/data_structures/get_read_groups.py
