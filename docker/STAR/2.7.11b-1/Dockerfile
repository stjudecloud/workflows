FROM debian:10 as python

RUN apt-get update \
    && apt-get --yes install --no-install-recommends python3 \
    && rm -rf /var/lib/apt/lists/*

FROM quay.io/biocontainers/star:2.7.11b--h43eeafb_3

COPY --from=python /lib/x86_64-linux-gnu/libexpat.so.1 /lib/x86_64-linux-gnu/
COPY --from=python /lib/x86_64-linux-gnu/libz.so.1 /lib/x86_64-linux-gnu/
COPY --from=python /usr/lib/python3.7/ /usr/local/lib/python3.7/
COPY --from=python /usr/bin/python3 /usr/local/bin/

COPY --from=scripts star/sort_star_input.py /home

ENTRYPOINT ["STAR"]
