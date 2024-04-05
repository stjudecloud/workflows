from perl:5.38.2 as perl 

#from quay.io/biocontainers/pairix:0.3.7--py38h7cf9df2_5

#COPY --from=perl /usr/local/bin/perl  /usr/bin/perl
#COPY --from=perl /usr/local/lib/perl5 /usr/lib/perl5
#COPY --from=perl /usr/lib/x86_64-linux-gnu/perl /usr/lib/x86_64-linux-gnu/perl
#COPY --from=perl /usr/lib/x86_64-linux-gnu/libperl.so.5.36.0 /usr/lib/x86_64-linux-gnu/libperl.so.5.36.0
#COPY --from=perl /usr/lib/x86_64-linux-gnu/libperl.so.5.36 /usr/lib/x86_64-linux-gnu/libperl.so.5.36

RUN apt update && apt install -y samtools

WORKDIR /
RUN git clone https://github.com/4dn-dcic/pairix && cd pairix && git checkout 0.3.8 && make
ENV PATH="${PATH}:/pairix/bin:/pairix/util:/pairix/util/bam2pairs"


