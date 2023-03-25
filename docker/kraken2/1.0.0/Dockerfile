FROM quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2

RUN sed -i 's;$full_path =~ s#^ftp://${qm_server}${qm_server_path}/##;$full_path =~ s#^https://${qm_server}${qm_server_path}/##;g' /usr/local/libexec/rsync_from_ncbi.pl