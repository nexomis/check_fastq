FROM python:3.12
 
COPY ./src/check_fastq.py /usr/local/bin/check_fastq.py

RUN chmod 755 /usr/local/bin/check_fastq.py \
  && pip install psutil

ENTRYPOINT []