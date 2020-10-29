FROM clearlinux/numpy-mp:latest
COPY . /bench
RUN pip install -r requirements.txt
ENTRYPOINT ["docker-entrypoint.sh", "./rngbias.py"]
