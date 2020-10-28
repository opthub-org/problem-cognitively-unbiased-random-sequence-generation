FROM clearlinux/numpy-mp
COPY . /bench
RUN pip install -r requirements.txt
ENTRYPOINT ["docker-entrypoint.sh", "python3", "rngbias.py"]
