FROM python:3.12-bookworm

#Install necessary packages for daetools
RUN apt-get update && \
    apt-get install --yes libgl1-mesa-glx libgfortran5

#Download and install daetools
RUN pip install https://sourceforge.net/projects/daetools/files/daetools/2.3.0/daetools-2.3.0-gnu_linux-x86_64.zip
