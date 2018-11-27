FROM continuumio/miniconda3

# Grab requirements.txt. In case of using miniconda, have the following line uncommented:
#ADD ./webapp/requirements_py2.txt /tmp/requirements.txt

#Necessary to pip install vamp:
RUN apt-get -y update && \
    apt-get -y install --reinstall build-essential=12.3 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN conda install -c conda-forge librosa

#Grab requirements: In case of miniconda, with ./webapp/requirements_py.txt /tmp/requirements.txt; in case of miniconda3, with ./webapp/requirements_py3.txt /tmp/requirements.txt
ADD ./webapp/requirements_py3.txt /tmp/requirements.txt

RUN pip install --upgrade pip && \
    pip install numpy==1.15.4 && \
    pip install -r /tmp/requirements.txt

#Grab vamp plugins
add ./vamp_plugins/ /usr/local/lib/vamp/

# Grab code and songs
ADD ./webapp /opt/webapp/
ADD ./Songs /opt/Songs/
WORKDIR /opt/webapp

CMD gunicorn --bind 0.0.0.0:$PORT wsgi