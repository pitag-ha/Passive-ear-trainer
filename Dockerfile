FROM continuumio/miniconda3

# Grab requirements.txt. In case of using miniconda, have the following line uncommented:
ADD ./webapp/requirements_py2.txt /tmp/requirements.txt

#In case of using miniconda3, have the following line uncommented:
ADD ./webapp/requirements_py3.txt /tmp/requirements.txt

# Install dependencies
RUN apt-get -y update

#Necessary to pip install vamp:
RUN apt-get -y install --reinstall build-essential=12.3

RUN pip install --upgrade pip
RUN pip install numpy==1.15.4
RUN pip install -r /tmp/requirements.txt

#Grab vamp plugins
add ./vamp_plugins/ /usr/local/lib/vamp/

# Grab code and songs
ADD ./webapp /opt/webapp/
ADD ./Songs /opt/Songs/
WORKDIR /opt/webapp

RUN conda install scikit-learn
RUN conda install -c conda-forge librosa

CMD gunicorn --bind 0.0.0.0:$PORT wsgi