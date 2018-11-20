FROM continuumio/miniconda

# Grab requirements.txt.
ADD ./webapp/requirements.txt /tmp/requirements.txt

#RUN /bin/sh -c pip install -qr /tmp/requirements.txt

# Install dependencies
#RUN pip install --upgrade pip
#RUN pip install --upgrade setuptools
RUN pip install -qr /tmp/requirements.txt

#ADD ./vamp_setup/VAMP-0.9.0 /tmp/vamp_setup
#RUN python /tmp/vamp_setup/setup.py install

# Add our code
ADD ./webapp /opt/webapp/
WORKDIR /opt/webapp

RUN conda install scikit-learn
RUN conda install -c conda-forge librosa

CMD gunicorn --bind 0.0.0.0:$PORT wsgi