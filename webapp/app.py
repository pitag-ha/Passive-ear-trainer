from flask import Flask, jsonify
from sklearn import datasets, svm
import numpy as np
import vamp
import librosa
import os 

os.environ['VAMP_PATH'] = "/opt/webapp/vamp"

app = Flask(__name__)

# Load Dataset from scikit-learn.
digits = datasets.load_digits()

@app.route('/')
def hello():
    # clf = svm.SVC(gamma=0.001, C=100.)
    # clf.fit(digits.data[:-1], digits.target[:-1])
    # prediction = clf.predict(digits.data[-1:])

    # return jsonify({'prediction': repr(prediction)})

    vamp_plugins = vamp.list_plugins()

    #song = librosa.load('/home/sonja/Dropbox/RC/Passive ear trainer/Songs/Ska-p Los hijos bastardos de la globalizacion con Letra-upnPasIYeMc.m4a')
    #some_samples = song[0][100:120]
    return f'{vamp_plugins}'

if __name__ == '__main__':
    app.run(host='0.0.0.0')
