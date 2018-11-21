from flask import Flask, jsonify
from flask import render_template
from sklearn import datasets, svm
import numpy as np
import vamp
import librosa
#import os 

#os.environ['VAMP_PATH'] = "/opt/webapp/vamp"

app = Flask(__name__)

# Load Dataset from scikit-learn.
digits = datasets.load_digits()

@app.route('/')
def hello():
    # clf = svm.SVC(gamma=0.001, C=100.)
    # clf.fit(digits.data[:-1], digits.target[:-1])
    # prediction = clf.predict(digits.data[-1:])

    # return jsonify({'prediction': repr(prediction)})

    #vamp_plugins = vamp.list_plugins()

    song = librosa.load('/opt/Songs/Kuh muht-TdheW61w4Co.m4a')
    some_samples = song[0][10000:10020]
    vamp_plugins = vamp.list_plugins()
    #return 'Some samples from my song: {} /n My vamp plugins: {}'.format(some_samples, vamp_plugins)
    #return 'hello'

    return render_template('test.html', samples=some_samples, plugins=vamp_plugins)

if __name__ == '__main__':
    app.run(host='0.0.0.0')
