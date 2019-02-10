#%%
#%matplotlib inline
import librosa, librosa.display
import IPython.display as ipd
import vamp
import matplotlib.pyplot as plt
import numpy as np
import scipy.fftpack as fft
import re
import collections
import pprint
pp = pprint.PrettyPrinter(indent=4)
pprint = pp.pprint
import vamp

#%%
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/internationale.mp4'
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/Autumn Leaves - Yenne Lee plays 2004 Pepe Romero Jr.-HxGT5z6d-GA.m4a'
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/Beatles - I wanna hold your hand - live.m4a'
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/Beatles - I wanna hold your hand -studio.m4a'
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/Dumb Ways to Die-IJNR2EpS0jw.m4a'
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/Queen - I Want to Break Free (Official Lyric Video)-WUOtCLOXgm8.m4a'
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/Selena Gomez, Marshmello - Wolves (Official Music Video)-cH4E_t3m3xM.m4a'
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/Ska-p Los hijos bastardos de la globalizacion con Letra-upnPasIYeMc.m4a'
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/The Killers - Human-RIZdjT1472Y.m4a'
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/THE MUFFIN SONG (asdfmovie feat. Schmoyoho)-LACbVhgtx9I.m4a'
#path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/blowing_in_the_wind.m4a'
path = '/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/applaus_sportfreundestiller.m4a'

song = librosa.load(path)
song_analysis = vamp.collect(song[0], song[1], "nnls-chroma:chordino")['list']
pprint(list(enumerate(song_analysis)))

