#!/usr/bin/env python3
# coding: utf-8


# %%
import sys
sys.path.append(".")
print(sys.path)



#%%
%load_ext autoreload
%autoreload 2




# %%


# import librosa, librosa.display
import IPython.display as ipd
import pprint
pp = pprint.PrettyPrinter(indent=4)
pprint = pp.pprint

from tonality import find_tonality
from generation_audio_to_add import generate_output
from song_paths import PATH
from song_analysis import get_song_info




# %%


song_array, song_sample_rate, song_analysis = get_song_info(PATH)
pprint(list(enumerate(song_analysis)))





# %%


scale_deg = 1  # 1,2,3,4,5,6 or 7
timber = 'high_pitch_synth'  # 'low_pitch_synth', 'high_pitch_synth' or 'cow'
amp = 1  # 1/4, 1/2, 1, 2 

tonality = find_tonality(song_analysis)
# if not certainty:
#     print('The tonality is probably wrong, sorry for that!')
print(tonality)


# %%


output = generate_output(scale_deg, timber, amp, song_array, song_sample_rate, song_analysis, tonality)
print(output[:20])
ipd.Audio(output, rate=song_sample_rate)


# # %%


# output = song[0] + to_add
# print(len(to_add))
# print(len(song[0]))
# print(tonality)
# bounds = find_bounds_of_chord(song_analysis, chord_regex = re.escape(tonality)+r'.*')
# beginning, end = bounds[4]
# ipd.Audio(output, rate=song[1])


