from note_generator import synth_generator
from note_generator import cow_generator

import numpy as np
import librosa
import re

from constants import BASS_FREQ
from constants import HIGH_PITH_FREQ
from constants import ENCR_NOTES
from constants import DECR_FLAT_NOTES
from constants import DECR_SHARP_NOTES
from constants import SCALE_DEG_IN_SEMITONES


def find_bounds_of_chord(song_analysis, chord_regex, rate):
    #bounds_chord = [librosa.time_to_samples(t, song[1]) for t in timestamps_chord]

    pairs = list(zip(song_analysis[:-1], song_analysis[1:]))
    time_bounds = [(float(chord['timestamp']), float(next_chord['timestamp'])) for chord, next_chord in pairs if re.match(chord_regex, chord['label'])]
    bounds = [librosa.time_to_samples(t, rate)for t in time_bounds]
    return bounds


def generate_audio_to_add(song_analysis, scale_duration, chord, timber, amp, rate):
    if timber == 'cow':
        note_generator = cow_generator(amp, rate)
        if not note_generator:
            return 'Sorry, the cow is tired of mooing. Please, choose a different timber!'
    else:
        if timber == 'low_pitch_synth':
            freq = BASS_FREQ[chord]
        else:
            freq = HIGH_PITH_FREQ[chord]
        note_generator = synth_generator(freq, amp, rate)
    
    bounds_chord = find_bounds_of_chord(song_analysis, re.escape(chord)+r'.*', rate)
    add_to_song = np.zeros(scale_duration)
    for beginning, end in bounds_chord:
        if end > scale_duration:
            end = scale_duration
        add = note_generator(end - beginning)
        add_to_song[beginning:end] = add
        
    return add_to_song


#given a tonality and a scale degree, find the chord corresponding to the scale degree in the given tonality
def absolute_chord(tonality, scale_deg):
    encr_chord = (ENCR_NOTES[tonality] + SCALE_DEG_IN_SEMITONES[scale_deg]) % 12
    if len(tonality) > 1 and tonality[1] == 'b':
        chord = DECR_FLAT_NOTES[encr_chord]
    else:
        chord = DECR_SHARP_NOTES[encr_chord]
    return chord


def find_amplitude_rel_song(song_array, strength):
    song_intensity = max(song_array)
    return song_intensity*strength


def generate_output(scale_deg, timber, strength, song_array, song_sample_rate, song_analysis, scales):
    for scale in scales:
        if scale.mode = "maj":
            highlighted_chord = absolute_chord(scale.maj, scale_deg)
        else:
            highlighted_chord = absolute_chord(scale.min, scale_deg)
        length = len(song_array)
        amp = find_amplitude_rel_song(song_array, strength)
        to_add = generate_audio_to_add(song_analysis, length, highlighted_chord, timber, amp, song_sample_rate)
    output_array = song_array + to_add
    return output_array
