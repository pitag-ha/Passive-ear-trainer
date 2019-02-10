import librosa
import numpy as np

cow = librosa.load('/home/sonja/Dropbox/RC/passive_ear_trainer/Songs/Kuh muht-TdheW61w4Co.m4a')


def extend_short_rec(short_rec, desired_num_samples, desired_vol):
    bounds = [1113687, 1148687] #FIXME
    current_vol = max(short_rec[0])
    
    sample = last_not_0_before_par(short_rec, bounds[1])
    for _ in range(50):
        sample = sign_change(short_rec, sample-1)
    end_copy = sample

    for _ in range(1000):
        sample = sign_change(short_rec, sample-1)
    beginning_copy = sample

    extension = short_rec[0][bounds[0]:end_copy]
    while len(extension) < desired_num_samples:
        extension = np.concatenate((extension, short_rec[0][beginning_copy:end_copy]))
        
    extension = extension/current_vol*desired_vol
    
    return extension[:desired_num_samples]


def cow_generator(amp, rate):
    if rate != cow[1]:
        return False
    def f(num_samples):
        return extend_short_rec(cow, num_samples, amp)
    return f 


def last_not_0_before_par(rec, par):
    ind = par
    while rec[0][ind] == 0:
        ind -= 1
    return ind


def first_0_before_par(rec, par):
    ind = par
    while rec[0][ind] != 0:
        ind -= 1
    return ind


def sign_change(rec, par):
    ind = par
    while rec[0][ind] * rec[0][par] > 0:
        ind -= 1
    return ind


def note_by_num_samples(freq, num_samples, amp, rate):
    data = np.array([np.sin(2*np.pi*freq*librosa.samples_to_time(sample, rate))*amp for sample in range(num_samples)], dtype=np.float32)
    return data


def synth_generator(freq, amp, rate):
    def f(num_samples):
        return note_by_num_samples(freq, num_samples, amp, rate)
    return f


