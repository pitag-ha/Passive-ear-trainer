ENCR_NOTES = {'C':0, 'C#':1, 'Db':1, 'D':2, 'D#':3, 'Eb':3, 'E':4, 'Fb':4, 'E#':5, 'F':5, 'F#':6, 'Gb':6, 'G':7, 'G#':8, 'Ab':8, 'A':9, 'A#':10, 'Bb':10, 'B':11, 'Cb':11, 'B#':0}

DECR_FLAT_NOTES = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab', 'A', 'Bb', 'B']
DECR_SHARP_NOTES = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
DECR_NOTES_PAIRS = [('C'), ('Db', 'C#'), ('D'), ('Eb', 'D#'), ('E'), ('F'), ('Gb', 'F#'), ('G'), ('Ab', 'G#'), 'A', ('Bb', 'A#'), ('B')]

def get_major_scale(root):
    return (map(lambda rel_note: (root + rel_note) % 12, [0, 2 ,4 ,5 ,7 ,9, 11]))

MAJOR_SCALES = {tuple(get_major_scale(ENCR_NOTES[root])) : root for root in ENCR_NOTES.keys()}

SCALE_DEG_IN_SEMITONES = {1:0, 2:2, 3:4, 4:5, 5:7, 6:9, 7:11}

BASS_FREQ = {'C':65.41, 'C#':69.30, 'D':73.42, 'D#':77.78, 'E':82.41, 'F':87.31, 'F#':92.50, 'G':98.00, 'G#':103.83, 'A':110.00, 'A#':116.54, 'B':123.47}
HIGH_PITH_FREQ = {'C':130.81, 'C#':138.59, 'D':146.83, 'D#':155.56, 'E':164.81, 'F':174.61, 'F#':185.00, 'G':196.00, 'G#':207.65, 'A':220.00, 'A#':233.08, 'B':246.94}