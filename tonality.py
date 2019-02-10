# start interactive IPython testing via (ear_trainer) sonja@Naoyo:~/Dropbox/RC/passive_ear_trainer$   ipython -i start_interactive_tonality_test.py

import collections
import itertools

from constants import ENCR_NOTES, DECR_FLAT_NOTES, DECR_SHARP_NOTES
from constants import MAJOR_SCALES

from chord_annotations import extract_root, extract_rootnote
from chord_annotations import extract_third_mod, extract_fifth_mod
from chord_annotations import extract_seventh_mod
from chord_annotations import get_notes
from datastructures import Note, ScaleInfo


def tonality_from_dominant(dominant):
    if len(dominant) > 1 and dominant[1] == 'b':
        decr = DECR_FLAT_NOTES
    else:
        decr = DECR_SHARP_NOTES
    arithm_repr_root = (ENCR_NOTES[dominant] - 7) % 12
    return decr[arithm_repr_root]


Scale = collections.namedtuple('Scale', ['root', 'arithmetic_scale', 'scale'])
PossibleScales = collections.namedtuple('PossibleScales', ['Scale', 'beginning', 'end'])


def find_major_third(chord):
    root, modalities = extract_root(chord)
    third, rest = extract_third_mod(modalities)
    fifth, rest = extract_fifth_mod(rest)
    if fifth == "perfect":
        if third == "major":
            base_of_third = root
        if third == "minor":
            base_of_third = DECR_FLAT_NOTES[(ENCR_NOTES[root] + 3) % 12]
        return Note(base_of_third)
    else:
        return None

# def interval(note1, note2):
#     note1_encr = ENCR_NOTES[note1]
#     note2_encr = ENCR_NOTES[note2]
#     interval_in_seminotes = (note2_encr - note1_encr) % 12
#     return interval_in_seminotes


def fifth(note, direction):
    note_encr = ENCR_NOTES[note]
    if direction == "above":
        return (DECR_SHARP_NOTES[(note_encr + 7) % 12])
    if direction == "below":
        return (DECR_SHARP_NOTES[(note_encr - 7) % 12])


def fifth_stacking(notes):
    if len(notes) == 1:
        return notes
    # if len(notes) == 2:
    #     interval = notes[0].distance_to_note(notes[1])
    #     if interval == 7:
    #         return notes
    #     elif interval == 5:
    #         return list(reversed(notes))
    #     else:
    #         return None
    for i, note in enumerate(notes):
        substack = fifth_stacking(notes[:i] + notes[i+1:])
        if substack:
            if note.distance_to_note(substack[0]) == 7:
                return [note] + substack
            elif substack[-1].distance_to_note(note) == 7:
                return substack + [note]
            else:
                return None


def third_analysis(list_chords):
    first_scale = ScaleInfo([], 0)
    scales = [first_scale]
    for i, chord in enumerate(list_chords):
        last_scale = scales[-1]
        if i < 4 or last_scale.beginning <= i-4:
            scale = last_scale
            possible_modulation = None
        else:
            possible_modulation = last_scale
            scale = scales[-2]
        if chord == 'N':
            continue
        new = find_major_third(chord)
        if new:
            if new in scale.coherent_maj_thirds:
                if new.sign:
                    scale.sign = new.sign
                continue
            matches_key = scale.update_key(new)  # note the side-effects of update_key
            if possible_modulation:
                matches_modulation = possible_modulation.update_key(new)  # note the side-effects of update_key
                if matches_modulation:
                    if new.sign:
                        possible_modulation.sign = new.sign
                else:
                    del scales[-1]
                    if not matches_key:
                        return False
            else:
                if not matches_key:
                    possible_modulation = ScaleInfo([new], i)
                    scales.append(possible_modulation)
    scales = list(filter(lambda scale: scale.beginning < (len(list_chords) - 4), scales))
    for scale in scales:
        scale.min_root = (scale.maj_root).stack_semitones(-3, scale.sign)
    partition = [scale.beginning for scale in scales] + [(len(list_chords) - 1)]
    bounds = list(zip(partition[:-1], partition[1:]))
    for i, span in enumerate(bounds):
        chords = [chord for chord in list_chords[span[0]: span[1]] if chord != 'N']
        roots = list(map(lambda chord: extract_root(chord)[0], chords))
        num_appearances = collections.Counter(roots)
        scale = scales[i]
        scale.end = span[1]
        if num_appearances[scale.min_root.annotation] > num_appearances[scale.maj_root.annotation]:
            scale.mode = "min"
        else:
            scale.mode = "maj"
    return scales
    return [(scale.beginning, scale.maj_root) for scale in scales if scale.beginning < (len(list_chords) - 4)]


def sign(roots):
    appearing_signs = {}
    for root in roots:
        if len(root) == 2:
            appearing_signs.add(root[1])
    if len(appearing_signs) == 0:  # if no sign appears, it doesn't matter, so we just choose sharp, why not?
        return "#"
    if len(appearing_signs) == 1:
        sign, = appearing_signs
        return sign
    if len(appearing_signs) > 1:
        return False

def get_root_of_correspondent_mode(root_maj_scale, list_chords):
    root_min_scale = root_maj_scale.stack_semitones(-3)  # FIXME: check sign!
    min_root = root_min_scale.annotation
    maj_root = root_maj_scale.annotation
    list_chords = filter(lambda chord: chord != 'N', list_chords)
    list_roots = list(map(lambda chord: extract_root(chord)[0], list_chords))
    maj_min_roots = list(filter(lambda root: root in [maj_root, min_root], list_roots))
    occurrences = collections.Counter(maj_min_roots)
    return occurrences.most_common(1)[0][0]


def find_tonality(song_analysis):

    # find last chord
    last_chord_index = -1
    while song_analysis[last_chord_index]['label'] == 'N':
        last_chord_index -= 1
    last_chord = song_analysis[last_chord_index]['label']
    last_root = extract_root(last_chord)[0]
    
    # return last_root

    # find the chord that is played the most
    roots_info = [extract_root(chord['label']) for chord in song_analysis if extract_root(chord['label'])]
    roots = [root for root, _ in roots_info]
    numb_appearance_root = collections.Counter(roots)
    most_appearing_root = numb_appearance_root.most_common(1)[0][0]
    
    # return most_appearing_root
    
    # find tonality via dominant7 chords:
    maj_dominant7_chords = {
        root 
        for root, modalities in roots_info
        if (
            extract_third_mod(modalities)[0] == "major"
            and extract_seventh_mod(modalities) == "minor"
            )
        }
    
    potential_min_dominant7_chords = {root for root, modalities in roots_info if extract_third_mod(modalities)[0] == "minor" and extract_seventh_mod(modalities) == "minor"}

    tonalities = {}
    if maj_dominant7_chords:
        tonalities = {tonality_from_dominant(dominant) for dominant in maj_dominant7_chords if tonality_from_dominant(dominant) in roots}
    # return tonalities

    list_chords = [chord_info['label'] for chord_info in song_analysis]

    # major thirds method:
    return third_analysis(list_chords)
    keys = third_analysis(list_chords)
    maj_scale_roots = [info[1] for info in keys]
    indices = [info[0] for info in keys] + [len(list_chords)]
    spans = list(zip(maj_scale_roots, indices[:-1], indices[1:]))
    return [(beginning, get_root_of_correspondent_mode(maj_scale_root, list_chords[beginning:end])) for maj_scale_root, beginning, end in spans]

    # find possible scales:
    scales = []
    first_index = 0
    possible_scales = MAJOR_SCALES
    for i, chord in enumerate(list_chords):
        if chord == 'N':
            continue
        possible_new_scales = {scale: root for scale, root in possible_scales.items() if get_notes(chord).issubset(set(scale))}
        if not possible_new_scales:
            s = set()# ipd.Audio(output, rate=song[1])



            for arithmetic_scale in possible_scales:
                s.add(Scale(possible_scales[arithmetic_scale], tuple(arithmetic_scale), tuple([DECR_SHARP_NOTES[note] for note in arithmetic_scale])))

            scales.append(PossibleScales(s, first_index, i-1))
            first_index = i
            possible_scales = {scale: root for scale, root in MAJOR_SCALES.items() if get_notes(chord).issubset({note for note in scale})}
        else:
            possible_scales = possible_new_scales

    return scales
    

        

#     sign = ''
#     if len(most_appearing_chord) > 1:
#         if most_appearing_chord[1] == '#' or most_appearing_chord[1] == 'b':
#             sign = most_appearing_chord[1]

#     tonality = most_appearing_chord[0] + sign
#    tonality = most_appearing_root
    
#     certainty = True
#     if last_chord != most_appearing_chord:
#         certeinty = False
        



def print_scales(scales, song_analysis):
    list_chords = [chord_info['label'] for chord_info in song_analysis]
    for ss in scales:
        print(f"From chord nr {ss.beginning} till chord nr {ss.end} the possible scales are:")
        for s in ss.Scale:
            print("Scale: ", s.scale) 
        next_chord = list_chords[ss.end + 1]
        print(f"The next chord is {next_chord}. Its (basic) notes are {[DECR_SHARP_NOTES[note] for note in get_notes(next_chord)]}")
        print("---------------")