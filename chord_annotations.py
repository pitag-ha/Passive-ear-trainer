from constants import ENCR_NOTES
from datastructures import Note

import re

def extract_rootnote(chord):
    root = chord[:2] if (len(chord) > 1 and chord[1] in ['b', '#']) else chord[0]
    if root == "N":
        return
    next = len(root)
    if root[0] in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
        return Note(root), chord[next:]
    else:
        raise Exception("This chord doesn't have a root note.")


def extract_root(chord):
    root = chord[:2] if (len(chord) > 1 and chord[1] in ['b', '#']) else chord[0]
    if root == "N":
        return
    next = len(root)
    if root[0] in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
        return root, chord[next:]
    else:
        raise Exception("This chord doesn't have a root note.")


def extract_third_mod(piece_chord):
    third = "major"
    length_token = 0
    if piece_chord:
        if piece_chord[0] == 'm' and piece_chord[:3] != "maj":
            third = 'minor'
            length_token = 1
        elif piece_chord[:3] == "sus":
            third = piece_chord[: 4]
            length_token = 4
    return third, piece_chord[length_token:]


def extract_fifth_mod(piece_chord):
    fifth = "perfect"
    length_token = 0
    if piece_chord and piece_chord[:3] in ["dim", "aug"]:
        fifth = piece_chord[:3]
        length_token = 3
    return fifth, piece_chord[length_token:]


def extract_bass(root, piece_chord):
    bass_regex = re.compile(r"\/([A-G](?:#|b)?)$")
    bass_specification = re.findall(bass_regex, piece_chord)
    if bass_specification:
        bass = bass_specification[0]
        rest = piece_chord[: -(len(bass) + 1)]
    else:
        bass = root
        rest = piece_chord
    return bass, rest


def extract_seventh_mod(piece_chord):
    seventh = None 
    re_others = re.compile(r"(?:maj)?(?:7)")
    seventh = re.findall(re_others, piece_chord)
    meaning_abbreviations = {"maj7": "major", "7": "minor"}
    if seventh:
        return meaning_abbreviations[seventh[0]]
    else:
        return
    

def get_notes(chord, with_bass=False):
    seminotes_third = {"major": 4, "minor": 3, "sus2": 2, "sus4": 5}
    seminotes_fifth = {"perfect": 7, "dim": 6, "aug": 8}
    seminotes_seventh = {"major": 11, "minor": 10}
    root, rest = extract_root(chord)
    third_mod, rest = extract_third_mod(rest)
    fifth_mod, rest = extract_fifth_mod(rest)
    bass, rest = extract_bass(root, rest)
    seventh_mod = extract_seventh_mod(rest)
    notes = {ENCR_NOTES[root], (ENCR_NOTES[root] + seminotes_third[third_mod]) % 12, (ENCR_NOTES[root] + seminotes_fifth[fifth_mod]) % 12}
    if with_bass and bass:
        print("BASS")
        notes.add(ENCR_NOTES[bass])
    if seventh_mod:
        notes.add(((ENCR_NOTES[root] + seminotes_seventh[seventh_mod])) % 12)
    return notes


def _examples():
    print("extract_seventh_mod of bamaj7b10/G#: ", extract_seventh_mod("bamaj7b10/G#"))
    print("extract_bass of D and sus2maj7b10/G#: ", extract_bass("D", "sus2maj7b10/G#"))
    print("extract_fifth_mod of maj7b10/G#: ", extract_fifth_mod('maj7b10/G#'))
    print("get_notes of Bbmmaj7 : ", get_notes("Bbmmaj7"))


if __name__ == "__main__":
    _examples()

