from constants import ENCR_NOTES, DECR_FLAT_NOTES, DECR_SHARP_NOTES


class Note:
    def __init__(self, annotation):
        self.annotation = annotation
        self.encr = ENCR_NOTES[annotation]
        self.sign = None
        if len(annotation) == 2:
            self.sign = annotation[1]

    def __str__(self):
        return self.annotation

    def __repr__(self):
        return self.annotation

    def __eq__(self, other):
        return self.encr == other.encr

    def stack_semitones(self, num_semitones, sign = None):
        if sign == "#":
            decr = DECR_SHARP_NOTES
        elif sign == "b":
            decr = DECR_FLAT_NOTES
        else:
            if self.sign:
                decr = DECR_SHARP_NOTES if self.sign == "#" else DECR_FLAT_NOTES
            else:
                decr = DECR_SHARP_NOTES if num_semitones >= 0 else DECR_FLAT_NOTES  # FIXME
        return Note(decr[self.encr + num_semitones])

    def distance_to_note(self, note):
        interval = (note.encr - self.encr) % 12
        return interval


class ScaleInfo:
    def __init__(self, thirds, beginning):
        self.sign = None
        for third in thirds:
            if third.sign:
                self.sign = third.sign
        self.coherent_maj_thirds = thirds
        self.incoherent_maj_thirds = set()
        self.beginning = beginning
        self.end = None
        self.min_root = None
        if len(thirds) == 3:
            self.maj_root = thirds[1]
        else:
            self.maj_root = None
        self.mode = None

    def __repr__(self):
        return (f"(maj_root: {self.maj_root}, min_root: {self.min_root}, mode: {self.mode}, "
                f"coherent major thirds: {self.coherent_maj_thirds}, "
                f"incoherent major thirds: {self.incoherent_maj_thirds}, "
                f"beginning: {self.beginning}, end: {self.end})")

    def _is_match(self, new_third):
        if self.incoherent_maj_thirds:
            lower_notes_thirds = self.incoherent_maj_thirds
            for first, second in itertools.combinations(self.incoherent_maj_thirds, 2):
                fifth_stacking = stack_as_fifth(first, second, new_third)
                if fifth_stacking:
                    return fifth_stacking
            return None       
        elif self.coherent_maj_thirds:
            lower_notes_thirds = self.coherent_maj_thirds
            num_thirds = len(lower_notes_thirds)
            if num_thirds == 1:
                old = lower_notes_thirds[0]
                interval = old.distance_to_note(new_third)
                if interval == 2:
                    scale_root = old.stack_semitones(7)
                    stacked_lower_notes = [old, scale_root, new_third]
                elif interval == 10:
                    scale_root = new_third.stack_semitones(7)
                    stacked_lower_notes = [new_third, scale_root, old]
                elif interval == 7:
                    stacked_lower_notes = [old, new_third]
                elif interval == 5:
                    stacked_lower_notes = [new_third, old]
                else:
                    stacked_lower_notes = None
                return stacked_lower_notes
            elif num_thirds == 2:
                lower = lower_notes_thirds[0]
                higher = lower_notes_thirds[1]
                if new_third.distance_to_note(lower) == 7:
                    stacked_lower_notes = [new_third] + lower_notes_thirds
                elif higher.distance_to_note(new_third) == 7:
                    stacked_lower_notes = lower_notes_thirds + [new_third]
                else:
                    stacked_lower_notes = None
                return stacked_lower_notes
            else:
                if new_third in lower_notes_thirds:
                    return lower_notes_thirds
                else:
                    return None
        else:
            return [new_third]

    def update_key(self, new_third):
        thirds = self._is_match(new_third)
        if not thirds:
            if self.maj_root or len(self.incoherent_maj_thirds) > 2:
                return False
            else:
                self.coherent_maj_thirds = []
                self.incoherent_maj_thirds = set((self.coherent_maj_thirds + [new_third]))
        else:
            self.coherent_maj_thirds = thirds
            if len(thirds) == 3:
                self.maj_root = thirds[1]
            return True

