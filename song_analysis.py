import librosa
import vamp

def get_song_info(path):
    song_array, song_sample_rate = librosa.load(path)
    song_analysis = vamp.collect(song_array, song_sample_rate, "nnls-chroma:chordino")['list']
    return song_array, song_sample_rate, song_analysis
