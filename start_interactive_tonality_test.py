from song_paths import PATH
from pprint import pprint
from song_analysis import get_song_info

_, _, song_analysis = get_song_info(PATH)
pprint(list(enumerate(song_analysis)))  

from tonality import find_tonality
from tonality import print_scales