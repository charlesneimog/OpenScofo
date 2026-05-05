import OpenScofo
import librosa

scofo = OpenScofo.OpenScofo(48000, 2048, 512)
config = scofo.get_configuration()

config.chroma_size = 24

scofo.set_configuration(config)
config2 = scofo.get_configuration()
print(config2.chroma_size)
