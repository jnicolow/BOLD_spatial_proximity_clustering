import os
import pandas as pd
import requests

def get_sequences(species):
    print(species)
    url = f'https://www.boldsystems.org/index.php/API_Public/sequence?taxon={species}'
    response = requests.get(url)
    sequnceInfos = response.text.split('>')[1:]
    sequenceDf = pd.DataFrame(columns=['Species', 'BOLDProcessId', 'Sequence'])
    for sequenceInfo in sequnceInfos:
        dataSplit = sequenceInfo.split('|') 
        species = dataSplit[1]
        ID = dataSplit[0]
        # marker = dataSplit[2] # sometimes after the marker there is no | to solit on so in that case the markjer will have the sequence appended to it
        sequenceRaw = dataSplit[-1] # some times the marker is missing so best to just take the last thing
        sequenceDf = pd.concat([sequenceDf, pd.DataFrame({'Species':[species], 'BOLDProcessId':[ID], 'Sequence':[sequenceRaw.split('\n')[1].strip('\r')]})])
    return sequenceDf