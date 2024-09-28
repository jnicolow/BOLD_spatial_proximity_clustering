import os
import pandas as pd
import requests

def species_from_species_namestxt():
    with open(os.path.join('bold', 'bold_queries', 'species_names.txt')) as f: lines = f.readlines()
    return(lines[0].split(','))

def get_sequences(species):
    print(f'downloading: {species}')
    url = f'https://www.boldsystems.org/index.php/API_Public/sequence?taxon={species}'
    response = requests.get(url)

    sequence_df = pd.DataFrame(columns=['species', 'BOLD_process_id', 'country', 'sequence'])
    if 'exceeded your allowed request quota' in response.text:
        print(response.text) # if you request to much they will say this
        return sequence_df
    sequence_infos = response.text.split('>')[1:] # the first element is just an xml header
    
    for sequence_info in sequence_infos:
        data_split = sequence_info.split('|') 
        species = data_split[1]
        ID = data_split[0]
        country = get_country(ID)
        # marker = dataSplit[2] # sometimes after the marker there is no | to solit on so in that case the markjer will have the sequence appended to it
        sequence_raw = data_split[-1] # some times the marker is missing so best to just take the last thing
        sequence_df = pd.concat([sequence_df, pd.DataFrame({'species':[species], 'BOLD_process_id':[ID], 'country':[country], 'sequence':[sequence_raw.split('\n')[1].strip('\r')]})])
    return sequence_df


def get_country(process_id):
    # return country from specimen api
    url = f'https://www.boldsystems.org/index.php/API_Public/specimen?ids={process_id}'


    response = requests.get(url)

    country_string = '</country'
    filtered_list = [item for item in response.text.split('>')[1:] if country_string in item] # [1:] because the first is always xml header
    country = filtered_list[0].replace(country_string, '')
    if country == '': country = 'Unknown'
    return(country)
