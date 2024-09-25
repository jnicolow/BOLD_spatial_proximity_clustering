import os
import pandas as pd
import numpy as np
import torch
import esm

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

# Get embeddings
def get_embeddings_batch(sequences):
    try:
        _, _, batch_tokens = batch_converter(list(enumerate(sequences)))

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True) # use the representation from the 33rd layer

        sequenceEmbedings = {}

        for i, seq in enumerate(sequences):
            trimmedRep = results['representations'][33][i][1:1+len(seq)] # the first and last embeding are padding adn every thing after the length of the sequence is padding
            sequenceEmbedings[seq] = trimmedRep.mean(0) # mean pooling for representation of entire protein sequence
    except KeyError as e:
        print('Error in one of the sequences so doing the batch one at a time to figure it out')
        # clear memory (this should be garbadge collected here anyways)
        sequenceEmbedings = {}
        for i, seq in enumerate(sequences):
           # try:
             #   _, _, batch_tokens = batch_converter(list(enumerate([seq])))
            #    with torch.no_grad():
              #      results = model(batch_tokens, repr_layers=[33], return_contacts=True)
               # trimmedRep = results['representations'][33][i][1:1+len(seq)]
               # sequenceEmbedings[seq] = trimmedRep.mean(0)
           # except KeyError as e:
            #    print(f"Skipping sequence {i} due to KeyError: {e}")
                sequenceEmbedings[seq] = None

    return pd.Series(sequenceEmbedings)



# Define chunk size
chunk_size = 20
N = 5

# Read and process data in chunks
trainingDfPath = os.path.join('data', 'training_data.csv')
savePath = os.path.join('data', 'embeddings2.csv')
npzFilePath = os.path.join('data', 'arrays.npz')

for i, dfCombinedSmall in enumerate(pd.read_csv(trainingDfPath, chunksize=chunk_size)):
    # Get embeddings for each group of N rows
    embeddingDf = dfCombinedSmall.groupby(dfCombinedSmall.index // N)['AminoAcidSeq'].apply(get_embeddings_batch).reset_index(level = 0, drop=True)

    # Reset and rename the index
    embeddingDf = embeddingDf.reset_index().rename(columns={'index':'AminoAcidSeq2'})
    embeddingDf = embeddingDf.rename(columns={'AminoAcidSeq2':'AminoAcidSeq', 'AminoAcidSeq':'embeddings'})

    # just add numt to the embeddings df
    merged_df = pd.merge(embeddingDf, dfCombinedSmall[['AminoAcidSeq', 'numt']], on='AminoAcidSeq', how='left')
    # y = merged_df['numt'].to_numpy() # for saving in np format for using with models
    # X = np.array([embedding.numpy() for embedding in merged_df['embeddings']])
    del embeddingDf, dfCombinedSmall

    # Save the results to a CSV file, appending if the file already exists
    if i == 0:
        print(f'incremental save {i}')
        # np.savez(, X=X, y=y) # NO EASY WAY TO APPEND
        merged_df.to_csv(savePath, mode='w', index=False)
    else:
        print(f'incremental save {i}')
        merged_df.to_csv(savePath, mode='a', header=False, index=False)

    # dont really need this as long as it is outputed because then its saved
    # with open(os.path.join('data', 'last_processed_chunk_index.txt'), 'w') as f:
    #     f.write(str(i)) #

