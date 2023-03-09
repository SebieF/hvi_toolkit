import ankh
import h5py
import torch

from Bio import SeqIO
from tqdm import tqdm

fasta_path = "../merged_dataset_complete.fasta"
output_path = "reduced_embeddings_file_ankh.h5"

protein_sequences = {seq.id: list(seq.seq) for seq in sorted(list(SeqIO.parse(fasta_path, "fasta")),
                                                             key=lambda seq: len(seq.seq),
                                                             reverse=True)}
print(f"Embedding {len(protein_sequences)} protein sequences!")

model, tokenizer = ankh.load_large_model()
model.eval()
device = "cuda:0" if torch.cuda.is_available() else "cpu"
model.to(device)

ankh_large_dim = 1536

with torch.no_grad():
    embeddings_dict = {}
    for seq_id, sequence in tqdm(protein_sequences.items()):
        outputs = tokenizer.batch_encode_plus([sequence],
                                              add_special_tokens=True,
                                              padding=True,
                                              is_split_into_words=True,
                                              return_tensors="pt")
        embeddings = model(input_ids=outputs['input_ids'].to(device),
                           attention_mask=outputs['attention_mask'].to(device))
        embedding = embeddings[0].detach().cpu()
        per_protein_embedding = torch.mean(embedding.squeeze(), dim=0)
        embeddings_dict[seq_id] = per_protein_embedding

with h5py.File(output_path, "w") as output_file:
    idx = 0
    for seq_id, embedding in embeddings_dict.items():
        output_file.create_dataset(str(idx), data=embedding, compression="gzip", chunks=True,
                                   maxshape=ankh_large_dim)
        output_file[str(idx)].attrs["original_id"] = seq_id
        idx += 1

# Verify h5 file
created_file = h5py.File(output_path, 'r')

for idx, embedding in created_file.items():
    original_sequence_id = created_file[idx].attrs["original_id"]
    assert embedding.shape[0] == ankh_large_dim, "New dimension is not correct"
    assert embedding[0] == embeddings_dict[original_sequence_id][0]

assert len(created_file.keys()) == len(protein_sequences.keys())
created_file.close()
