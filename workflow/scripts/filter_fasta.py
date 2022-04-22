from Bio import SeqIO

rs = [
    r
    for r in SeqIO.parse(snakemake.input[0], "fasta")
    if r.id in snakemake.params["filt"]
]

SeqIO.write(rs, snakemake.output[0], "fasta")
