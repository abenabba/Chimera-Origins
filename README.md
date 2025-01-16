# Chimera-Origins
Analysis and visualization of mutagenized library sequence alignments

This repository documents the bioinformatics pipeline I used to analyze a mutagenized fluorescent protein library. The library was created using DNA shuffling to fragment parental variants, followed by random re-annealing through primerless PCR. The resulting mutagenized library was transformed in *E. coli*, flow-sorted based on brightness. The sequencing data from this sort was later used as a training set for a ProtGPT2 language model. The full library preparation process is depicted in the graphic below:

<img width="934" alt="image" src="https://github.com/user-attachments/assets/f875a5ae-e5f8-4687-8eda-1678d592339b" />

After long-read sequencing of the libraries, I aligned the shuffled sequences with their parental variants to identify which sections of the chimera were derived from each parent. Using BLAST, I selected variants with the highest bit score (indicates stronger similarity) and lowest e value (lower E-value signifies that the alignment is less likely to be a random occurrence and is therefore more likely to be biologically meaningful). From these top hits, I used [EMBOSS NEEDLE](https://www.ebi.ac.uk/jdispatcher/psa/emboss_needle) to determine precisely where the chimera aligns with the parents. A summary of these findings is included [in this doc](https://docs.google.com/document/d/16SDklMaBfDLj5UsnRXUyHSzPwHZA0g1PTXB5TC_d980/edit?usp=sharing).

After determining which parents align to which sections of the chimera, I used (alphafold)[https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=kOblAo-xetgx] to predict its structure. I then uploded the structures on to Pymol and colored the different parental sections of the chimera.
In this repository, I've included the fasta files for one of the variants mentioned in the doc. I've also included the Fasta file of the library sequencing for both the parents and shuffled library. Below is the annotated pymol image I made of this variant:
![final structures for figure2](https://github.com/user-attachments/assets/f75c4097-894c-4703-9fdc-50638111f284)


# Reproducing the Analysis
1. Set Up the Environment
Initialize Conda and activate the find_chimera environment:

```
conda init
conda activate find_chimera
```

Install BLAST:

```
conda install -n find_chimera -c bioconda blast -y
```

2. Prepare the BLAST Databases
Create BLAST databases for the parental library and shuffled library sequences:

```
# Create a BLAST database for the parental library

conda run -n find_chimera makeblastdb -in caa.fasta -dbtype nucl -out caa.db

# Create a BLAST database for the shuffled sequences

conda run -n find_chimera makeblastdb -in BC2_after_consesnus_R.fasta -dbtype nucl -out BC2_after_consesnus_R.db
```

3. Perform the BLAST Alignments
Align a shuffled sequence against the parental library database:

```
conda run -n find_chimera blastn -query BC2_after_consesnus_R.fasta -db caa.db -evalue 1e-6 -outfmt 6 -out BC2_alignment_results.tsv
```

Align the parental library sequences against the shuffled database:

```
conda run -n find_chimera blastn -query All_parents_codon_labeled_nobrac.fasta -db caa.db -evalue 1e-6 -outfmt 6 -out BC2_alignment_results2.tsv
```

4. Changing the Alignment Direction
To align the entire parental library (caa.fasta) with the shuffled library (BC2_after_consesnus_R.fasta), you can swap the roles of the query and database files. Here’s the adjusted command:
```
# Align the parental library against the shuffled library
conda run -n find_chimera blastn -query caa.fasta -db BC2_after_consesnus_R.db -evalue 1e-6 -outfmt 6 -out Parental_vs_Shuffled_alignment.tsv
This ensures that each sequence from the parental library is aligned to all sequences in the shuffled library.
```
5. Output Details
```
The -outfmt 6 flag specifies the tabular format of the alignment output, including details like bit scores, E-values, and alignment length
```
