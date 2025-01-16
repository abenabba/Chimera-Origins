# Chimera-Origins
Analysis and visualization of mutagenized library sequence alignments

This repository documents the bioinformatics pipeline I used to analyze a mutagenized fluorescent protein library. The library was created using DNA shuffling to fragment parental variants, followed by random re-annealing through primerless PCR. The resulting mutagenized library was transformed in *E. coli*, flow-sorted based on brightness. The sequencing data from this sort was later used as a training set for a ProtGPT2 language model. The full library preparation process is depicted in the graphic below:

<img width="934" alt="image" src="https://github.com/user-attachments/assets/f875a5ae-e5f8-4687-8eda-1678d592339b" />

After long-read sequencing of the libraries, I aligned the shuffled sequences with their parental variants to identify which sections of the chimera were derived from each parent. Using BLAST, I selected variants with the highest bit score (indicates stronger similarity) and lowest e value (lower E-value signifies that the alignment is less likely to be a random occurrence and is therefore more likely to be biologically meaningful). From these top hits, I used [EMBOSS NEEDLE](https://www.ebi.ac.uk/jdispatcher/psa/emboss_needle) to determine precisely where the chimera aligns with the parents. A summary of these findings is included [in this doc](https://docs.google.com/document/d/16SDklMaBfDLj5UsnRXUyHSzPwHZA0g1PTXB5TC_d980/edit?usp=sharing).

After determining which parents align to which sections of the chimera, I used (alphafold)[https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=kOblAo-xetgx] to predict its structure. I then uploded the structures on to Pymol and colored the different parental sections of the chimera.
In this repository, I've included the fasta files for one of the variants mentioned in the doc. I've also included the Fasta file of the library sequencing for both the parents and shuffled library. Below is the annotated pymol image I made of this variant:
![final structures for figure2](https://github.com/user-attachments/assets/f75c4097-894c-4703-9fdc-50638111f284)

