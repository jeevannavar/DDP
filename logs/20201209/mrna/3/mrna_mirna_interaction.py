#### ------------------------------------------------------------------------
#### This file is for selecting interactions and making them into a csv file
#### ------------------------------------------------------------------------

from csv import reader
import pandas as pd

# listing the miRNAs in the dataset
with open("BRCA/BRCA_mirna_preproc.csv") as f:
    line = f.readline().strip()

mirnas = line.split(",")[1:]

# listing the mRNAs in the dataset
with open("BRCA/BRCA_mrna_preproc.csv") as f:
    line = f.readline().strip()

mrnas = line.split(",")[1:]
converter = {mrna.split("|")[0]: mrna for mrna in mrnas}
mrnas = [mrna.split("|")[0] for mrna in mrnas]

# Using the mir2gene database to select the edges of the miRNA-mRNA network
mir2gene = {}
with open("Interaction_Network/mir2gene.csv") as f:
    file = reader(f)
    header = next(file)
    i=0
    for row in file:
        i+=1
        if row[0] in mirnas and row[2] in mrnas:
            if row[0] in mir2gene:
                if row[2] not in mir2gene[row[0]]:
                    mir2gene[row[0]].append(row[2])
            else:
                mir2gene[row[0]] = [row[2]]

#print(mir2gene)
print("Number of edges present in mir2gene network = ", i)
print("Number of edges selected =", sum([len(x) for x in mir2gene.values()]))

# Generating a csv file of the selected interactions
miRNAs = pd.read_csv('BRCA/BRCA_mirna_preproc.csv', index_col=0)
mRNAs = pd.read_csv('BRCA/BRCA_mrna_preproc.csv', index_col=0)
miRNA_mRNA = pd.DataFrame(index=miRNAs.index)

for miRNA in mir2gene.keys():
    for mRNA in mir2gene[miRNA]:
        interaction = miRNA+"_"+mRNA
        miRNA_mRNA[interaction] = miRNAs[miRNA]*mRNAs[converter[mRNA]]

miRNA_mRNA.to_csv("BRCA/BRCA_mrna_mirna_interaction.csv", index=True)