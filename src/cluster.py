import pandas
import hdbscan

df = pandas.read_csv("data/Ecoli_K12_MG1655.3UTR.mRNA.seq.cdhit.1.csv")
print(df.head())