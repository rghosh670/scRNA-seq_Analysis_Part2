import pandas as pd

f = open('wardClusters/cluster0.txt', 'r')
genes = f.read().splitlines()
f.close()

eset = pd.read_csv('data/rawReadsWithCellNames.csv',index_col=[0], usecols=genes)

print(eset)