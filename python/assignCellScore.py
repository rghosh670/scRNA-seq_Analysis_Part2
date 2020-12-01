import pandas as pd

cellNames = pd.read_csv('data/t6.csv',index_col=[0], usecols=[0]).index.tolist()
cellScore = pd.DataFrame(index = cellNames)

cellSum = pd.read_csv('data/t6sum.csv')['sum'].tolist()

for j in range(15):
    f = open('leidenClusters/cluster' + str(j) + '.txt', 'r')
    genes = f.read().splitlines()
    genes.append('cell')
    f.close()

    eset = pd.read_csv('data/t6.csv',index_col=[0], usecols=genes)

    sum = eset.sum(axis = 1).tolist()
    scores = []

    for index, value in enumerate(sum):
        scores.append(value/cellSum[index])

    cellScore['cluster' + str(j)] = scores

cellScore.to_csv('leidenClusters/cellBinLeidenScore.csv')
# exec(open('python/sendEmail.py').read())