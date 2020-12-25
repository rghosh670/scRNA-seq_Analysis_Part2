import pandas as pd
import importlib
from cellNames import getCellName

NUM_CELLS = 37148

cells = []
for i in range(NUM_CELLS):
    cells.append(getCellName(i))

result = pd.factorize(cells)[0].tolist()
result = [int(i) + 1 for i in result]

with open('data/factorizedCellList', 'w') as f: # write out results to text file
    for item in result:
        f.write("%s\n" % item)

f.close()

