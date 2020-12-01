import scipy
from scipy.cluster.hierarchy import fcluster
from scipy import stats
import scipy.cluster.hierarchy as hac
import statistics
from statistics import mode
from math import floor
import random
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import matplotlib
from cellNames import getCellName
import datetime
from numpy import genfromtxt

dist = genfromtxt('data/betterCondensedDistance.csv', delimiter=',')
dist = dist.tolist()
print(dist)
print(len(dist))

Z = hac.linkage(dist,  method='ward')
np.savetxt('data/betterWardLinkage.csv', Z, delimiter = ',')
exec(open('python/sendEmail.py').read())