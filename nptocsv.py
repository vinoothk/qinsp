import numpy as np
import pandas as pd
# data = np.load('A.npy')
# np.savetxt('A.txt',data, delimiter = ',')
# print(data[0][0])
# print(data.shape,type(data),len(data))


df = pd.read_csv('B.dat', sep=",", header=None)
print(df.info)
