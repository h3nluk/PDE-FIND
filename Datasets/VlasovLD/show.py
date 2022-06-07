import numpy as np
import matplotlib.pyplot as plt

f_e = np.load('f_e.npy')
x = np.load('x.npy')
v = np.load('v.npy')
E = np.load('E.npy')
t = np.load('t.npy')


print('fe shape:', np.shape(f_e))
print('x shape:', np.shape(x))
print('v shape:', np.shape(v))
print('E shape:', np.shape(E))
print('t shape:', np.shape(t))

# ~ print(x)
# ~ print(v)
# ~ print(t)



