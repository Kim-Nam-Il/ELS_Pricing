import numpy as np
import matplotlib.pyplot as plt
N = 180
S = np.zeros([N,1])
S[0] = 100
vol = 0.079
r = 0.0165
T = 1
dt = T/N
t = np.linspace(0,T,N)
z = np.random.normal(0,1,N)
for i in range(N-1):
    S [i+1, 0] = S [i, 0] * np.exp((r-0.5*vol**2)*dt+vol*z[i]*np.sqrt(dt))
plt.plot(t,S[:,0],'ko-')
plt.xlabel('Time')
plt.ylabel('Stock Process')
plt.show()
