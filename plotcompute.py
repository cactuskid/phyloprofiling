import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit


df = pd.read_csv( 'observations.txt', names = ['size', 'time'] )
df.sort_values('size')

fig = plt.figure()
ax = plt.gca()
print(df)




newX = np.logspace(.5, 3, base=10)  # Makes a nice domain for the fitted curves.
                                   # Goes from 10^0 to 10^3
                                   # This avoids the sorting and the swarm of lines.

# Let's fit an exponential function.  
# This looks like a line on a lof-log plot.
def myExpFunc(x, a):
    return np.power(a, x)

popt, pcov = curve_fit(myExpFunc, df['size'] , df['time'])
print( popt)
plt.plot(newX, myExpFunc(newX, *popt), 'r-', 
         label="({0:.3f} ^ x)".format(*popt))

residuals = df['time']- myExpFunc(df['size'] , popt[0]) #, popt[1] )
ss_res = np.sum(residuals**2)
ss_tot = np.sum((df['time']-np.mean(df['time']))**2)


r_squared = 1 - (ss_res / ss_tot)
print (r_squared
)

ax.grid(b='on')
plt.plot([], [], ' ', label="R^2 :"+str(r_squared))
plt.legend(loc='lower right' )

ax.scatter(df['size'] , df['time'] , c='blue', alpha=0.5, edgecolors='none')
ax.set_yscale('log')
plt.show()

