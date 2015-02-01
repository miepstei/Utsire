from dcprogs.likelihood import Log10Likelihood
import numpy as np

bursts = [  [0.1, 0.2, 0.1],                  # 1st burst 
            [0.2],                            # 2nd burst
            [0.15, 0.16, 0.18, 0.05, 0.1] ]   # 3rd burst
""" List of bursts.  

    Each burst is a list of observed open and shut intervals. 
    There should always be an odd number of intervals, since bursts end in a shut states.
"""

likelihood = Log10Likelihood(bursts, nopen=1, tau=2.5e-5, tcritical=1e-4)
matrix = [[ -100, 100 ], 
          [ 5e6 * 1e-6 , -5e6 * 1e-6]] 

print(matrix)
result = likelihood(matrix)
print("Computation for Julia missed events likelihood: {:.14f}".format(result*np.log(10)))

likelihood = Log10Likelihood(bursts, nopen=1, tau=1e-4, tcritical=1e-2)
result = likelihood(matrix)
print("Computation for Julia missed events likelihood: tres=1e-4 {:.14f}".format(result*np.log(10)))
likelihood = Log10Likelihood(bursts, nopen=1, tau=1e-5, tcritical=1e-2)
result = likelihood(matrix)
print("Computation for Julia missed events likelihood tau=1e-5: {:.14f}".format(result*np.log(10)))
likelihood = Log10Likelihood(bursts, nopen=1, tau=1e-4, tcritical=1e-3)
result = likelihood(matrix)
print("Computation for Julia missed events likelihood tres=1e-3: {:.14f}".format(result*np.log(10)))
likelihood = Log10Likelihood(bursts, nopen=1, tau=1e-4, tcritical=-1e-2)
result = likelihood(matrix)
print("Computation for Julia missed events likelihood tcrit=-1e-2 (noChs): {:.14f}".format(result*np.log(10)))

bursts.append([0.10,0.18,0.02,0.1,0.03])
likelihood = Log10Likelihood(bursts, nopen=1, tau=1e-4, tcritical=1e-2)
result = likelihood(matrix)
print("Computation for additional burst [0.10,0.18,0.02,0.1,0.03]: {:.14f}".format(result*np.log(10)))

