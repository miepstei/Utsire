#!/usr/bin/python

from dcprogs.likelihood import Log10Likelihood,IdealG,MissedEventsG, QMatrix, ExactSurvivor
import numpy as np

print "Generating numbers for 1996 CHS test case"

print "Setting up parameters\n"


tres = 5e-5
tcrit = 1e-4
nopen = 2
k = 5

Q = QMatrix([ [ -3050, 50,           0,      3000,   0 ],
      [ 2/3,   -(500 + 2/3), 500,    0,      0 ],
      [ 0,     15000,        -19000, 4000,   0 ],
      [ 15,    0,            50,     -2065,  2000],
      [ 0,     0,            0,      10,     -10] ], nopen)

dcpOptions = {'nmax': 2, 'xtol': 1e-12, 'rtol': 1e-12, 'itermax': 100,
            'lower_bound': -1e6, 'upper_bound': 0}

bursts = [  [0.1, 0.2, 0.1],                  # 1st burst 
            [0.2],                            # 2nd burst
            [0.15, 0.16, 0.18, 0.05, 0.1] ]   # 3rd burst

t = [1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 1e-4]

print "*** Test 1: Generating likelihood ***"

likelihood = Log10Likelihood(bursts, nopen, tres, tcrit)
result = likelihood(Q)
print("Computation for missed events likelihood: {:.14f}".format(result*np.log(10)))

print "*" * 20 + "\n"

print "*** Test 2: IdealG ***"
idealg = IdealG(Q)

print idealg.af(t[9])
print idealg.fa(t[9])

print "*" * 20 + "\n"

print "*** Test 3: MissedEventsG ***"
missedeventsG = MissedEventsG(Q,tau=tres)

print missedeventsG.af(t[9])
print missedeventsG.fa(t[9])

print "*" * 20 + "\n"

print "*** Test 4: Equilibrium Occupancies"
print "\ta) Ideal occupancies"
print idealg.initial_occupancies
print idealg.final_occupancies

print "\tb) Missed Events occupancies"
print missedeventsG.initial_occupancies
print missedeventsG.final_occupancies

print "*" * 20 + "\n"

print "*** Test 5: CHS occupancies"
print missedeventsG.initial_CHS_occupancies(tcrit)
print missedeventsG.final_CHS_occupancies(tcrit)
print "*" * 20 + "\n"

print "*** Test 6: Exact Survivor Recursive Matrices"
survivor = ExactSurvivor(Q, tres);
print(survivor.recursion_af(0, 0, 0))
print "*" * 20 + "\n"



