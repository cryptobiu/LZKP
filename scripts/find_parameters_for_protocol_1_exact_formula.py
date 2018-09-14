#!/usr/bin/env python

"""find_parameters_for_protocol_1.py: Finds optimal parameters for the "Cut-and-Choose" protocol. Using exact formula."""

__author__ = 'Roee Sefi'
__credits__ = ['Ariel Nof', 'Carsten Baum']
__email__ = 'roee.sefi@gmail.com'

import operator as op

# Implementation from https://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python
def nCr(n, r):
    r = min(r, n - r)
    numer = reduce(op.mul, xrange(n, n - r, -1), 1)
    denom = reduce(op.mul, xrange(1, r + 1), 1)

    return numer // denom

print '{:<5}{:<5}{:<5}{:<25}{:<20}{:<25}{:<20}'.format('M', 'N', 'q', 'Optimal tau for 2^(-40)', 'M - tau (40)', 'Optimal tau for 2^(-80)', 'M - tau (80)')
for M in range(5, 200 + 1, 10):
	for N in [2, 4, 8, 16, 32, 64]:
		for q in [15, 31, 59, 61]:
			tau_80 = 0	# Finds best tau for 2^(-80) security
			tau_40 = 0  # Finds best tau for 2^(-40) security

			for tau in range(1, M):

				flag_80 = True
				flag_40 = True

				# c1 = The number of off-line incorrect squares generated by P
				# c2 = The number of on-line emulations where P cheats
				for c1 in range(M - tau + 1):
					if not flag_40:
						break
						
					for c2 in range(tau + 1):
						numerator = nCr(M - c1, tau)
						denominator = (nCr(M, tau) * (N ** c2) * (q ** (M - tau - c1 - c2)))
						
						# error = 1.0 * nCr(M - c1, tau) / (nCr(M, tau) * (N ** (M - tau - c1))) - Avoid floating point

						# if error > 2 ** (-80):
						if numerator * (2 ** 80) > denominator:
							flag_80 = False

						# if error > 2 ** (-40):
						if numerator * (2 ** 40) > denominator:
							flag_40 = False

							break

				if flag_80:
					tau_80 = tau

				if flag_40:
					tau_40 = tau

			print '{:<5}{:<5}{:<5}{:<25}{:<20}{:<25}{:<20}'.format(M, N, q, tau_40 if tau_40 > 0 else 'None', M - tau_40 if tau_40 > 0 else '', 
																	tau_80 if tau_80 > 0 else 'None', M - tau_80 if tau_80 > 0 else '')


