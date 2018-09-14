#!/usr/bin/env python

import xlwt
import sys


# From https://stackoverflow.com/questions/15389768/standard-deviation-of-a-list
def mean(data):
    """Return the sample arithmetic mean of data."""
    n = len(data)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return sum(data)/float(n) # in Python 2 use sum(data)/float(n)

def _ss(data):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def stddev(data, ddof=0):
    """Calculates the population standard deviation
    by default; specify ddof=1 to compute the sample
    standard deviation."""
    n = len(data)
    # if n < 2:
        # raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/(n-ddof)
    return pvar**0.5

if len(sys.argv) != 4:
	print 'usage: ./parseMeasurements.py <execServers-file> <measurments-file> <results-file>'
	sys.exit(1)

config = open(sys.argv[1], 'r').readlines()
q = [c.replace('\n', '') for c in config if c.startswith('q=(')]
n = [c.replace('\n', '') for c in config if c.startswith('n=(')]
m = [c.replace('\n', '') for c in config if c.startswith('m=(')]
N = [c.replace('\n', '') for c in config if c.startswith('N=(')]
M = [c.replace('\n', '') for c in config if c.startswith('M=(')]
tau = [c.replace('\n', '') for c in config if c.startswith('tau=(')]
num_trials = int([c.replace('\n', '') for c in config if c.startswith('NUM_TRIALS=')][0].split('=')[1])

measurements = open(sys.argv[2], 'r').readlines()

commit_measurements = measurements[0].replace('\n', '')
measurements = [measurement.replace('\n', '') for measurement in measurements]

results = open(sys.argv[3], 'r').readlines()

commit_results = results[0].replace('\n', '')
results = [result.replace('\n', '') for result in results]

assert (commit_measurements == commit_results)
print '{}'.format(commit_measurements)
print 'Execution time: {}'.format(measurements[1])

m_id = 2
r_id = 2

print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<70}{:<10}{:<10}'.format('q', 'n', 'm', 'N', 'M', 'tau', 'MT?', 'times','AVG', 'STD')
for i_protocol, protocol in enumerate(['1', '2']):
	print measurements[m_id]
	assert measurements[m_id] == 'protocol {}'.format(protocol)
	m_id += 1

	for i in range(len(q[i_protocol].split('=')[1][1:-1].replace('"', '').split())):
		qq = q[i_protocol].split('=')[1][1:-1].replace('"', '').split()[i]
		nn = n[i_protocol].split('=')[1][1:-1].replace('"', '').split()[i]
		mm = m[i_protocol].split('=')[1][1:-1].replace('"', '').split()[i]

		NN = N[i_protocol].split('=')[1][1:-1].replace('"', '').split()[i]
		MM = M[i_protocol].split('=')[1][1:-1].replace('"', '').split()[i]

		for _ in range(num_trials):
			assert(results[r_id].split(',')[0] == '1')
			r_id += 1

		try:
			TT = tau[i_protocol].split('=')[1][1:-1].replace('"', '').split()[i]
		except:
			TT = ''

		for MMM, TTT in zip(MM.split(':'), TT.split(':')):
			mmm = map(int, [mmmm.split(',')[0] for mmmm in measurements[m_id].split()])
			print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<70}{:<10.4f}{:<10.4f}'.format(qq, nn, mm, NN, MMM, TTT, '', measurements[m_id], mean(mmm), stddev(mmm))
			# print qq, nn, mm, NN, MMM, measurements[m_id]
			m_id += 1
			print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<70}{:<10.4f}{:<10.4f}'.format(qq, nn, mm, NN, MMM, TTT, 'X', measurements[m_id], mean(mmm), stddev(mmm))
			# print qq, nn, mm, NN, MMM, measurements[m_id]
			m_id += 1