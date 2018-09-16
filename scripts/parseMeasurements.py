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
config = [c.split('#')[0] for c in config]
q = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('q=(')]
n = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('n=(')]
m = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('m=(')]
N = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('N=(')]
M = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('M=(')]
tau = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('tau=(')]
num_trials = int([c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('NUM_TRIALS=')][0].split('=')[1])

measurements = open(sys.argv[2], 'r').readlines()

commit_measurements = measurements[0].replace('\n', '').replace('\t', '').strip()
measurements = [measurement.replace('\n', '').replace('\t', '').strip() for measurement in measurements]

results = open(sys.argv[3], 'r').readlines()

commit_results = results[0].replace('\n', '').replace('\t', '').strip()
results = [result.replace('\n', '').replace('\t', '').strip() for result in results]

assert (commit_measurements == commit_results)
print '{}'.format(commit_measurements)
print 'Execution time: {}'.format(measurements[1])

m_id = 2
r_id = 2

print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<55}{:<10}{:<10}{:<55}{:<10}{:<10}'.format('q', 'n', 'm', 'N', 'M', 'tau', 'MT?', 'total times', 'AVG', 'STD', 'computation times', 'AVG', 'STD')
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
		print MM
		for _ in range(num_trials):
			assert(results[r_id].split(',')[-1] == '1')
			r_id += 1

		try:
			TT = tau[i_protocol].split('=')[1][1:-1].replace('"', '').split()[i]
		except:
			TT = ':'.join(['' for _ in range(len(MM.split(':')))])

		for MMM, TTT in zip(MM.split(':'), TT.split(':')):
			mmm = map(int, [mmmm.split(',')[0] for mmmm in measurements[m_id].split()])
			ccc = map(int, [cccc.split(',')[1] for cccc in measurements[m_id].split()])
			print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<55}{:<10.4f}{:<10.4f}{:<55}{:<10.4f}{:<10.4f}'.format(qq, nn, mm, NN, MMM, TTT, '',
			 	  ','.join(map(str, mmm)), mean(mmm), stddev(mmm), ','.join(map(str, ccc)), mean(ccc), stddev(ccc))
			# print qq, nn, mm, NN, MMM, measurements[m_id]
			m_id += 1
			mmm = map(int, [mmmm.split(',')[0] for mmmm in measurements[m_id].split()])
			ccc = map(int, [cccc.split(',')[1] for cccc in measurements[m_id].split()])
			print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<55}{:<10.4f}{:<10.4f}{:<55}{:<10.4f}{:<10.4f}'.format(qq, nn, mm, NN, MMM, TTT, 'X',
			      ','.join(map(str, mmm)), mean(mmm), stddev(mmm), ','.join(map(str, ccc)), mean(ccc), stddev(ccc))
			# print qq, nn, mm, NN, MMM, measurements[m_id]
			m_id += 1