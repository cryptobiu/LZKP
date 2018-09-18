#!/usr/bin/env python

from xlutils.copy import copy
import xlrd
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

try:
	rb = xlrd.open_workbook('stats.xls', formatting_info=True)
	book = copy(rb)
except:
	book = xlwt.Workbook()

sheet = book.add_sheet("v1.7", cell_overwrite_ok=True)

config = open(sys.argv[1], 'r').readlines()
config = [c.split('#')[0] for c in config]
q = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('q=(')]
n = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('n=(')]
m = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('m=(')]
N = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('N=(')]
M = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('M=(')]
tau = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('tau=(')]
X = [c.replace('\n', '').replace('\t', '').strip() for c in config if c.startswith('X=(')]
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

row_id = 0
col_id = 3
for c, v in enumerate(['q', 'n', 'm', 'N', 'M', 'tau', 'MT?', 'total times', 'AVG', 'STD', 'computation times', 'AVG', 'STD', 'equation #1 times', 'AVG', 'STD']):
	if v == 'total times' or v == 'computation times' or v == 'equation #1 times':
		sheet.write_merge(row_id, row_id, col_id, col_id + num_trials - 1, v)
		col_id += num_trials
	else:
		sheet.write(row_id, col_id, v)	
		col_id += 1
row_id += 1


print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<75}{:<12}{:<12}{:<75}{:<12}{:<12}{:<75}{:<12}{:<12}{:<75}{:<12}{:<12}'.format('q', 'n', 'm', 'N', 'M', 'tau', 'MT?', 'total times', 'AVG', 'STD', 'computation times', 'AVG', 'STD', 'equation #1 times', 'AVG', 'STD', 'cut-and-choose times', 'AVG', 'STD')
for i_protocol, protocol in enumerate(['1', '2']):
	print measurements[m_id]
	sheet.write(row_id, 0, measurements[m_id])
	row_id += 1

	assert measurements[m_id] == 'protocol {}'.format(protocol)
	m_id += 1

	for i in range(len(q[i_protocol].split('=')[1][1:-1].replace('"', '').split())):
		qq = q[i_protocol].split('=')[1][1:-1].replace('"', '').split()[i]
		nn = n[i_protocol].split('=')[1][1:-1].replace('"', '').split()[i]
		mm = m[i_protocol].split('=')[1][1:-1].replace('"', '').split()[i]

		for j in range(len(N[i_protocol].split('=')[1][1:-1].replace('"', '').split())):
			NN = N[i_protocol].split('=')[1][1:-1].replace('"', '').split()[j]
			MM = M[i_protocol].split('=')[1][1:-1].replace('"', '').split()[j]

			for _ in range(num_trials):
				assert(results[r_id].split(',')[-1] == '1')
				r_id += 1

			try:
				TT = tau[i_protocol].split('=')[1][1:-1].replace('"', '').split()[j]
			except:
				TT = ':'.join(['' for _ in range(len(MM.split(':')))])

			for MMM, TTT in zip(MM.split(':'), TT.split(':')):
				mmm = map(int, [mmmm.split(',')[0] for mmmm in measurements[m_id].split()])
				ccc = map(int, [cccc.split(',')[1] for cccc in measurements[m_id].split()])
				eq1 = map(int, [eq11.split(',')[2] for eq11 in measurements[m_id].split()])

				if i_protocol == 0:
					cac = map(int, [cacc.split(',')[3] for cacc in measurements[m_id].split()])

				if i_protocol == 0:
					print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<75}{:<12.4f}{:<12.4f}{:<75}{:<12.4f}{:<12.4f}{:<75}{:<12.4f}{:<12.4f}{:<75}{:<12.4f}{:<12.4f}'.format(qq, nn, mm, NN, MMM, TTT, '',
				 	  	','.join(map(str, mmm)), mean(mmm), stddev(mmm), ','.join(map(str, ccc)), mean(ccc), stddev(ccc), ','.join(map(str, eq1)), mean(eq1), stddev(eq1), ','.join(map(str, cac)), mean(cac), stddev(cac))
				else:
					print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<75}{:<12.4f}{:<12.4f}{:<75}{:<12.4f}{:<12.4f}{:<75}{:<12.4f}{:<12.4f}'.format(qq, nn, mm, NN, MMM, TTT, '',
				 	  	','.join(map(str, mmm)), mean(mmm), stddev(mmm), ','.join(map(str, ccc)), mean(ccc), stddev(ccc), ','.join(map(str, eq1)), mean(eq1), stddev(eq1))

				m_id += 1

				col_id = 3

				if i_protocol == 0:
					for c, (v1, v2) in enumerate(zip(['q', 'n', 'm', 'N', 'M', 'tau', 'MT?', 'total times', 'AVG', 'STD', 'computation times', 'AVG', 'STD', 'equation #1 times', 'AVG', 'STD', 'cut-and-choose times', 'AVG', 'STD'],
											         [qq, nn, mm, NN, MMM, TTT, '', ','.join(map(str, mmm)), mean(mmm), stddev(mmm), ','.join(map(str, ccc)), mean(ccc), stddev(ccc), ','.join(map(str, eq1)), mean(eq1), stddev(eq1), ','.join(map(str, cac)), mean(cac), stddev(cac)])):
						if v1 == 'total times' or v1 == 'computation times' or v == 'equation #1 times' or v == 'cut-and-choose times':
							for vv2 in v2.split(','):
								sheet.write(row_id, col_id, vv2)
								col_id += 1
						else:
							sheet.write(row_id, col_id, v2)
							col_id += 1
				else:
					for c, (v1, v2) in enumerate(zip(['q', 'n', 'm', 'N', 'M', 'tau', 'MT?', 'total times', 'AVG', 'STD', 'computation times', 'AVG', 'STD', 'equation #1 times', 'AVG', 'STD'],
											         [qq, nn, mm, NN, MMM, TTT, '', ','.join(map(str, mmm)), mean(mmm), stddev(mmm), ','.join(map(str, ccc)), mean(ccc), stddev(ccc), ','.join(map(str, eq1)), mean(eq1), stddev(eq1)])):
						if v1 == 'total times' or v1 == 'computation times' or v == 'equation #1 times':
							for vv2 in v2.split(','):
								sheet.write(row_id, col_id, vv2)
								col_id += 1
						else:
							sheet.write(row_id, col_id, v2)
							col_id += 1

				row_id += 1

	qq = q[i_protocol].split('=')[1][1:-1].replace('"', '').split()[-1]
	nn = n[i_protocol].split('=')[1][1:-1].replace('"', '').split()[-1]
	mm = m[i_protocol].split('=')[1][1:-1].replace('"', '').split()[-1]
	NN = N[i_protocol].split('=')[1][1:-1].replace('"', '').split()[-1]
	MM = M[i_protocol].split('=')[1][1:-1].replace('"', '').split()[-1]

	try:
		TT = tau[i_protocol].split('=')[1][1:-1].replace('"', '').split()[-1]
	except:
		TT = ':'.join(['' for _ in range(len(MM.split(':')))])

	XX = X[i_protocol].split('=')[1][1:-1].replace('"', '').split()

	for XXX in XX:
		for MMM, TTT in zip(MM.split(':'), TT.split(':')):
			for _ in range(num_trials):
				assert(results[r_id].split(',')[-1] == '1')
				r_id += 1

			mmm = map(int, [mmmm.split(',')[0] for mmmm in measurements[m_id].split()])
			ccc = map(int, [cccc.split(',')[1] for cccc in measurements[m_id].split()])
			eq1 = map(int, [eq11.split(',')[2] for eq11 in measurements[m_id].split()])

			if i_protocol == 0:
				cac = map(int, [cacc.split(',')[3] for cacc in measurements[m_id].split()])

			if i_protocol == 0:
				print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<75}{:<12.4f}{:<12.4f}{:<75}{:<12.4f}{:<12.4f}{:<75}{:<12.4f}{:<12.4f}{:<75}{:<12.4f}{:<12.4f}'.format(qq, nn, mm, NN, MMM, TTT, XXX,
			 	  	','.join(map(str, mmm)), mean(mmm), stddev(mmm), ','.join(map(str, ccc)), mean(ccc), stddev(ccc), ','.join(map(str, eq1)), mean(eq1), stddev(eq1), ','.join(map(str, cac)), mean(cac), stddev(cac))
			else:
				print '\t\t{:<5}{:<5}{:<7}{:<5}{:<5}{:<5}{:<5}{:<75}{:<12.4f}{:<12.4f}{:<75}{:<12.4f}{:<12.4f}{:<75}{:<12.4f}{:<12.4f}'.format(qq, nn, mm, NN, MMM, TTT, XXX,
			 	  	','.join(map(str, mmm)), mean(mmm), stddev(mmm), ','.join(map(str, ccc)), mean(ccc), stddev(ccc), ','.join(map(str, eq1)), mean(eq1), stddev(eq1))

			m_id += 1

			col_id = 3
			if i_protocol == 0:
				for c, (v1, v2) in enumerate(zip(['q', 'n', 'm', 'N', 'M', 'tau', 'MT?', 'total times', 'AVG', 'STD', 'computation times', 'AVG', 'STD', 'equation #1 times', 'AVG', 'STD'],
											     [qq, nn, mm, NN, MMM, TTT, XXX, ','.join(map(str, mmm)), mean(mmm), stddev(mmm), ','.join(map(str, ccc)), mean(ccc), stddev(ccc), ','.join(map(str, eq1)), mean(eq1), stddev(eq1), ','.join(map(str, cac)), mean(cac), stddev(cac)])):
					if v1 == 'total times' or v1 == 'computation times' or v == 'equation #1 times' or v == 'cut-and-choose times':
						for vv2 in v2.split(','):
							sheet.write(row_id, col_id, vv2)
							col_id += 1
					else:
						sheet.write(row_id, col_id, v2)
						col_id += 1
			else:
				for c, (v1, v2) in enumerate(zip(['q', 'n', 'm', 'N', 'M', 'tau', 'MT?', 'total times', 'AVG', 'STD', 'computation times', 'AVG', 'STD', 'equation #1 times', 'AVG', 'STD'],
											     [qq, nn, mm, NN, MMM, TTT, XXX, ','.join(map(str, mmm)), mean(mmm), stddev(mmm), ','.join(map(str, ccc)), mean(ccc), stddev(ccc), ','.join(map(str, eq1)), mean(eq1), stddev(eq1)])):
					if v1 == 'total times' or v1 == 'computation times' or v == 'equation #1 times':
						for vv2 in v2.split(','):
							sheet.write(row_id, col_id, vv2)
							col_id += 1
					else:
						sheet.write(row_id, col_id, v2)
						col_id += 1

			row_id += 1

book.save('stats.xls')