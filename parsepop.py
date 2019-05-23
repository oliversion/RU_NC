import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv

def get_ca(filename):
	data = list()
	with open(filename) as file:
		reader = csv.reader(file)
		for row in reader:
			data.append(list(map(int,row)))
	x_max = max(row[0] for row in data)
	y_max = max(row[1] for row in data)
	x_min = min(row[0] for row in data)
	y_min = min(row[1] for row in data)
	datarr = [[0]*(x_max-x_min+1) for i in range(y_max-y_min+1)]
	for row in data:
		datarr[row[1]-y_min][row[0]-x_min] = row[2]
	return datarr

def print_map(ca,mapname):
	plt.imshow(ca)
	plt.colorbar()
	plt.savefig('{}.png'.format(mapname))
	plt.clf()

ca_100x100 = get_ca('100x100data.csv')
print_map(ca_100x100, '100x100map')
ca_500x500 = get_ca('500x500data.csv')
print_map(ca_500x500, '500x500map')
