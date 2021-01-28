# this script convert the full dataset mitdb (data and annotatiosn) to text files

from os import listdir, mkdir, system
from os.path import isfile, isdir, join, exists

#dir = 'ediagnostic/wfdb/'#'mitdb/'
dir = 'mitdb/'
#Create folder
dir_out = dir + 'csv2/'
if not exists(dir_out):
	mkdir(dir_out)

records = [f for f in listdir(dir) if isfile(join(dir, f)) if(f.find('.dat') != -1)]
#print records

for r in records:

	command = 'rdsamp -r ' + dir + r[:-4] + ' -c -H -f 0 -v >' + dir_out + r[:-4] + '.csv'
	print(command)
	system(command)

	command_annotations = 'rdann -r ' + dir + r[:-4] +' -f 0 -a atr -v > ' + dir_out + r[:-4] + '.ann'
	print(command_annotations)
	system(command_annotations)
	
	write_annotations = 'wrann -r ' + dir + r[:-4] + ' -a qrs.txt ' + dir_out + r[:-4] + '.atr'
	print(write_annotations)
	system(write_annotations)
