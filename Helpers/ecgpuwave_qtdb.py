# this script detect P, QRS, T onset/offset waves from MITDB data and save 
# in text files

from os import listdir, mkdir, system
from os.path import isfile, isdir, join, exists

dir = 'qtdb/'
#Create folder
dir_out = dir + 'qrs3/'
if not exists(dir_out):
	mkdir(dir_out)

records = [f for f in listdir(dir) if isfile(join(dir, f)) if(f.find('.dat') != -1)]
#print records

for r in records:
	command = 'ecgpuwave -r ' + dir +r[:-4] + ' -i atr -a qrs'
	print(command)
	system(command)

	command_annotations = 'rdann -r ' + dir +r[:-4] +' -f 0 -a qrs -v >' + dir_out + r[:-5] + '.txt'
	print(command_annotations)
	system(command_annotations)
