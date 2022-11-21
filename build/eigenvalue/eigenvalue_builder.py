#! /usr/bin/env python

import subprocess
import numpy
import os

partition_info=['CMT',16] # = [partition,ncores]
# partition_info=['debug',16] # = [partition,ncores]
time_str='7-00:00:00'
project_name=os.getcwd().split('/')[-3]
myemail=os.environ["MYEMAIL"]

#Log the current submission
logstr='''eigenvalue:

***write me***

Codebase: 
Cluster: ganymede
gittag: '''

#Do the git-ing
cmd='git commit -a -m "Commit before run eigenvalue" > /dev/null'
print(cmd)
subprocess.call(cmd,shell=True)

cmd='gittag.py > temp.temp'
subprocess.call(cmd,shell=True)
logstr+=open('temp.temp','r').readline()
subprocess.call('rm temp.temp',shell=True)

open('eigenvalue.log','w').write(logstr)

#Setup the versionmap and qsub files
vmap_file=open('versionmap.dat','w')
vmap_file.write('vnum\tL\n')

task_file=open('eigenvalue.task','w')
template_file='eigenvalue.template'
template_contents=open(template_file,'r').read()

vnum=0
deltas = [0.01,0.02,0.03,0.05,0.06,0.08,0.1,0.12,0.13,0.15,0.16,0.18,0.2,0.21,0.23,0.25]
for L in xrange(16):
	qsub_file=template_file.replace('.template','_'+str(vnum)+'.qsub')
	fout=open(qsub_file,'w')

	contents=template_contents.replace('###',str(vnum))
        contents=contents.replace('*project*',project_name)
	contents=contents.replace('*111*',str(deltas[L]))
	out_file_base = 'data_'+str(L)+'_*lll*.out'
	contents=contents.replace('*ppp*',out_file_base.replace('*lll*','julia'))
	vmap_file.write(str(vnum)+'\t'+str(L)+'\n')
	task_file.write('bash eigenvalue_'+str(vnum)+'.qsub\n')
	fout.write(contents)
	fout.close()
	
	vnum+=1
	

n_nodes=int(numpy.ceil(float(vnum)/partition_info[1]))

# Pad to an even multiple of cores per node
n_cores=n_nodes*partition_info[1]
for j in range(vnum,n_cores):
        task_file.write('echo "Fake run"\n')

# Finally output sbatch file
contents=open('eigenvalue.sbatch.template','r').read()
contents=contents.replace('*nnn*',str(n_cores)) # The total number of processors
contents=contents.replace('*NNN*',str(n_nodes)) # The total number of nodes
contents=contents.replace('*ttt*',time_str) # The wall clock time per processor
contents=contents.replace('*partition*',partition_info[0]) # Partition
contents=contents.replace('*myemail*',myemail) # Email
contents=contents.replace('*project*',project_name) # Project name
open('eigenvalue.sbatch','w').write(contents)
