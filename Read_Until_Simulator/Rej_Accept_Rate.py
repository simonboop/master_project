# a script to look through the log files of Read_Until_Simulator.py
# and make logfile with the current rejection rata at the specific time stamp
# of each line
#need to find regex patterns
import re
#used to convert dictionary strings to dictionaries
import ast
import sys
#used to do some maths to find out rejection rate
import numpy as np


#names prefix of names of agloithms file names
algor = ["1b","NR","T","P","std_res"]

#loop through algorithms and find the file in the folder sys.argv specified
for index, value in enumerate(algor):
#doesthis test use std _restrict
	if value == "std_res" and sys.argv[2] != "yes":
		pass
	else:
#there we're 3 repeats for each algor find each one
		for i in range(1,4):
			#make a file of that specific test with sepcific algortihm and repeat
			with open(value+str(i)+"_rates_"+str(sys.argv[1])+".txt", 'w') as line:
				line.write("Time,RejRate,AcceptRate")
#open and loop through log file, of that specific test and algorithm
			with open("./"+str(sys.argv[1])+"/"+str(value)+str(i)+"counts.txt") as f:
				last = None
#loop through file and look at each line store time, seq, obs and rej counts
				for line in f:
					time_pattern = re.findall('Time:(\d*\.\d*)', line)
					time = time_pattern
					Obs_pattern = re.findall('Obs:({.*?})', line)
					obsdict = ast.literal_eval(Obs_pattern[0])
					Seq_pattern = re.findall('Seq:({.*?})', line)
					seqdict = ast.literal_eval(Seq_pattern[0])
					Rej_pattern = re.findall('Rej:({.*?})', line)
					rejdict = ast.literal_eval(Rej_pattern[0])
#work out this lines rates for each bin add it to new log file
					seq = sum(seqdict.values())
					obs = sum(obsdict.values())
					rej = sum(rejdict.values())
					with open(str(value)+str(i)+"_rates_"+str(sys.argv[1])+".txt", "a") as line:
						line.write("\n"+str(time[0])+","+str(float(rej)/float(obs))+","+str(float(seq)/float(obs)))
