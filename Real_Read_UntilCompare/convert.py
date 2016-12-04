# a script log files from aReaduntil.py with a time format
#that starts from zero
import re
import ast
import sys

#new file to hold conversion
with open("./"+str(sys.argv[1])+"_convert", 'w') as line:
    pass

#loop through log file line by line change the time stamp print
with open(sys.argv[1]) as f:
    count = 0
    smallest_time = 0
    #find the smallest time stamp
    for line in f:
        time_pattern = re.findall('Time:(\d*\.\d*)', line)
        if count == 0:
            smallest_time = float(time_pattern[0])
        elif smallest_time > float(time_pattern[0]):
            smallest_time = float(time_pattern[0])
        count += 1
    print(smallest_time)
#take away from each time stamp
with open(sys.argv[1]) as f:
    #minus it from each stamp
    for line in f:
        time_pattern = re.findall('Time:(\d*\.\d*)', line)
        #minus smalles time convert to seconds
        newtime = (float(time_pattern[0]) - float(smallest_time))/float(1000)
        #find the rest of the string
        Obs_pattern = re.findall('(Obs:{.*?})', line)
        Seq_pattern = re.findall('(Seq:{.*?})', line)
        Rej_pattern = re.findall('(Rej:{.*?})', line)
        Channel_pattern = re.findall('Channel:\d*',line)
        Read_pattern = re.findall('Read:\d*',line)
        #append change into new file
        with open("./"+str(sys.argv[1])+"_convert", 'a') as line:
            line.write("INFO:"+Channel_pattern[0]+","+Read_pattern[0]+",Time:"+str(newtime)+","+Obs_pattern[0]+","+str(Rej_pattern[0])+","+str(Seq_pattern[0])+"\n")



