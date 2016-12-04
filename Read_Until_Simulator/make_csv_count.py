import re
import ast
import sys

NR1 = ("./"+str(sys.argv[1])+"/NR1counts.txt")
NR2 = ("./"+str(sys.argv[1])+"/NR2counts.txt")
NR3 = ("./"+str(sys.argv[1])+"/NR3counts.txt")
P1 = ("./"+str(sys.argv[1])+"/P1counts.txt")
P2 = ("./"+str(sys.argv[1])+"/P2counts.txt")
P3 = ("./"+str(sys.argv[1])+"/P3counts.txt")
T1 = ("./"+str(sys.argv[1])+"/T1counts.txt")
T2 = ("./"+str(sys.argv[1])+"/T2counts.txt")
T3 = ("./"+str(sys.argv[1])+"/T3counts.txt")
b1 = ("./"+str(sys.argv[1])+"/1b1counts.txt")
b2 = ("./"+str(sys.argv[1])+"/1b2counts.txt")
b3 = ("./"+str(sys.argv[1])+"/1b3counts.txt")
std1 = ("./"+str(sys.argv[1])+"/std_res1counts.txt")
std2 = ("./"+str(sys.argv[1])+"/std_res2counts.txt")
std3 = ("./"+str(sys.argv[1])+"/std_res3counts.txt")

with open("./counts/"+str(sys.argv[1])+"_count", 'w') as line:
    line.write("Repeat,Algorithm,Time,Obs,Seq,Rej")

if sys.argv[2] == "yes":
    files = [NR1,NR2,NR3,P1,P2,P3,T1,T2,T3,b1,b2,b3,std1,std2,std3]
    files_names =["NR","NR","NR","P","P","P","T","T","T","b","b","b","std","std","std"]
else:
    files = [NR1,NR2,NR3,P1,P2,P3,T1,T2,T3,b1,b2,b3,std1,std2,std3]
    files_names =["NR","NR","NR","P","P","P","T","T","T","b","b","b"]
    
    
temp = 0
for index, count in enumerate(files):
    temp += 1
    time = 0
    obs = {}
    seq = {}
    rej = {}
    with open(count) as f:
        last = None
        for last in(line for line in f if line.rstrip('\n')):
                time_pattern = re.findall('Time:(\d*\.\d*)', last)
                time = time_pattern
                Obs_pattern = re.findall('Obs:({.*?})', last)
                obs = ast.literal_eval(Obs_pattern[0])
                Seq_pattern = re.findall('Seq:({.*?})', last)
                seq = ast.literal_eval(Seq_pattern[0])
                Rej_pattern = re.findall('Rej:({.*?})', last)
                rej = ast.literal_eval(Rej_pattern[0])
        with open("./counts/"+str(sys.argv[1])+"_count", "a") as line:
            line.write("\n"+str(temp)+","+files_names[index]+","+str(time[0])+","+str(sum(obs.values()))+","+str(sum(seq.values()))+","+str(sum(rej.values())))

