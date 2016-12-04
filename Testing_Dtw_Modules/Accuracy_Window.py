# a script to go through benchtests window log files and determine the window/position that a modules has predicted the squiggle to match
#boundries
amplicons = [[52,1737],[1738,3443],[3444,5114],[5115,6745],[6746,8408],[8409,10052],[10053,11723],[11724,13454],[13456,15091],[15092,16774],[16774,19150]]

#update composite time, lowest dist, loc of each window for each amplicon in the format [time,dist,start,end]
composite_results_mlpy = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
composite_results_mlpy_sub = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
composite_results_ftw = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
composite_results_ucr = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
composite_results_cdtw = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
composite_results_dtw = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
composite_results = {"mlpy":composite_results_mlpy,"mlpy_sub":composite_results_mlpy_sub,"ftw":composite_results_ftw,"ucr":composite_results_ucr,"cdtw":composite_results_cdtw,"dtw":composite_results_dtw}

#old log file with predicted locations
with open("bench_log.txt", 'r') as textfile:
    for index,line in enumerate(textfile):
        #miss headers
        if index != 0:
            #make list for each line strip newline character
            newline = line.split(",")
            if '\n' in newline[-1]:
                newline[-1] = newline[-1][:-1]
            #add time taken so far for amplicon
            #newline2 == algorithm name
            #newline0 == amplicon number
            #newline-3 == time
            #newline1 == window num
            #newline3 == dist
            #insertime
            composite_results[newline[2]][int(newline[0])-1][0] += float(newline[-3])
            #update with distance when a lower distance is found
            if float(composite_results[newline[2]][int(newline[0])-1][1]) == 0:
                #insert dist
                composite_results[newline[2]][int(newline[0])-1][1] = float(newline[3])
                #insert start pos
                composite_results[newline[2]][int(newline[0])-1][2] = float(newline[-2])
                #inser end pos
                composite_results[newline[2]][int(newline[0])-1][3] = float(newline[-1])
            #if the new distance is less than update dictionary
            elif float(composite_results[newline[2]][int(newline[0])-1][1]) > float(newline[3]):
                composite_results[newline[2]][int(newline[0])-1][1] = float(newline[3])
                composite_results[newline[2]][int(newline[0])-1][2] = float(newline[-2])
                composite_results[newline[2]][int(newline[0])-1][3] = float(newline[-1])

#new file to store Algorithm Amplicon Prediction Accuracy
with open("Algorithm_Accuracy.txt", 'w') as textfile:
    textfile.write("Algorithm,Amp1,Amp2,Amp3,Amp4,Amp5,Amp6,Amp7,Amp8,Amp9,Amp10,Amp11")

#newfile to store Algorithm times to predict location
with open("Algorithm_Time.txt", 'w') as textfile:
    textfile.write("Algorithm,Amp1,Amp2,Amp3,Amp4,Amp5,Amp6,Amp7,Amp8,Amp9,Amp10,Amp11")

#loop through location predicted for each amplicon by each algorithm
for algor in composite_results:
    for amp in range(0,11):
        #find mid point of location which amplicon is it in
        if algor == "ucr":
            if amplicons[amp][0] <= int(composite_results[algor][amp][-1]) <= amplicons[amp][1]:
                result = "Correct"
            else:
                result = "Incorrect"
        else:
            #find mid point by doing ((endpoint - startpoint)/2) + startpoint
            if amplicons[amp][0] <= (((composite_results[algor][amp][-1]-composite_results[algor][amp][-2])/2)+composite_results[algor][amp][-2]) <= amplicons[amp][1]:
                result = "Correct"
            else:
                result = "Incorrect"
        #fill accuracy log
        with open("Algorithm_Accuracy.txt", 'a') as textfile:
            if amp == 0:
                textfile.write("\n"+str(algor)+","+result)
            else:
                textfile.write(","+result)
        #fill time log
        with open("Algorithm_Time.txt", 'a') as textfile:
            if amp == 0:
                textfile.write("\n"+str(algor)+","+str(composite_results[algor][amp][0]))
            else:
                textfile.write(','+str(composite_results[algor][amp][0]))












