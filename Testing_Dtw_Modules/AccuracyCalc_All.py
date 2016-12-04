# a script to go through benchtests overall log files and determine if method correctly predicted amplicon
#boundries of amplicons
amplicons = [[52,1737],[1738,3443],[3444,5114],[5115,6745],[6746,8408],[8409,10052],[10053,11723],[11724,13454],[13456,15091],[15092,16774],[16774,19150]]
#log file to store new data
with open("bench_logallnew.txt", 'w') as textfile:
    textfile.write("Amplicon,Module,Dist,Time,Start,End,CorrectlyPredicted"+"\n")

#old log file with predicted locations
with open("bench_logall.txt", 'r') as textfile:
    for index,line in enumerate(textfile):
        #miss headers
        if index != 0:
            #make list for each line strip newlie character
            newline = line.split(",")
            if '\n' in newline[-1]:
                newline[-1] = newline[-1][:-1]
            #if its matched the query to ref globally or to big instant fail
            if int(newline[-1]) - int(newline[-2]) > 1000:
                result = "Incorrect"
            #otherwise see where mid point of that path is to determine which amplicon it has match it to
            else:
                correct = amplicons[int(newline[0])-1]
                #ucr handled different no path given
                if newline[1] == "ucr":
                    if correct[0] <= int(newline[-1]) <= correct[1]:
                        result = "Correct"
                    else:
                        result = "Incorrect"
                else:
                    mid = ((int(newline[-1]) - int(newline[-2]))/2)+int(newline[-2])
                    if correct[0] <= mid <= correct[1]:
                        result = "Correct"
                    else:
                        result = "Incorrect"
            newline.append(result)
            #write result in new log
            with open("bench_logallnew.txt", 'a') as textfile:
                textfile.write(','.join(newline)+"\n")







