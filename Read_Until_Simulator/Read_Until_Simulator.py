# A port of the sum_readutil.R argument addition
# A simulation of read normalisation program to be used
# in read until sequencing
#imports to help display progress can be removed alongside print statements
import sys
import time
#used to make arrays and sampling, used for random genration
import numpy as np
import random
#makeing graphs
import matplotlib.pyplot as plt
#makeing gif
from wand.image import Image as Image_wand
from wand.api import library
#argument parser
import argparse
#make directories / files
import os
#regex for makeing graphs
import re
#convert str to dict
import ast


#makes the array which represents the the sequencing events as reads pass randomly through the pore
class amplicon_events:
    def __init__(self, num_amp=3, amp_prob=[0.01,0.2,0.1], event_num=100000, rand_len=True):
        #potential bases used
        self.bases = ["A", "C", "G", "T"]
        self.lengths = rand_len
        self.num_amp = num_amp
        #make length list if none given
        if type(rand_len) is bool:
            self.lengths = [0] * self.num_amp
            for x in range(len(self.lengths)):
                self.lengths[x] = random.randint(500, 2000)
        #list each type amplicon with its randomly assigned sequence
        self.types = [0]*self.num_amp
        for x in range(len(self.types)):
            self.types[x] = np.random.choice(self.bases, self.lengths[x], p=[0.25, 0.25, 0.25, 0.25])
        self.event_num = event_num
        #make a numpy array of prob list so probs sum to 1 so it can be used later on to generate array
        self.nor_prob = np.array(amp_prob)
        self.nor_prob /= self.nor_prob.sum()
        print("Normalised probs "+str(self.nor_prob))
        #make event array with data given
        self.reads = np.random.choice(self.num_amp, self.event_num, replace=True, p=self.nor_prob)

# loops through array events attempts to normalise amplicon levels as it simulates mapping to genome
class Read_Until:
    def __init__(self, event_array, num_amp,stop,amplicon_lengths, algorithm, channels, time_para):
        #time paras
        self.time_para = time_para
        print(self.time_para.reject_t)
        #channels use only really to make text file count.txt look more like the real thing
        self.channels = channels
        #amplicon lengths
        self.amplicon_lengths = amplicon_lengths
        #coverage stop point
        if type(stop) is int:
            self.stop = [stop]
        else:
            self.stop = stop
        #for Min_proportion to deal with custom coverages
        self.factor = [float(cov)/float(min(self.stop)) for cov in self.stop]
        #bin to store seen events
        self.amp_bin = [1]*num_amp
        self.prob = [0]*num_amp
        self.accept_bin = [0]*num_amp
        self.reject_bin = [0]*num_amp
        #the value of the max freq of accept reads
        self.acpt_max = 0
        self.regc_max = 0
        #list to record amplicon freq for both accept and all for every 1000
        self.accept_record = []
        self.all_record = []
        self.rejected_record = []
        #path of text file to record data overwrite old data
        self.path = "./output/counts.txt"
        #record of time taken
        self.time_taken = 0
        with open(self.path, "w") as file:
            pass
        #choice of algorithm
        self.algorithm = algorithm
        counter = 0
        coverage_mean = self.stop[counter]
        customfactor = [1]
        for i in range(0, len(event_array)):
            event = event_array[i]
            #deal with custom depth
            if len(self.stop) > 1:
                counter = event
                coverage_mean = sum(self.stop)/len(self.stop)
                customfactor = []
                for cov in self.stop:
                    customfactor.append(float(cov)/float(coverage_mean))
            #time to reload
            self.time_taken += self.time_para.reload_t
            #time taken to intial sequence and match read
            self.time_taken += self.time_para.match_t + (250/self.time_para.seq_t)
            # add too bin
            self.amp_bin[event] = self.amp_bin[event]+1
            #update probabilities and inverse as a new read is added
            self.prob = [float(sum(self.amp_bin)) / float(x) for x in self.amp_bin]
            #in event read is too small to reject in time
            if (((self.amplicon_lengths[event] * 2) / self.time_para.seq_t)-(self.time_para.match_t+(250/self.time_para.seq_t))) < 0:
                self.accept_bin[event] = self.accept_bin[event] + 1
                #take away for time saved with small reads
                self.time_taken += (((float(self.amplicon_lengths[event]) * float(2)) / float(self.time_para.seq_t))-float((self.time_para.match_t+(250/self.time_para.seq_t))))
            #else algorithm
            #no read until
            elif self.algorithm == "no_readuntil":
                    self.accept_bin[event] = self.accept_bin[event] + 1
                    #time taken to sequence the rest of the read
                    self.time_taken += (((float(self.amplicon_lengths[event]) * float(2)) / float(self.time_para.seq_t))-float((self.time_para.match_t+(250/self.time_para.seq_t))))
            #when over threshold
            elif self.accept_bin[event] >= self.stop[event]:
                self.reject_bin[event] = self.reject_bin[event] + 1
                #time taken to reject a read
                self.time_taken += self.time_para.reject_t
            #prob algorithm
            elif self.algorithm == "prob":
                #work out how likely we should add the amplicon to the tally
                self.accept_prob = float(self.prob[event]) / float(max(self.prob))
                #roll a fake dice if the accept_prob is greater add to the accepted bin
                #random beat
                print(self.accept_prob)
                if ((self.accept_prob * customfactor[counter])> random.uniform(0, 1)):
                #hardset beat doesnt work
                #if (self.accept_prob > 0.5):
                    self.accept_bin[event] = self.accept_bin[event] + 1
                    #time taken to sequence the rest of the read
                    self.time_taken += (((float(self.amplicon_lengths[event]) * float(2)) / float(self.time_para.seq_t))-float((self.time_para.match_t+(250/self.time_para.seq_t))))
                else:
                    self.reject_bin[event] = self.reject_bin[event] + 1
                #time taken to reject a read
                    self.time_taken += self.time_para.reject_t
            # default threshold algorithm
            elif self.algorithm == "threshold":
                self.accept_bin[event] = self.accept_bin[event] + 1
                #time taken to sequence the rest of the read
                self.time_taken += (((float(self.amplicon_lengths[event]) * float(2)) / float(self.time_para.seq_t))-float((self.time_para.match_t+(250/self.time_para.seq_t))))
            # Min_Proportion algorithm coverage is not aloud to be higher 1 > than
            # the lowest coverage
            elif self.algorithm == "Min_Proportion":
                #deal with custom depth
                if min(self.stop)+1 != min(self.accept_bin)+1:
                    #for x in range(0,len(self.accept_bin)):
                        #if self.accept_bin[event] > (self.accept_bin[x])+1 and self.stop[event] == min(self.stop):
                            #temp_decision= "rej"
                    #if self.accept_bin[self.stop.index(min(self.stop))] == 0:
                        #for x in range(0,len(self.accept_bin)):
                            #if self.accept_bin[event] > (self.accept_bin[x])+1 and self.stop[event] == min(self.stop):
                                #temp_decision= "rej"
                    #else:
                    temp_decision = "acc"
                    #temp coverage limits based on the current minimum amplcion count
                    temp_coverage = [self.factor[index] * self.accept_bin[self.factor.index(min(self.factor))] for index, value in enumerate(self.accept_bin)]
                    #if choosing for minmum coverage amplicon
                    if event == self.factor.index(min(self.factor)):
                        #check whether amplicon accepted bin of other amplicons(not min) are reached
                        for index, value in enumerate(temp_coverage):
                            if index == event:
                                pass
                            #if not reject
                            elif value > self.accept_bin[index]:
                                temp_decision= "rej"
                    #if choosing for other coverage amplicons
                    else:
                        #have the other amplicons reached temp coverage
                        for index, value in enumerate(temp_coverage):
                            if value < self.accept_bin[index]:
                                temp_decision= "rej"
                        #is this bin in question lower than the desire amplicon
                        if self.accept_bin[event] <= temp_coverage[event]:
                            temp_decision = "acc"
                    if temp_decision == "acc":
                        self.accept_bin[event] += 1
                        #time taken to sequence the rest of the read
                        self.time_taken += (((float(self.amplicon_lengths[event]) * float(2)) / float(self.time_para.seq_t))-float((self.time_para.match_t+(250/self.time_para.seq_t))))
                    else:
                        self.reject_bin[event] = self.reject_bin[event] + 1
                        #time taken to reject a read
                        self.time_taken += self.time_para.reject_t
                else:
                    self.accept_bin[event] += 1
                    #time taken to sequence the rest of the read
                    self.time_taken += (((float(self.amplicon_lengths[event]) * float(2)) / float(self.time_para.seq_t))-float((self.time_para.match_t+(250/self.time_para.seq_t))))
            # reads only added if it does not increase standard deviation above a certain level
            # standard deviation restrict
            elif self.algorithm == "std_restrict":
                #deal with custom depth
                if min(self.stop)+1 != min(self.accept_bin)+1:
                    temp_decision = "acc"
                    accpet_copy = list(self.accept_bin)
                    accpet_copy[event] += 1
                    accpet_copy = np.asarray(accpet_copy)
                    if np.std(accpet_copy) > self.stop[0]/20:
                        temp_decision = "rej"
                    if temp_decision == "acc":
                        self.accept_bin[event] += 1
                        #time taken to sequence the rest of the read
                        self.time_taken += (((float(self.amplicon_lengths[event]) * float(2)) / float(self.time_para.seq_t))-float((self.time_para.match_t+(250/self.time_para.seq_t))))
                    else:
                        self.reject_bin[event] = self.reject_bin[event] + 1
                        #time taken to reject a read
                        self.time_taken += self.time_para.reject_t
                else:
                    self.accept_bin[event] += 1
                    #time taken to sequence the rest of the read
                    self.time_taken += (((float(self.amplicon_lengths[event]) * float(2)) / float(self.time_para.seq_t))-float((self.time_para.match_t+(250/self.time_para.seq_t))))

            #appending text file
            with open(self.path, "a") as text_file:
                temp_obs = {index+1:value-1 for index, value in enumerate(self.amp_bin)}
                temp_acc = {index+1:value for index, value in enumerate(self.accept_bin)}
                temp_rej = {index+1:value for index, value in enumerate(self.reject_bin)}
                text_file.write("INFO:Channel:"+str(random.randint(0+1,self.channels))+",Read:"+str(i)+",Time:"+str(self.time_taken)+".0,Obs:"+str(temp_obs)+",Rej:"+str(temp_rej)+",Seq:"+str(temp_acc)+"\n")
            #break loop when coverage reached in each amplicon
            if self.stop is not None:
                if len(self.stop) == 1:
                    if all(amplicon >= self.stop for amplicon in self.accept_bin):
                        #readjust for the 1 given to amp bin for free at the begining
                        self.amp_bin = [x-1 for x in self.amp_bin]
                        break
                elif len(self.stop) > 1:
                    if all(amplicon >= self.stop[index] for index,amplicon in enumerate(self.accept_bin)):
                        #readjust for the 1 given to amp bin for free at the begining
                        self.amp_bin = [x-1 for x in self.amp_bin]
                        break

# Graph makeing class (png+gif) of Read_Until Data
class Graphy:
    #option is type a graph wanted
    def __init__(self,filename,time_interval):
        self.filename = filename
        self.time_interval = time_interval
        self.eventtimes = []
        self.obs_dict_history = []
        self.seq_dict_history = []
        self.rej_dict_history = []
        #regex to get data needed to make graphs taken from text file so this class can be used in real data
        with open(self.filename, 'r') as logline:
            lines = logline.readlines()
            for line in lines:
                time_pattern = re.findall('Time:(\d*\.\d*)', line)
                self.eventtimes.append(time_pattern[0])
                Obs_pattern = re.findall('Obs:({.*?})', line)
                self.obs_dict_history.append(ast.literal_eval(Obs_pattern[0]))
                Seq_pattern = re.findall('Seq:({.*?})', line)
                self.seq_dict_history.append(ast.literal_eval(Seq_pattern[0]))
                Rej_pattern = re.findall('Rej:({.*?})', line)
                self.rej_dict_history.append(ast.literal_eval(Rej_pattern[0]))
        #take records nearest to time divisions
        time_counter = 0
        self.obs_dict_needed =[]
        self.seq_dict_needed =[]
        self.rej_dict_needed =[]
        for index,value in enumerate(map(float, self.eventtimes)):
            if value > time_counter:
                time_counter += self.time_interval
                #get obs rej seq event that corrospond to index
                self.obs_dict_needed.append(self.obs_dict_history[index])
                self.rej_dict_needed.append(self.rej_dict_history[index])
                self.seq_dict_needed.append(self.seq_dict_history[index])
        #stats needed to make graph
        self.max_seq = max([i for i in self.seq_dict_needed[-1].values()])
        self.max_rej = max([i for i in self.rej_dict_needed[-1].values()])
    #makeing graphs
    def png(self, option,acc_rej):
        #barcharts
        if option == "Bar":
            if acc_rej == "acc":
                pool = self.seq_dict_needed
                max = self.max_seq
                name = "Acceptance"
            else:
                pool = self.rej_dict_needed
                max = self.max_rej
                name = "Rejection"
            for index, event in enumerate(pool):
                plt.xlabel('Amplicon')
                plt.ylabel('Count')
                plt.axis([1, len(pool[0])+1, 0, max])
                plt.title(name + " Amplicon Freq After Read Until")
                temp = []
                for key, value in event.iteritems():
                    temp.append(value)
                plt.bar(range(1,len(pool[0])+1), temp,label='Bars', width=0.9)
                plt.savefig("./output/"+name + "_Bar" + str(index) + ".png")
                plt.close()
            self.gif(name + "_Bar", len(pool))
        #dot plot
        if option == "Dot":
            if acc_rej == "acc":
                pool = self.seq_dict_needed
                name = "Acceptance"
            else:
                pool = self.rej_dict_needed
                name = "Rejection"
            rate_index = []
            count_index =[]
            for index, event in enumerate(self.obs_dict_needed):
                rate = float(sum(pool[index].values()))/float(sum(self.obs_dict_needed[index].values()))
                rate_index.append(rate)
                count= sum(self.obs_dict_needed[index].values())
                count_index.append(count)
            rate_temp = []
            count_temp = []
            count_max = count_index[-1]
            for index, event in enumerate(self.obs_dict_needed):
                for i in range(0,index+1):
                    rate_temp.append(rate_index[i])
                    count_temp.append(count_index[i])
                plt.plot(count_temp, rate_temp, 'ro')
                plt.xlabel('nReads')
                plt.ylabel(name+"Rate")
                plt.axis([1, (float(count_max)), 0, 1])
                plt.savefig("./output/"+ name + "_Dot" + str(index) + ".png")
                plt.close()
    #make gifs
    def gif(self,target,lastimgnum):
        with Image_wand() as new_gif:
            file_names = []
            for x in range(0, lastimgnum):
                file_names.append("./output/"+str(target) + str(x) + ".png")
            library.MagickSetOption(new_gif.wand, 'dispose', 'background')
            for x in file_names:
                library.MagickReadImage(new_gif.wand, x)
            new_gif.save(filename=("./output/"+str(target)+".gif"))

#class to holder the time parameters for this run to be passed to other objects
class timeholder():
    def __init__(self, reload_t, reject_t, match_t, seq_t):
        self.reload_t = reload_t
        self.reject_t = reject_t
        self.match_t = match_t
        self.seq_t = seq_t




#intialise program
if __name__ == "__main__":
    #make directory if does not exist to store images
    if not os.path.exists("./output"):
        os.makedirs("./output")
    #parse arguments
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--number", help ="(int) of the types of amplicons", type=int, required=True)
    parser.add_argument("-p", "--probs", help ="(list(float/int)) of the probs/ratios of each amplicon, e.g 1 2 3...", type=float, required=True, nargs='+')
    parser.add_argument("-rp", "--read_pool", help ="(int) size of the read_pool simulation will be sampling", type=int, required=True)
    parser.add_argument("-al", "--amplicon_lengths", help="(list(int)) the lengths of amplicons, e.g 1 2 3...", type=int, required=True, nargs='+')
    parser.add_argument("-c", "--coverage", help="(list) coverage desired across amplicons, coverage can be applied e,g 200 100 200.Custom depth will not work with std_restrict", required=True, nargs='+')
    parser.add_argument("-ac", "--algorithm_choice", help="choose an algorithm e.g threshold, prob, 1_by_1, no_readuntil, std_restrict", type=str, required=True)
    parser.add_argument("-nc", "--number_channels", help="(int/float) the number of channels in simulation", type=int, required=True)
    parser.add_argument("-mt", "--match_time", help="time taken to match 250bp squiggle to ref, default 1s", type=float)
    parser.add_argument("-st", "--seq_time", help="(float)seq speed default 250bp per second", type=float)
    parser.add_argument("-rt", "--reload_time", help="time taken to reload new read, deafult = 1s", type=float)
    parser.add_argument("-rjt", "--reject_time", help="time taken to reject a read, default = 1s", type=float)
    args = parser.parse_args()
    #checking user input
    pass_ticket = "pass"
    #optional arguments
    if args.reload_time is not None:
        reload_time = args.reload_time
    else:
        reload_time = 1
    if args.reject_time is not None:
        reject_time = args.reject_time
    else:
        reject_time = 1
    if args.match_time is not None:
        match_time = args.match_time
    else:
        match_time = 1
    if args.seq_time is not None:
        seq_time = args.seq_time
    else:
        seq_time = 250
    time_para = timeholder(reload_time,reject_time,match_time,seq_time)
    algorithms = ["threshold", "prob", "Min_Proportion", "no_readuntil", "std_restrict"]
    if type(args.coverage) is int:
        pass
    if type(args.coverage) is list:
        if args.number != len(args.coverage):
            pass_ticket="fail"
            print("number of amplicons must == len(coverage) when custom depth applied")
        try:
            args.coverage = map(int, args.coverage)
        except:
            print("unable to convert custom depth to int a item/s in custom depth list are not numbers")
            pass_ticket = "fail"
    else:
        print("invalid coverage type must be list or int")
        pass_ticket = "fail"
    if args.algorithm_choice not in algorithms:
        pass_ticket = "fail"
        print("algorithm not found")
    if args.number != len(args.probs) or args.number != len(args.amplicon_lengths):
        pass_ticket = "fail"
        print("number of amplicons must == len(probs) and len(amplicon_lengths)")
    #warnings
    for amp in args.amplicon_lengths:
        if (((amp * 2) / 250)-2) < 0:
            print("Warning! read length of "+ str(amp)+ " is to small to be rejected by read until")
    if pass_ticket == "pass":
        print("Correct Formatting of Input... \nInitializing...\nMakeing event array(readpool) for read until to sample from...")
        events = amplicon_events(args.number,args.probs, args.read_pool, args.amplicon_lengths)
        print("Complete Array Made...")
        print("Simulating Read Until...")
        run = Read_Until(events.reads, args.number, args.coverage, args.amplicon_lengths,args.algorithm_choice,args.number_channels, time_para)
        print("Simulation Complete...")
        print("Results:")
        print("Obs Reads bin "+ str(run.amp_bin))
        print("Accepted Reads bin"+ str(run.accept_bin))
        print("Rejected Reads bin"+ str(run.reject_bin))
        print("Time Taken in seconds == "+str(run.time_taken)+"s")
        #ask user whether to make graphs
        proceed = "stop"
        while True:
            user_input = str(raw_input("Make Graphs (y/n)?:"))
            if user_input == "y":
                proceed = "pass"
                break
            elif user_input == "n":
                proceed = "stop"
                print("bye")
                break
        #makeing graphs
        if proceed == "pass":
            proceed = "stop"
            while True:
                user_input = str(raw_input("Make a graph every x seconds of the simulation...\ninput x: "))
                try:
                    user_input = int(user_input)
                    if type(user_input) is int:
                        if user_input < run.time_taken:
                            time4graph = user_input
                            proceed = "pass"
                            break
                        else:
                            print("x must be < time_taken which is "+str(run.time_taken))
                    else:
                        print("x must be an int")
                except:
                    print("x must be an int")

            if proceed == "pass":
                print("Makeing Graphs...")
                graph = Graphy("./output/counts.txt",time4graph)
                graph.png("Bar","acc")
                print("Bar Acc made")
                graph.png("Dot","acc")
                print("Dot Acc made")
                graph.png("Bar","rej")
                print("Bar Rej made")
                graph.png("Dot","rej")
                print("Dot Rej made")
                print("Makeing gifs")
                graph.gif("Acceptance_Bar",len(graph.obs_dict_needed))
                print("1/4 done")
                graph.gif("Rejection_Bar",len(graph.obs_dict_needed))
                print("2/4 done")
                graph.gif("Acceptance_Dot",len(graph.obs_dict_needed))
                print("3/4 done")
                graph.gif("Rejection_Dot",len(graph.obs_dict_needed))
                print("4/4 done")
                print("All done bye!")












