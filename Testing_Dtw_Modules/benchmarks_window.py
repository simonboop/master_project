#test bench to test speed and accuracy of different modules of dtw with a window search
import numpy as np
#time
from datetime import datetime

#modules
#mlpy
import mlpy
#ftw
from fastdtw import fastdtw
from scipy.spatial.distance import cityblock
#ucrdtw
import _ucrdtw
#pydtw
from pydtw import dtw1d
#cdtw
from cdtw import pydtw
#dtw
from dtw import dtw

#log file to store timings and distances
with open("bench_log.txt", "w") as text_file:
    text_file.write("Amplicon,Window,Module,Dist,Time,Start,End")
    pass

#importing text files and turning them back into numpy arrays
Amp1F,Amp2F,Amp3F,Amp4F,Amp5F,Amp6F,Amp7F,Amp8F,Amp9F,Amp10F,Amp11F,RefF = 0,0,0,0,0,0,0,0,0,0,0,0
queries = ["Amp1F.txt","Amp2F.txt","Amp3F.txt","Amp4F.txt","Amp5F.txt","Amp6F.txt","Amp7F.txt","Amp8F.txt","Amp9F.txt","Amp10F.txt","Amp11F.txt","RefF.txt"]
querynames = [Amp1F,Amp2F,Amp3F,Amp4F,Amp5F,Amp6F,Amp7F,Amp8F,Amp9F,Amp10F,Amp11F,RefF]

#amplicons
for index,value in enumerate(queries):
    squigglefile = open(value, "r")
    querynames[index] = squigglefile.read().split('\n')
    del(querynames[index][-1])
    querynames[index] = map(float, querynames[index])
    querynames[index] = np.asarray(querynames[index])

#window
count = 0
for window in range(0,(len(querynames[-1])/250)):
    y = querynames[-1][count:count+500]
    print(querynames[-1][count:count+500])
#test each amplicon against F reference
    for amp in range(0,11):
        #amplicons 1F
        x = querynames[amp]
        #mlpy
        timeb = datetime.now()
        mlpystddist, mlpystdcost, mlpystdpath = mlpy.dtw_std(x, y, dist_only=False)
        timet = datetime.now() - timeb
        print("mlpy complete on amp "+str(amp+1))
        with open("bench_log.txt", "a") as text_file:
            text_file.write("\n"+str(amp+1)+","+str(window+1)+",mlpy,"+str(mlpystddist)+","+str(timet.microseconds)+','+str(mlpystdpath[1][0]+count)+','+str(mlpystdpath[1][-1]+count))
        path1 = np.savetxt('paths/'+"amp_"+str(amp+1)+"_window_"+str(window+1)+'_query_mlpy.txt',mlpystdpath[0],delimiter=',')
        path2 = np.savetxt('paths/'+"amp_"+str(amp+1)+"_window_"+str(window+1)+'_ref_mlpy.txt',mlpystdpath[1],delimiter=',')
        timeb = datetime.now()
        mlpysubdist, mlpysubcost, mlpysubpath = mlpy.dtw_subsequence(x, y)
        timet = datetime.now() - timeb
        print("mlpy sub complete on amp "+str(amp+1))
        with open("bench_log.txt", "a") as text_file:
            text_file.write("\n"+str(amp+1)+","+str(window+1)+",mlpy_sub,"+str(mlpysubdist)+","+str(timet.microseconds)+','+str(mlpysubpath[1][0]+count)+','+str(mlpysubpath[1][-1]+count))
        path1 = np.savetxt('paths/'+"amp_"+str(amp+1)+"_window_"+str(window+1)+'_query_mlpysub.txt',mlpysubpath[0],delimiter=',')
        path2 = np.savetxt('paths/'+"amp_"+str(amp+1)+"_window_"+str(window+1)+'_ref_mlpysub.txt',mlpysubpath[1],delimiter=',')
        #ftw
        timeb = datetime.now()
        ftwdistance, ftwpath = fastdtw(x, y, dist=cityblock)
        timet = datetime.now() - timeb
        print("ftw complete on amp "+str(amp+1))
        ftwpath = zip(*ftwpath)
        ftwpath = tuple(map(tuple, ftwpath))
        ftwpath = np.asarray(ftwpath)
        with open("bench_log.txt", "a") as text_file:
            text_file.write("\n"+str(amp+1)+","+str(window+1)+",ftw,"+str(ftwdistance)+","+str(timet.microseconds)+','+str(ftwpath[1][0]+count)+','+str(ftwpath[1][-1]+count))
        path1 = np.savetxt('paths/'+"amp_"+str(amp+1)+"_window_"+str(window+1)+'_query_ftw.txt',ftwpath[0],delimiter=',')
        path2 = np.savetxt('paths/'+"amp_"+str(amp+1)+"_window_"+str(window+1)+'_ref_ftw.txt',ftwpath[1],delimiter=',')
        #ucrdtw
        timeb = datetime.now()
        ucrloc, ucrdist = _ucrdtw.ucrdtw(y, x, 0.05)
        timet = datetime.now() - timeb
        print("ucr complete on amp "+str(amp+1))
        with open("bench_log.txt", "a") as text_file:
            text_file.write("\n"+str(amp+1)+","+str(window+1)+",ucr,"+str(ucrdist)+","+str(timet.microseconds)+','+str(ucrloc+count)+','+str(ucrloc+count))
        with open('paths/'+str(amp+1)+","+str(window+1)+'_loc_ucr.txt', "w") as text_file:
            text_file.write(str(ucrloc))
        #cydtw
        timeb = datetime.now()
        cdtw_master = pydtw.dtw(x,y,pydtw.Settings(step = 'p0sym',           #Sakoe-Chiba symmetric step with slope constraint p = 0
                                                   window = 'palival',       #type of the window
                                                   param = 2.0,              #window parameter
                                                   norm = False,             #normalization
                                                   compute_path = True))
        timet = datetime.now() - timeb
        print("cdtw complete on amp "+str(amp+1))
        cdtwpath1, cdtwpath2 = zip(*cdtw_master.get_path())
        cdtwpath1 = np.asarray(cdtwpath1)
        cdtwpath2 = np.asarray(cdtwpath2)
        with open("bench_log.txt", "a") as text_file:
            text_file.write("\n"+str(amp+1)+","+str(window+1)+",cdtw,"+str(cdtw_master.get_dist())+","+str(timet.microseconds)+','+str(cdtwpath2[0]+count)+','+str(cdtwpath2[-1]+count))
        path1 = np.savetxt('paths/'+"amp_"+str(amp+1)+"_window_"+str(window+1)+'_query_cdtw.txt',cdtwpath1,delimiter=',')
        path2 = np.savetxt('paths/'+"amp_"+str(amp+1)+"_window_"+str(window+1)+'_ref_cdtw.txt',cdtwpath2,delimiter=',')
        #dtw
        timeb = datetime.now()
        dtwdist, dtwcost, dtwacc, dtwpath = dtw(x,y, dist=cityblock)
        timet = datetime.now() - timeb
        print("dtw complete on amp "+str(amp+1))
        dtwpath1 = np.asarray(dtwpath[0])
        dtwpath2 = np.asarray(dtwpath[1])
        with open("bench_log.txt", "a") as text_file:
            text_file.write("\n"+str(amp+1)+","+str(window+1)+",dtw,"+str(dtwdist)+","+str(timet.microseconds)+','+str(dtwpath2[0]+count)+','+str(dtwpath2[-1]+count))
        path1 = np.savetxt('paths/'+"amp_"+str(amp+1)+"_window_"+str(window+1)+'_query_dtw.txt',dtwpath1,delimiter=',')
        path2 = np.savetxt('paths/'+"amp_"+str(amp+1)+"_window_"+str(window+1)+'_ref_dtw.txt',dtwpath2,delimiter=',')
    count += 250






