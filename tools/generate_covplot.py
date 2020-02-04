#!/home/migrau/env/.uma2/bin/python

import os,sys
import glob
import csv
import fileinput

#PARAM 1. Coverage file
covfile = sys.argv[1] #N10031924_S23-covX.txt
nameref = os.path.basename(covfile).split("-covX")[0] #N10031924_S23
dirref = os.path.dirname(covfile)

with open(dirref+"/"+nameref+".conf", "w") as f:
  for line in fileinput.input("/home/yi-mo/pipelines/miguel/paired/tools/hist.conf", inplace=False):
    f.write(line.replace("00000", dirref+"/"+nameref).replace("XXXXXX", "10000")) #.replace("dirpath", plot_dir))

print ("circos -conf " + dirref+"/"+nameref+".conf")  
circos_cmd = "circos -conf " + dirref+"/"+nameref+".conf" 
os.system(circos_cmd)
