#!/home/migrau/env/.uma2/bin/python
##############################################################
#  Script     : createGraphfiles.py
#  Author     : uma sangumathi
#  Date       : 28/04/2015
#  Last Edited: 24/05/2015, uma
#  Description: Reformat the files to generate graphs
##############################################################
# Purpose:
#  Specific for Influenza virus with 8 segments...
# Requirements:
#  1. Biopython module
#  2. Samtools

#############################################################

#HOW TO CALL
#python ~/pipelines/miguel/paired/tools/createGraphfiles_Full_17Nov.py /home/migrau/pipelines/miguel/paired/out_cutadapt_iseq/run3_2Dec19_updatedDB/S23/mig/N10031924_S23.fa N10031924_S23.R1.paired.fastq N10031924_S23.R2.paired.fastq

import os
import sys
import glob
from Bio import SeqIO

def findGene(s):
  if "(HA)" in s:
    return "HA"
  elif "(MP)" in s:
    return "MP"
  elif "(NA)" in s:
    return "NA"
  elif "(NS)" in s:
    return "NS"
  elif "(NP)" in s:
    return "NP"
  elif "(PA)" in s:
    return "PA"
  elif "(PB1)" in s:
    return "PB1"
  elif "(PB2)" in s:
    return "PB2"
  elif "(NEP)" in s:
    return "NS"
  elif "(M1)" in s:
    return "MP"

# To create coverage graph....
subgraph = ["HA", "NA", "MP", "PA", "PB1", "PB2", "NP", "NS"]
graph_order = [("PB2", "hs15"), ("PB1", "hs13"), ("PA", "hs9"), ("HA", "hs8"), ("NP", "hs7"), ("NA","hs17"), ("MP","hs21"), ("NS", "hs1")]
gene_size = {"PB2": 2300, "PB1" : 2300, "PA": 2230, "NP": 1570, "NA":1460, "HA":1760 , "MP": 1030 , "NS":870 }

# if len(sys.argv)==2:
#   xd = sys.argv[-1]
#   print xd
# else:
#   print "Incorrect arguments: enter outout directory"
#   sys.exit(0)
#for x in glob.glob(xd + "/*-sort.bam"):
ref_file = sys.argv[1] #FASTA REF RENAMED
R1=sys.argv[2] #R1
R2=sys.argv[3] #R2
ref_dict = {}
#ref_file = x.replace("-sort.bam",".fa")
print (ref_file)
dirref = os.path.dirname(ref_file)
nameref = os.path.basename(ref_file).split(".")[0]
print(nameref)
newfa=""
for r in SeqIO.parse(ref_file, "fasta"):
  print(r.id)
  #nid = r.description.split("(")[1]+"|"+findGene(r.description)+"|"
  nid = r.description.replace("_","|")
  #print(nid)
  #ref_dict[r.id.split("|")[1]] = [r.id, r.id.split("|")[1] , str(1) , str(len(str(r.seq)))]
  ref_dict[nid.split("|")[1]] = [nid, nid.split("|")[1] , str(1) , str(len(str(r.seq)))]
  newfa+=">"+nid+"\n"+str(r.seq)+"\n"
with open (dirref+"/"+nameref+"_reftemp.fa",'w') as fi:
  fi.write(newfa)
ref_genes = [k for k in ref_dict.keys() if k in subgraph]
#redo alignment with the corrected-names reference
os.system("bwa index "+dirref+"/"+nameref+"_reftemp.fa")
bwa_cmd = "bwa mem "+dirref+"/"+nameref+"_reftemp.fa -t 2 " + R1 +" "+ R2 + "> " + dirref +"/"+nameref+ "_depth.sam"
#bwa_cmd = "bwa aln "+dirref+"/"+nameref+"_reftemp.fa -t 2 -1 " + R1 +" -2 "+ R2 + "> " + dirref +"/"+nameref+ "_depth.sam"
print(bwa_cmd)
os.system(bwa_cmd)
sambam="samtools view -bT "+dirref+"/"+nameref+"_reftemp.fa " + dirref +"/"+nameref+ "_depth.sam | samtools sort > " + dirref +"/"+nameref+ "_depth.bam"
print(sambam)
os.system(sambam)
x=dirref+"/"+nameref+"_depth.bam"

# Write coverage per position
#sam_cmd = "samtools mpileup -f " + ref_file + " " + x + "  >  " + x.replace(".bam", ".pileup")
#miguel edit: count-orphans and no max-depth:
sam_cmd = "samtools mpileup -d 0 -A -f " + dirref+"/"+nameref+"_reftemp.fa" + " " + dirref+"/"+nameref+"_depth.bam  >  " + x.replace(".bam", ".pileup")
print (sam_cmd)
if not os.path.isfile(x.replace(".bam", ".pileup")):os.system(sam_cmd)
cov_cmd = " awk '{print $1, $2, $2, $4}'  " + x.replace(".bam", ".pileup")  + "  >  " + x.replace(".bam", ".coverage")
print (cov_cmd)
os.system(cov_cmd)
for ll in ref_dict:  # Assumin the gene name is 'B/Sydney/13/2014|NS|'
  lll = ref_dict[ll][0].split("|")
  cmdl = "sed -e 's/" + lll[0].replace("/", "\/") + "//g' " + x.replace(".bam", ".coverage")  + " | sed -e 's/|//g'  > " + x.replace(".bam", "-coverage.txt")
  print (cmdl)
  os.system(cmdl)

# find mean coverage
for i in ref_genes:
  cmd = "grep '" + i + "' " +  x.replace(".bam", ".coverage") + " | awk '{SUM += $4 }  END { print SUM/ NR}' > "+dirref+"/"+nameref+"_reftemp.txt"
  os.system("grep '" + i + "' " +  x.replace(".bam", ".coverage") + " | awk '{SUM += $4 }  END { print SUM, NR, SUM/ NR}'")
  print (i)
  os.system(cmd)
  if os.path.isfile(dirref+"/"+nameref+"_reftemp.txt"):
    for k in open(dirref+"/"+nameref+"_reftemp.txt", "r"):
      ref_dict[i].append(k.strip("\n") + "x")
    print ("Sum: " , k)
    os.system("rm "+dirref+"/"+nameref+"_reftemp.txt")

for i in gene_size:
  f = 0
  for j in ref_dict:
    if i == j: f = 1
  if f == 0:
    ref_dict[i] = [i, i, str(1), str(gene_size[i]), "0x"]
print (ref_dict)


# Write ideogram flu coords
with open( x.replace("_depth.bam","-flu_coords.txt"), "w") as f:
  for j in graph_order:
    for k in ref_dict:
      if j[0] == k:
        f.write("chr - " + ref_dict[k][1] + "\t" + j[0] + "\t" + ref_dict[k][2] + "\t" + ref_dict[k][3] + "\t" + j[1] + "\r\n")


# Write mean coverage file...
with open(x.replace("_depth.bam","-covX.txt"), "w") as f:
  for k in ref_dict:
    print (ref_dict[k])
    f.write( ref_dict[k][1] + "\t" +  ref_dict[k][2] + "\t" +  ref_dict[k][3] + "\t" +  ref_dict[k][4] + "\r\n")
