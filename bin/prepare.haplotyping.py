#!/usr/bin/env python

from __future__ import division
import sys

dpfile = sys.argv[1]
gtfile = sys.argv[2]
sampleInfoFile = sys.argv[3]

#1, parse sampleInfoFile
fsampleinfo = open(sampleInfoFile)
embryonames = []
for line in fsampleinfo:
	if line.startswith("Identity") or (not line):
		continue
	col = line.strip().split(",")
	if col[0].strip() == "Father":
		fathername = col[2].strip()
	elif col[0].strip() == "Mother":
		mothername = col[2].strip()
	elif col[0].strip() == "AffectedSon": 
		if col[1].strip() == "Yes":
			childname = col[2].strip()
			affected = "son"
		else:
			continue
	elif col[0].strip() == "AffectedDaughter":
		if col[1].strip() == "Yes":
			childname = col[2].strip()
			affected = "duaghter"
		else:
			continue
	elif col[0].strip() == "Embryo":
		embryonames.append(col[2].strip())
	else:
		print(line.strip()+"\n")
		print("Identity not matched")
		exit()
		
print("Father:", fathername)
print("Mother:", mothername)
print("Affected child:", childname, affected)
print("Embryos:", embryonames)

fproband = open("intermediate/proband", "w")
fproband.write(childname+"\n")
fproband.close()
#1, end


#2, save each sample's depth
dp = {}
fdp = open(dpfile)
for line in fdp:
	col = line.strip().split("\t")
	chrom = col[0]
	pos = col[1]
	rs = col[2]
	if rs == ".":
		continue
	ref = col[3]
	alt = col[4]
	dp[rs]={}
	for i in range(5,len(col)):
		s=col[i].split("=")[0] #sample name
		d=col[i].split("=")[1] # depth
		if d == ".":
			dp[rs][s] = 0
		else:
			dp[rs][s]= int(d)
		
fdp.close()
#2, end


#3, trio-phase
fhaploidPat = open("intermediate/father.snp.het.filter","w")
fhaploidMat = open("intermediate/mother.snp.het.filter","w")

flag_haplo_header = 1

gt = {}
fgt = open(gtfile)
for line in fgt:
	col = line.strip().split("\t")
	chrom = col[0]
	pos = int(col[1])
	rs = col[2]	
	if rs == ".":
		continue
	ref = col[3]
	alt = col[4]

	# filter father and mother sequencing depth
	if dp[rs][fathername]<3:
		continue
	if dp[rs][mothername]<3:
		continue
	gt[rs]={}
	
	for i in range(5,len(col)):
		s=col[i].split("=")[0] # sample name
		g=col[i].split("=")[1] # genotype
		if g == "./.":
			gt[rs][s] = "NA"
		else:
			gt[rs][s]= g
			
	if gt[rs][fathername] == gt[rs][mothername]:
		continue
	if gt[rs][fathername][0] == gt[rs][fathername][2] and gt[rs][mothername][0]==gt[rs][mothername][2]:
		continue

############## phasing ######################################
	husband = gt[rs][fathername]
	wife = gt[rs][mothername]
	husband_bases = husband.split("/")
	wife_bases = wife.split("/")
	
	if flag_haplo_header ==1:
		outline_haplo_headerPat = "#chr\tpos\thet1\tLH1\thet1Num\thet2\tLH2\thet2Num\tAll"
		outline_haplo_headerMat = "#chr\tpos\thet1\tLH1\thet1Num\thet2\tLH2\thet2Num\tAll"
		for i in range(5,len(col)):
			s=col[i].split("=")[0] #sample name
			if (s in embryonames) or (s==childname) :
				outline_haplo_headerPat = outline_haplo_headerPat +"\t"+ s +"\t"+ "LH"
				outline_haplo_headerMat = outline_haplo_headerMat +"\t"+ s +"\t"+ "LH"			
		fhaploidPat.write(outline_haplo_headerPat+"\n")
		fhaploidMat.write(outline_haplo_headerMat+"\n")
		flag_haplo_header = 0
	het1NumPat = 0
	het1NumMat = 0
	het2NumPat = 0
	het2NumMat = 0
	outline_haploPat = ""
	outline_haploMat = ""
	for i in range(5,len(col)):
		s=col[i].split("=")[0] #sample name
		if s == fathername:
			continue		
		elif s == mothername:
			continue
		elif (s in embryonames) or (s==childname):
			if dp[rs][s]<3:
				phase_to_husband="N"
				phase_to_wife="N"
				LH = 0				
			else:
				embryo = gt[rs][s] # embryo genotype
				if embryo == husband and embryo == wife :
					phase_to_husband="N"
					phase_to_wife="N"
					LH = 0
					
				elif (embryo[0] not in husband_bases) and (embryo[2] not in husband_bases):
					phase_to_husband="N"
					phase_to_wife="N"
					LH = 0
					
				elif (embryo[0] not in wife_bases) and (embryo[2] not in wife_bases):
					phase_to_husband="N"
					phase_to_wife="N"
					LH = 0
					
				elif (embryo[0] not in husband_bases) and (embryo[0] not in wife_bases):
					phase_to_husband="N"
					phase_to_wife="N"
					LH = 0
					
				elif (embryo[2] not in husband_bases) and (embryo[2] not in wife_bases):
					phase_to_husband="N"
					phase_to_wife="N"
					LH = 0
					
				elif (embryo[0] != embryo[2]): # embryo is het 
					if (embryo[0] in husband_bases):
						if (embryo[0] not in wife_bases):
							if embryo[2] in wife_bases:
								phase_to_husband = embryo[0]
								phase_to_wife = embryo[2]
								LH = 1
							else:
								#print line.strip()+"\tyan6"
								phase_to_husband="N"
								phase_to_wife="N"
								LH = 0
						elif (embryo[0] in wife_bases):
							if (embryo[2] in husband_bases):
								phase_to_husband = embryo[2]
								phase_to_wife = embryo[0]
								LH = 1
							elif (embryo[2] in wife_bases):
								phase_to_husband = embryo[0]
								phase_to_wife = embryo[2]
								LH = 1
							else:
								#print line.strip()+"\tyan7"
								phase_to_husband="N"
								phase_to_wife="N"
								LH = 0
						else:
							#print line.strip()+"\tyan8"
							phase_to_husband="N"
							phase_to_wife="N"
							LH = 0
					elif (embryo[0] in wife_bases): # embryo[0] is not in husband
						if (embryo[2] in husband_bases):
							phase_to_husband = embryo[2]
							phase_to_wife = embryo[0]
							LH = 1
						else:
							#print line.strip()+"\tyan9"
							phase_to_husband="N"
							phase_to_wife="N"
							LH = 0

					else:
						#print line.strip()+"\tyan10"
						phase_to_husband="N"
						phase_to_wife="N"
						LH = 0
				else: #if embryo is homo 
					phase_to_husband = embryo[0]
					phase_to_wife = embryo[0]
					LH = 1
			outline_haploPat = outline_haploPat + "\t" + phase_to_husband+"\t"+str(LH)
			outline_haploMat = outline_haploMat + "\t" + phase_to_wife+"\t"+str(LH)
			
		if phase_to_husband == husband_bases[0]:
			het1NumPat = het1NumPat + 1
		elif phase_to_husband  == husband_bases[1]:
			het2NumPat = het2NumPat + 1
		if phase_to_wife == wife_bases[0]:
			het1NumMat = het1NumMat + 1
		elif phase_to_wife  == wife_bases[1]:
			het2NumMat = het2NumMat + 1
	AllPat = het1NumPat+het2NumPat
	AllMat = het1NumMat+het2NumMat
	if (husband[0] != husband[2]) and (AllPat>0):
		fhaploidPat.write(chrom+"\t"+str(pos)+"\t"+husband_bases[0]+"\t"+str(1)+"\t"+str(het1NumPat)+"\t"+husband_bases[1]+"\t"+str(1)+"\t"+str(het2NumPat)+"\t"+str(AllPat)+"\t"+outline_haploPat.strip()+"\n")		
	if (wife[0] != wife[2]) and (AllMat>0):
		fhaploidMat.write(chrom+"\t"+str(pos)+"\t"+wife_bases[0]+"\t"+str(1)+"\t"+str(het1NumMat)+"\t"+wife_bases[1]+"\t"+str(1)+"\t"+str(het2NumMat)+"\t"+str(AllMat)+"\t"+outline_haploMat.strip()+"\n")		

	
fgt.close()
fhaploidPat.close()
fhaploidMat.close()
