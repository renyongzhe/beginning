#!/bin/pyton
# convert SnpEff vcf result to tab file

import os
import sys
import re

#define transcript class to choose the satisfactory transcript
class Transcripts(object):
	def __init__(self,transList):
		self.length = len(transList)
		self.data = transList

	# filter harmful class
	def HarmfulClass(self):
		transList_out = []
		for i in self.data:
			if re.search("HIGH|MODERATE",i):
				transList_out.append(i)
		return transList_out

	# from annotation get NM or NR id
	def GetNM(self):
		transNames = [ re.search("\|(N[M|R]_.*?)\|",i).group(1) for i in self.data if re.search("N[M|R]_",i) ]
		return transNames

	# count the transcript length and pick up the max
	def ChooseMaxLength(self,transHash):
		transLength = []
		transList_out = []
		transNames = self.GetNM()

		for transID_pre in transNames:
			transID = transID_pre.split(".")[0] if len(transID_pre.split(".")) >1 else transID_pre
	# count the length of transcript
			if transID in transHash:
			# Default the transID includ in transHash	
				transLength.append(transHash[transID])
		transLengthMax = max(transLength)
		transLengthMaxCount = transLength.count(transLengthMax)
		for index in range(0, transLengthMaxCount):
			tmpIndex = transLength.index(transLengthMax)
			transList_out.append(self.data[tmpIndex])
			transLength[tmpIndex] =  transLengthMax - 1
		return transList_out

	# pick up the minimum transcript ID
	def ChooseMiniNM(self):
		transNames = self.GetNM()
		transMiniNum = map(lambda x:float(x.split("_")[1]),transNames)
		transMiniNumPos = transMiniNum.index(min(transMiniNum))
		return self.data[transMiniNumPos]

# count the transcript length from reference file
def GetNMLength():
			transLenHash = {}
			transNumHash = {}
			with open("/fastdisk/pipeline-programs/SnpEff/snpEff/data/hg19/genes.refseq") as fi:
				for line in fi:
					if line.startswith("#"):
						continue
					tmp = line.strip().split("\t")
					transID = tmp[1]
					first_exon_start = int(tmp[6])
					last_exon_end = int(tmp[7])
					numExons_origin= int(tmp[8])
					transExonStarts = map(int,tmp[9].split(",")[:-1])
					transExonEnds = map(int,tmp[10].split(",")[:-1])
					exonStart_positions_NEW = []
					exonEnd_positions_NEW = []
					numExons = 0

					# the following will trim UTR interval
					if first_exon_start != last_exon_end: # in this case no UTR assigned and replacement, else: pseudogene etc.
						for index in range(numExons_origin):
							if transExonStarts[index] < last_exon_end and transExonEnds[index] > first_exon_start:
								numExons += 1
								# replace first/last exons:
								if transExonStarts[index] < first_exon_start and transExonEnds[index] > first_exon_start and transExonEnds[index] <= last_exon_end:
									exonStart_positions_NEW.append(first_exon_start)
									exonEnd_positions_NEW.append(transExonEnds[index])
								elif transExonStarts[index] < last_exon_end and transExonEnds[index] > last_exon_end and transExonStarts[index] >= first_exon_start:
									exonStart_positions_NEW.append(transExonStarts[index])
									exonEnd_positions_NEW.append(last_exon_end)
								elif transExonEnds[index] > last_exon_end and transExonStarts[index] < first_exon_start:
									exonStart_positions_NEW.append(first_exon_start)
									exonEnd_positions_NEW.append(last_exon_end)
								else:
									exonStart_positions_NEW.append(transExonStarts[index])
									exonEnd_positions_NEW.append(transExonEnds[index])  
					else:
						exonStart_positions_NEW = transExonStarts
						exonEnd_positions_NEW = transExonEnds

					ExonInterval = map(lambda x:x[1]-x[0],zip(exonStart_positions_NEW,exonEnd_positions_NEW))
					transLength = sum(ExonInterval)
					#print transID,transLength
					if transID not in transLenHash:
						transLenHash[transID] = transLength
						transNumHash[transID] = 1
					else:
						transNumHash[transID] += 1 
						if transLength > transLenHash[transID]:
							transLenHash[transID] = transLength
			for j in transNumHash:
				print "\t".join([j,str(transNumHash[j])])
			return transLenHash

def ConvertAA():
	three = "Gly Ala Val Leu Ile Phe Trp Tyr Asp Asn Glu Lys Gln Met Ser Thr Cys Pro His Arg"
	one = "G A V L I F W Y D N E K Q M S T C P H R"
	three_3 = three.split(" ")
	one_1 = one.split(" ")
	AAHash = dict(zip(three_3,one_1))
	return AAHash

def convert(filein):
	transHash = GetNMLength()
	AAHash = ConvertAA()
	NonHarmNum = 0
	AllItems = 0 
	transEqualLength = 0 
	fileout = filein + ".filter"
	with open(fileout,'w') as fo:
		with open(filein) as fi:
			for line in fi:
				if not re.match("#",line):
					AllItems += 1 
					tmp = line.strip().split("\t")
					transcripts = [ item for item in tmp[7].split(";") if re.match("ANN",item) ][0].split("=")[1].split(",")
					# all transcripts will be NonHarmNum after filtering following rules
					trans = Transcripts(transcripts)
					transHarm = trans.HarmfulClass()
			 		if len(transHarm) == 0:
						NonHarmNum += 1
					else:
						transHarmMaxLength = Transcripts(transHarm)
						transHarmMaxLength = transHarmMaxLength.ChooseMaxLength(transHash)
						if len(transHarmMaxLength) == 1:
							OutVirs = transHarmMaxLength[0].split("|")
						else:
							transEqualLength += 1
							transHarmMaxLengthMiniNM = Transcripts(transHarmMaxLength)
							OutVirs = transHarmMaxLengthMiniNM.ChooseMiniNM().split("|")

						ExonNum = OutVirs[8].split("/")[0]

						if len(OutVirs[10]) != 0:
							aa = re.search("p.([A-Z][a-z]{2})([0-9]*)(.*)",OutVirs[10])
							OriginAA = AAHash[aa.group(1)] if aa.group(1) in AAHash else aa.group(1)
							AltAA = AAHash[aa.group(3)] if aa.group(3) in AAHash else aa.group(3)
							AAchange = "p."+OriginAA+aa.group(2)+AltAA
						else:
							AAchange = "."
						OutVir = "\t".join([OutVirs[1],OutVirs[3],OutVirs[6],ExonNum,OutVirs[9],AAchange])
						Outmp = "\t".join(tmp[0:7])
						#print "\t".join([Outmp,OutVir])
						fo.write("\t".join([Outmp,OutVir])+"\n")
				# else:
				# 	fo.write(line.strip()+"\n")
		print "All items :  %s"%AllItems
		print "NonHarmNum items : %s " %(NonHarmNum)
		print "Equal transcript length items : %s" %transEqualLength
		print "Result saved in %s"%(fileout)

convert(sys.argv[1])