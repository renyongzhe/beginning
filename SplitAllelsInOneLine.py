#!/bin/python 
# This script used for split multiple alleles one line into multiple lines
# input vcf or vcf.gz

import sys
import re
import gzip

def GetInfo(info,index):
# get INFO by every allel
	#infoHash = {}
	infoOut = []
	for i in info.strip().split(";"):
		# some field is just a tag,do not have value,like STR 
		if re.search("=",i):
			KeyValue = i.split("=")
			Key = KeyValue[0]
			#print KeyValue
			Value = KeyValue[1].split(",")
			#infoHash[KeyValue[0]] = KeyValue[1].split(",")
			if len(Value) >1 :
				Field = Key + "=" + Value[index]
				infoOut.append(Field)
			else:
				infoOut.append(i)
		else:
			infoOut.append(i)

	return ";".join(infoOut)

def GetFormat(format,allelesNum,index,type):
	# get format information by every allel
	formatout = []
	tmp = format.strip().split(":")
	Genotypes = tmp[0]
	Genotype = "0/1" if re.match("0",Genotypes) else "1/1"
	VariableToSplitNum = len(tmp) -2 
	if type == "normal":
		Genotype ="0/0"
		VariableToSplitNum = len(tmp)
	formatout.append(Genotype)
	# the last two field is SA_MAP_AF ,SA_POST_PROB ,each has 3 values, all retained in every allel,
	# only in tumor format

	for i in tmp[1:VariableToSplitNum]:
		Values = i.split(",")
		if len(Values) > allelesNum:
			Value = Values[0] + "," + Values[index+1]
			formatout.append(Value)
		else:
			#print len(Values),index
			Value = Values[index]
			formatout.append(Value)
	if type == "tumor":
	# must add the last 2 field iin tumor
		formatout.extend(tmp[-2:])
	return ":".join(formatout)

def GetOutput(tmp,ALTSNum,ALTOrder,POS,REF,ALT):
	INFO = GetInfo(tmp[7],ALTOrder)
	TUMOR_FORMAT = GetFormat(tmp[9],ALTSNum,ALTOrder,"tumor")
	NORMAL_FORMAT = GetFormat(tmp[10],ALTSNum,ALTOrder,"normal")
	return "\t".join([tmp[0],str(POS),tmp[2],REF,ALT,tmp[5],tmp[6],INFO,tmp[8],TUMOR_FORMAT,NORMAL_FORMAT])

def formatVcf(filein):
	if filein.endswith("gz"):
		fi = gzip.open(filein)
	else:
		fi = open(filein)
	for line in fi:
		if re.match("#",line):
			print line.strip()
		else:
			tmp = line.strip().split("\t")
			if re.search(",",tmp[4]):
				REF = tmp[3]
				ALTS = tmp[4].split(",")
				ALTSNum = len(ALTS)
				POS = 0

				#print line
				for ALT in ALTS:
					REF = tmp[3]
					ALTOrder = ALTS.index(ALT)
				# deletion
					if len(REF) > len(ALT):
						if len(ALT) == 1:
							POS = tmp[1]
							print GetOutput(tmp,ALTSNum,ALTOrder,POS,REF,ALT)
						else:
							ALTLen = len(ALT)
							ALT = ALT[-1]
							subtractLen = ALTLen - len(ALT)
							REF = REF[subtractLen:]
							POS = int(tmp[1]) + subtractLen
							print GetOutput(tmp,ALTSNum,ALTOrder,POS,REF,ALT)
					
					# insertion
					elif len(REF) < len(ALT):
						REFLen = len(REF)
						REF = REF[-1]
						subtractLen = REFLen - len(REF)
						ALT = ALT[subtractLen:]
						POS = int(tmp[1]) + subtractLen
						print GetOutput(tmp,ALTSNum,ALTOrder,POS,REF,ALT)
					# SNP
					elif len(REF) == len(ALT):
						for index in range(0,len(REF)):
							if not REF[index] == ALT[index]:
								POS = int(tmp[1]) + index 
								REF = REF[index]
								ALT = ALT[index]
								print GetOutput(tmp,ALTSNum,ALTOrder,POS,REF,ALT)
								break
			else:
				print line.strip()
	fi.close()

formatVcf(sys.argv[1])