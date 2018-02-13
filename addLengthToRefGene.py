#!/bin/python
import re
import sys

class Trans(object):
	def __init__(self,transcript):
		tmp = transcript.strip().split("\t")
		self.ID = tmp[1]
		self.gene = tmp[12]
		self.strand = tmp[3]
		self.first_exon_start = int(tmp[6])
		self.last_exon_end = int(tmp[7])
		self.numExons_origin= int(tmp[8])
		self.transExonStarts = map(int,tmp[9].split(",")[:-1])
		self.transExonEnds = map(int,tmp[10].split(",")[:-1])
		#self.cDNALength = GetLength()

	def GetCDSLength(self):
		exonStart_positions_NEW = []
		exonEnd_positions_NEW = []
		numExons = 0
		transLen = 0
		if self.first_exon_start != self.last_exon_end: # in this case no UTR assigned and replacement, else: pseudogene etc.
			for index in range(self.numExons_origin):
				# the following will trim UTR interval
				if self.transExonStarts[index] < self.last_exon_end and self.transExonEnds[index] > self.first_exon_start:
					numExons += 1

					if self.transExonStarts[index] < self.first_exon_start and self.transExonEnds[index] > self.first_exon_start and self.transExonEnds[index] <= self.last_exon_end:
						exonStart_positions_NEW.append(self.first_exon_start)
						exonEnd_positions_NEW.append(self.transExonEnds[index])
					elif self.transExonStarts[index] < self.last_exon_end and self.transExonEnds[index] > self.last_exon_end and self.transExonStarts[index] >= self.first_exon_start:
						exonStart_positions_NEW.append(self.transExonStarts[index])
						exonEnd_positions_NEW.append(self.last_exon_end)
					elif self.transExonEnds[index] > self.last_exon_end and self.transExonStarts[index] < self.first_exon_start:
						exonStart_positions_NEW.append(self.first_exon_start)
						exonEnd_positions_NEW.append(self.last_exon_end)
					else:
						exonStart_positions_NEW.append(self.transExonStarts[index])
						exonEnd_positions_NEW.append(self.transExonEnds[index])  
		else:
			exonStart_positions_NEW = self.transExonStarts
			exonEnd_positions_NEW = self.transExonEnds
				
		ExonInterval = map(lambda x:x[1]-x[0],zip(exonStart_positions_NEW,exonEnd_positions_NEW))
		transLen = sum(ExonInterval)
		return transLen

	def GetcDNAlength(self):
		ExonInterval = map(lambda x:x[1]-x[0],zip(self.transExonStarts,self.transExonEnds))
		transLen = sum(ExonInterval)
		return transLen

def GetNMLength(filein):
			transLenHash = {}
			transNumHash = {}
			geneStrand = {}
			with open(filein) as fi:
				for line in fi:
					if line.startswith("#"):
						print "\t".join([line.strip(),"cDNA","CDS"])
					else:
						trans = Trans(line)
						geneStrand[trans.gene] = trans.strand
						cDNALength = trans.GetcDNAlength()
						CDSLength = trans.GetCDSLength()
						print "\t".join([line.strip(),str(cDNALength),str(CDSLength)])

GetNMLength(sys.argv[1])