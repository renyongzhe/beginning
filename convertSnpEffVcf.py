#!/bin/pyton
# convert SnpEff vcf result to tab file

import os
import sys
import re
import argparse

parse = argparse.ArgumentParser(
    description="This script used for picking up the acceptable transcript and formating output.")

parse.add_argument('-c', '--cDNA', action='store_true', help="Count the transcript's cDNA length")
parse.add_argument('-C', '--CDS', action='store_true', help="Count the transcript's CDS length.[default]")
parse.add_argument('-f', '--first', action='store_true', help="Pick up the first ANN annotation.")
parse.add_argument('-o', '--out',  help="output_file.")
parse.add_argument('-i', '--inp',  help="input_file.")

args = parse.parse_args()


# define transcript class to choose the satisfactory transcript
class Transcripts(object):
    def __init__(self, transList):
        self.length = len(transList)
        self.choose = transList

    # self.NMs = self.GetNM()

    # from annotation get NM or NR id
    def GetTransID(self):
        # only retain NR_022014.1(HMGN2P46), NR_001566.1(TERC) in /fastdisk/pipeline-programs/SnpEff/snpEff/data/hg19/genes.refseq
        transNames = [re.search("\|(NM_.*?|NR_022014.*?|NR_001566.*?)\|", i).group(1) for i in self.choose if
                      re.search("\|(NM_.*?|NR_022014.*?|NR_001566.*?)\|", i)]
        # print transNames
        return transNames

    # filter harmful class
    def HarmfulClass(self):
        transList_out = []
        # HIGH > MODERATE > LOW > MODIFIER,give priority to output by this order
        High_list = [high for high in self.choose if re.search("HIGH", high)]
        if len(High_list) != 0:
            self.choose = High_list
        else:
            Moderate_list = [moderate for moderate in self.choose if re.search("MODERATE", moderate)]
            if len(Moderate_list) != 0:
                self.choose = Moderate_list
            else:
                Low_list = [low for low in self.choose if re.search("LOW", low)]
                if len(Low_list) != 0:
                    self.choose = Low_list
                else:
                    Modifier_list = [modifier for modifier in self.choose if re.search("MODIFIER", modifier)]
                    self.choose = Modifier_list

    # if ANN included NM and NR, retain NM preferentially
    def GetNM(self):
        NMs = [i for i in self.choose if re.search("NM_", i)]
        if len(NMs) != 0:
            self.choose = NMs

    # count the transcript length and pick up the max
    def ChooseMaxLength(self, transHash):
        transLength = []
        transList_out = []
        transNMs = self.GetTransID()
        # print transHash
        if len(transNMs) == 0:
            self.choose = []
        # print "---",transNMs
        else:
            # print transNMs
            for transID_pre in transNMs:
                transID = transID_pre.split(".")[0] if len(transID_pre.split(".")) > 1 else transID_pre

                # count the length of transcript
                # print transID,"----"
                if transID in transHash:
                    # print transID,"******"
                    # Default the transID includ in transHash
                    transLength.append(transHash[transID])
            # print "error:::",transLength
            transLengthMax = max(transLength)
            transLengthMaxCount = transLength.count(transLengthMax)
            for index in range(0, transLengthMaxCount):
                tmpIndex = transLength.index(transLengthMax)
                transList_out.append(self.choose[tmpIndex])
                transLength[tmpIndex] = transLengthMax - 1
            self.choose = transList_out

    # pick up the minimum transcript ID
    def ChooseMiniNM(self):
        transNMs = self.GetTransID()
        transMiniNum = map(lambda x: float(x.split("_")[1]), transNMs)
        transMiniNumPos = transMiniNum.index(min(transMiniNum))
        # print "before : ",self.choose
        # print "mimnum : ",transMiniNum
        # print "mimnum index : ",transMiniNumPos
        self.choose = [self.choose[transMiniNumPos]]

    # print "after: ",self.choose


# count the transcript length from reference file
class Trans(object):
    def __init__(self, transcript):
        tmp = transcript.strip().split("\t")
        # NM_001005221.2
        self.ID = tmp[1]
        self.gene = tmp[12]
        self.strand = tmp[3]
        self.cDNA_Start = tmp[4]
        self.first_exon_start = int(tmp[6])
        self.last_exon_end = int(tmp[7])
        self.numExons_origin = int(tmp[8])
        self.transExonStarts = map(int, tmp[9].split(",")[:-1])
        self.transExonEnds = map(int, tmp[10].split(",")[:-1])

    # self.cDNALength = GetLength()

    def GetCDSLength(self):
        exonStart_positions_NEW = []
        exonEnd_positions_NEW = []
        numExons = 0
        transLen = 0
        if self.first_exon_start != self.last_exon_end:  # in this case no UTR assigned and replacement, else: pseudogene etc.
            for index in range(self.numExons_origin):
                # the following will trim UTR interval
                if self.transExonStarts[index] < self.last_exon_end and self.transExonEnds[
                    index] > self.first_exon_start:
                    numExons += 1

                    if self.transExonStarts[index] < self.first_exon_start and self.transExonEnds[
                        index] > self.first_exon_start and self.transExonEnds[index] <= self.last_exon_end:
                        exonStart_positions_NEW.append(self.first_exon_start)
                        exonEnd_positions_NEW.append(self.transExonEnds[index])
                    elif self.transExonStarts[index] < self.last_exon_end and self.transExonEnds[
                        index] > self.last_exon_end and self.transExonStarts[index] >= self.first_exon_start:
                        exonStart_positions_NEW.append(self.transExonStarts[index])
                        exonEnd_positions_NEW.append(self.last_exon_end)
                    elif self.transExonEnds[index] > self.last_exon_end and self.transExonStarts[
                        index] < self.first_exon_start:
                        exonStart_positions_NEW.append(self.first_exon_start)
                        exonEnd_positions_NEW.append(self.last_exon_end)
                    else:
                        exonStart_positions_NEW.append(self.transExonStarts[index])
                        exonEnd_positions_NEW.append(self.transExonEnds[index])
        else:
            exonStart_positions_NEW = self.transExonStarts
            exonEnd_positions_NEW = self.transExonEnds

        ExonInterval = map(lambda x: x[1] - x[0], zip(exonStart_positions_NEW, exonEnd_positions_NEW))
        transLen = sum(ExonInterval)
        return transLen

    def GetcDNAlength(self):
        ExonInterval = map(lambda x: x[1] - x[0], zip(self.transExonStarts, self.transExonEnds))
        transLen = sum(ExonInterval)
        return transLen


def GetNMLength():
    transLenHash = {}
    transNumHash = {}
    geneStrand = {}
    nmId_cDNAStart = {}
    #with open("/mnt/pipeline-programs/SnpEff/snpEff/data/hg19/genes.refseq") as fi:
    for line in open("/mnt/pipeline-programs/SnpEff/snpEff/data/hg19/genes.refseq"):
        if line.startswith("#"):
            continue
        trans = Trans(line)
        # geneStrand include strand information and cDNA start
        geneStrand[trans.gene] = trans.strand

        # convert NM_004448.3 to NM_004448, remove version number
        trans_ID = trans.ID.split(".")[0] if len(trans.ID.split(".")) > 1 else trans.ID
        nmId_cDNAStart[trans_ID] = trans.cDNA_Start
        # print nmId_cDNAStart
        if args.cDNA:
            transLength = trans.GetcDNAlength()
        else:
            transLength = trans.GetCDSLength()
        # print transID,transLength
        if trans_ID not in transLenHash:
            transLenHash[trans_ID] = transLength
            transNumHash[trans_ID] = 1
        else:
            transNumHash[trans_ID] += 1
            if transLength > transLenHash[trans_ID]:
                transLenHash[trans_ID] = transLength
    # for j in transLenHash:
    # 	print "\t".join([j,str(transLenHash[j])])
    return transLenHash, geneStrand, nmId_cDNAStart


def ConvertAA():
    three = "Gly Ala Val Leu Ile Phe Trp Tyr Asp Asn Glu Lys Gln Met Ser Thr Cys Pro His Arg"
    one = "G A V L I F W Y D N E K Q M S T C P H R"
    three_3 = three.split(" ")
    one_1 = one.split(" ")
    AAHash = dict(zip(three_3, one_1))
    return AAHash


def convert(filein):
    transHash, geneStrand, nmId_cDNAStart = GetNMLength()
 
    AAHash = ConvertAA()
    NonHarmNum = 0
    AllItems = 0
    NRsOrMultiNMs = 0
    transEqualLength = 0
    fileout = args.out
    with open(fileout, 'w') as fo:
        fo.write(
            "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	Depth_Tumor	Depth_Normal	AF_Tumor	AF_Normal	Mut_Reads_Tumor	Mut_Reads_Normal	FR_Score	DP5_Tumor	Gene	TranscriptID	Strand	Type_1	Type_2	cDNA_Change	Protein_Change	Exon_Number	cDNA_Position	cDNA_Start" + "\n")
        with open(filein) as fi:
            for line in fi:
                if not re.match("[#\"]", line):
                    AllItems += 1
                    tmp = line.strip().split("\t")
                    # get Type 2
                    if len(tmp[3]) == len(tmp[4]) == 2:
                        Type_2 = "DNM"
                    elif len(tmp[3]) == len(tmp[4]) == 3:
                        Type_2 = "TNM"
                    elif (len(tmp[3]) == len(tmp[4])) >= 4:
                        Type_2 = "ONM"
                    elif len(tmp[3]) > len(tmp[4]) and len(tmp[4]) == 1:
                        Type_2 = "DEL"
                    elif len(tmp[3]) < len(tmp[4]) and len(tmp[3]) == 1:
                        Type_2 = "INS"
                    else:
                        Type_2 = "SNM"

                    # get tag information

                    info = tmp[7].split(";")
                    # print tmp
                    for var in info:
                        if re.match("DP=", var):
                            Depth_T = var.split("=")[1]
                        elif re.match("DP_N=", var):
                            Depth_N = var.split("=")[1]
                        elif re.match("AF=", var):
                            AF_T = var.split("=")[1]
                        elif re.match("AF_N=", var):
                            AF_N = var.split("=")[1]
                        elif re.match("FR=", var):
                            FR_Score = var.split("=")[1]
                        elif re.match("DP5=", var):
                            DP5 = var.split("=")[1]
                    # print Depth_T,AF_T
                    Mut_Reads_T = round(float(Depth_T) * float(AF_T))
                    Mut_Reads_N = round(float(Depth_N) * float(AF_N))

                    # all transcripts will be NonHarmNum after filtering following rules
                    transcripts = [item for item in tmp[7].split(";") if re.match("ANN", item)][0].split("=")[1].split(
                        ",")

                    # chose the first ANN annotation without regard to other all rules
                    if args.first:
                        OutVirs = transcripts[0].split("|")
                    else:
                        # our rules for picking up transcript
                        trans = Transcripts(transcripts)
                        # trans.HarmfulClass()
                        trans.GetNM()
                        if len(trans.choose) == 1:  # only one satisfies required class
                            OutVirs = trans.choose[0].split("|")

                        else:
                            # trans.GetNM()
                            trans.HarmfulClass()
                            if len(trans.choose) == 1:  # only one is NM,other is NR
                                OutVirs = trans.choose[0].split("|")
                            else:
                                NRsOrMultiNMs += 1
                                # print transcripts
                                trans.ChooseMaxLength(transHash)
                                # print "choose",trans.choose
                                if len(trans.choose) == 0:
                                    continue
                                elif len(trans.choose) == 1:  # only one has the maximun transcript length
                                    OutVirs = trans.choose[0].split("|")
                                else:
                                    transEqualLength += 1
                                    trans.ChooseMiniNM()  # get minimum NM/NR ID, after all this steps,if still have same ID or equal length, catch the first
                                    OutVirs = trans.choose[0].split("|")

                    ExonNum = OutVirs[8].split("/")[0]
                    ExonNum = "-" if ExonNum == "" else ExonNum
                    Gene = OutVirs[3]
                    Strand = geneStrand[Gene] if Gene in geneStrand else "-"  # some special gene name is inconsistent
                    
                    feature_type = OutVirs[5]
                    TranscriptID = OutVirs[6] if re.search("transcript",feature_type) else "-"

               
                    cDNA_Position = OutVirs[11].split("/")[0]
                    cDNA_Position = "-" if cDNA_Position == "" else cDNA_Position
                    # Type_1 = OutVirs[1].split("&")[0]
                    Type_1 = OutVirs[1]
                    Type_1=Type_1.replace("_variant","")
                    if TranscriptID in nmId_cDNAStart:
                        cDNAStart = nmId_cDNAStart[TranscriptID]
                    else:
                        #print line
                        cDNAStart="0"
                        #ExonNum="0"
                        #cDNA_Position="0"
                        #raise KeyError("%s do not in genes.refseq" % TranscriptID)

                    # convert amino acid
                    if len(OutVirs[10]) != 0:
                    	AAchange = re.sub("[A-Z][a-z]{2}",lambda x:AAHash[x.group()],OutVirs[10])

                    else:
                        AAchange = "-"
                    OutVir = "\t".join(
                        [Depth_T, Depth_N, AF_T, AF_N, str(Mut_Reads_T), str(Mut_Reads_N), FR_Score, DP5, Gene,
                         TranscriptID, Strand, Type_1, Type_2, OutVirs[9], AAchange, ExonNum, cDNA_Position, cDNAStart])
                    Outmp = "\t".join(tmp[0:7])
                    # print "\t".join([Outmp,OutVir])
                    fo.write("\t".join([Outmp, OutVir]) + "\n")

        print "All items :  %s" % AllItems
        print "NRsOrMultiNMs items : %s " % (NRsOrMultiNMs)
        print "Equal transcript length items : %s" % transEqualLength
        print "Result saved in %s" % (fileout)


if __name__ == "__main__":
    convert(args.inp)
