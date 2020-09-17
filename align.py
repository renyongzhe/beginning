import sys
import numpy as np
import string

NUC44 = np.array([[5,-4,-4,-4,-2],
                  [-4,5,-4,-4,-2],
                  [-4,-4,5,-4,-2],
                  [-4,-4,-4,5,-2],
                  [-2,-2,-2,-2,-1]])
NBET='ATGCN'
def scoreMat(NUN44,NBET,seq1,seq2,gap=-8):
    len1,len2 = len(seq1),len(seq2)
    scoreMatrix = np.zeros((len1+1,len2+1),int)
    arrowMat = np.zeros((len1+1,len2+1),int)

    arrowMat[0,:] = np.ones(len2+1)
    arrowMat[1:,0] = np.zeros(len1)
    for i in range(1,len1+1):
        for j in range(1,len2+1):
            s_mat = np.zeros(4)
            s_mat[0] = scoreMatrix[i-1,j] + gap
            s_mat[1] = scoreMatrix[i,j-1] + gap
            n1,n2 = NBET.index(seq1[i-1]),NBET.index(seq2[j-1])
            s_mat[2] = scoreMatrix[i-1,j-1]+NUC44[n1,n2] 
            scoreMatrix[i,j] = np.max(s_mat)
            arrowMat[i,j] = np.argmax(s_mat)
    return scoreMatrix,arrowMat

def Align(scoreMatrix,arrow,seq1,seq2):
    aln_seq1,aln_seq2='',''
    flat_scorMat = np.ravel(scoreMatrix)
    v,h = divmod(np.argmax(flat_scorMat),len(seq2)+1)
    while True:
        if arrow[v,h] == 0:
            aln_seq1 += seq1[v-1]
            aln_seq2 += "-"
            v -= 1
        elif arrow[v,h] == 1:
            aln_seq1 += "-"
            aln_seq2 += seq2[h-1]
            h -= 1
        elif arrow[v,h] == 2:
            aln_seq1 += seq1[v-1]
            aln_seq2 += seq2[h-1]
            v -= 1
            h -= 1
        elif arrow[v,h] == 3:
            break
        if (v==0 and h==0) or scoreMatrix[v,h] == 0:
            break
    aln1 = aln_seq1[::-1]
    aln2 = aln_seq2[::-1]
    return aln1,aln2
sq1=input("Please input the sequence 1:")
sq2=input("Please input the sequence 2:")
scoreMatrix,arrowMatrix=scoreMat(NUC44,NBET,sq1,sq2,gap=-8)
aln1,aln2 = Align(scoreMatrix,arrowMatrix,sq1,sq2)
print('The score matrix is:')
print (scoreMatrix)
print ('The arrow matrix is:')
print (arrowMatrix)
print ('The aligned sequences are:')
print (aln1)
print (aln2)