#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import plotly.plotly as py
import sys






def open_file(position=0,file="test1.fasta"):
	print('file opened...')
	seqs=[]
	with open(file, "rU") as handle:
	    for record in SeqIO.parse(handle, "fasta"):
	       seqs.append(record.seq)
	    return seqs[position]


def show_h(freq):
	#print freq
	plt.hist(freq,30)
	plt.title("Tandem repeats distance histogram")
	plt.xlabel("repeat length")
	plt.ylabel("Frequency")
	
	plt.show()

	#fig = plt.gcf()


	# plt.bar(x, height= [1,2,3])
	# plt.xticks(x+.5, ['a','b','c']);

	# plot_url = py.plot_mpl(fig, filename='mpl-basic-histogram')


#sequence = 'ATACTAGCGGACTTGCTCGTGCTCATACTGTGCTCATACTGCCT'
#sequence = 'TCACAATAAGTAACCATTGTCACCGTTGTCACCGCAGTAGGTTCGAAC'
#sequence = 'CATGCGGCTACACAGTCGCGGCAACAGTCGCGGCAATCTACAA'
#sequence = 'TTCTCGCGGCTTTGCAATCAAGTCGCCGGTCAAGTCGCCGGTTGTCCACAA'
#sequence = 'TACTCGCGAGATTCGAGATTTCAAGCTTTGTCACCAATGGAA'
file="test1.fasta"
n_seq=0
sequence=''
#sequence ='GCTATTGGTTAGGGCGAGACCATCCTACATGTCCTCCAGTCATGCACACCAATCAACCTTTTTCGGTCACCCAACTCAGGGATTTAATCAAAAATAGTGGGTATCGACCGATCCCACTGAACTGGCCCTCACAGGTTACTCTTTCATGGTCAAGAGCAAAAAGAATGAGGGGTGTCGATTGATGATAAGTGGCGACCCCGACTTCGGCCATTCATAGAACGACTAATATGGCTGCTTGCACTAGGACATAAAGTCCGGCGGCTACGTTAAAGGAGTCGTAAAAGAATTGGGGGACGCACATGGGCTTACAGTCCGGCGGCTACGTTAAAGGAGTCGTAAAAGAATTGGGGGACGCACATGGGCTTACAAGAGCTATCCCCACCAATACCTCCTAGTGGTACCTGATGAGTAACACCGCTGCTATCGTC'
#sequence = 'AGCCCACCTATATACCACAGACAGCGGGGAGAGGCAAGCGGGGCAACGATGCCGAAATATCCATTTCCAAGGTTCGGCGCGATGCTTGAGGATAGTCCGCCAGACCCTCTCATTTACGAAATTCGCATTAGTGTACGTTCTAACTATGAAATGGTGATAAGGCTGTCACGTCTAGGGAGCGAGGGTGACAGTTCTTATTTTCCAGGCGCCGTAAATATAAACCGTGCTCAAGCAATCCCAAATCACACCGAGCTTCATAGGACGTATGTGTTGGCCGACTTTAACTTATTTTCCAGGCGCCGTAAATATAAACCGTGCTCAAGCAATCCCAAATCACACCGAGCTTCATAGGACGTATGTGTTGGCCGACTTTAAGTACGTGGGCATTATTTTGTCGCCGCCATCCTCTCCACTTTGCACTTTCCCTGCGCTTCGGTATGTTCAAGAAAACTAGACAAGTCGAAGCAGAAGTGCACCCGCGGCTTATTAACGGATACGTCTCATAGTTTCGAGTTCCTCGCGTCCTAGGCATCCCCTCCGGGACATCCGTCCCTTCTGTGAAGCCGCTCTGGGGTG'
#sequence = 'GTCGGATGTCTACGCTGGACCAAACCATACTAAGATATTCATTGCATTAGACCGCTATTGGCTTCAGCTTGCTCTCTCAAGCCAGCGTATAGCCTGCTATTGGCTTCAGCTTGCTCTCTCAAGCCAGCGTATAGCCTCACCTCAAGTTCGGTCGGTACACTTAAGAGGCGTCGTCTGTCGGGCACCTTCTGGGACGCTCCCTACTTCCATATAAATGACGTTGGACTACCGGGGTGTCTAACGCTCCTAACAATTAAGGCAGTTATCATGCAGCACTGCGTGATGGCGACGATGCCGATCACCTGACCCAGTCGGGAAATACAATTATTTGTCGTCTCAGTCCCCG'
#sequence = 'ATGATCCAAGGTTCTAGTGCTCGCGAGCTGAACGACGCCAGCGGGACGGAACGCAGGACTACTCCGCCAACACGGTCGTGCGTCATATGCTCGTATGCAACCCCTCAGCTTGGGCATCATCTATGACAGGACTACTCCGCCAACACGGTCGTGCGTCATATGCTCGTATGCAACCCCTCAGCTTGGGCATCATCTATGAAACTTCTCTGGACCCAGTCGACTAAGCAGCCAGCAAGGCTAGTAGAGCGCGTATAGTCTGTGTAGGTGCTTAACTTTCTGTACCGCGCGCGAAGCCCTTAGCGCACGGCCCAAGCGGTGTTTGGTTTATGACCCCGAAGCATGTGTACGAACTTATTTGTGAGTCGGGGCAGCTAGGTGATCGAGCTCTGGCTCGCGCAGACTATAATTGGACACATACATCTTGGATACCTGATTCCGCGATAACTACGTCCTAAACCCGTCTCTGTCCAGTAGGCAATGCGACCTGTC'
#sequence = 'CTTTATCCTAGAATAGGGTTTCTTGCAGCACAACAGTTTAGATCTAGAGCAGCTACTCTTCAATACATCGCGACAGATAATTTTTGTTACTAGTACGAAAAGACTATCCCCGAGAATCAATCGGCGGCACTAGTAGGTGTGTAGTGGGGTACCGGGACAGTAGACGTGCATCACAGTGGCCATGCGGAACTGACCAAATGCCAGGGGAATGCCTTAGAGCCAACGTAACGACAGACAATGCTGGTACCTGATGTGGTGAGGCATCGCACATCGAGCGCCATGTCTCGGTGGGATAAGAGGGATCACTGTCTGCTCATATGGACTTACCCGATAGCGACAATTTGCCATCACGGCATCGAGCGCCATGTCTCGGTGGGATAAGAGGGATCACTGTCTGCTCATATGGACTTACCCGATAGCGACAATTTGCCATCACGGAATGCACGTTGCACGCCGGCTATACTTTATTTCTAAATGTGTCGCTTACAGCTAACCTATGCGGCGTGACCATCGATTATCTCTTCATTATATTTCACACTTAAATTC'


#sequence = 'ATACACTAGCG'
#sequence='GTGGGACATACATAG'

def find_occurs(seq,word):
	cnt = seq.count(word)
	locs = []
	k = 0
	for i in range(cnt):
		k=seq.find(word,k)
		locs.append(k)
		k+=1
	return locs

def ex_words(seq,wsize):
	words = []
	L = len(seq)
	for i in range(L-wsize+1):
		words.append(seq[i:i+wsize])
	words = set(words)
	return list(words)

def distances(locs):
	N = len(locs)
	D = []
	for i in range(N):
		if(i<N-1):
			D.append(locs[i+1]-locs[i])
	D = np.asarray(D)
	return D 

def variat(words,start,size):
	for i in range(start,size+1):
		words=words+ex_words(sequence,i)
	return words


def histogram(upper,D,hist):
	nbins = upper 
	for i in range (len(D)):
		k=D[i]
		if(0 < k < nbins): 
			#hist[k] +=1
			hist.append(k)
	return hist

def iterative_levenshtein(s, t, costs=(1, 1, 1)):
    """ 
        iterative_levenshtein(s, t) -> ldist
        ldist is the Levenshtein distance between the strings 
        s and t.
        For all i and j, dist[i,j] will contain the Levenshtein 
        distance between the first i characters of s and the 
        first j characters of t
        
        costs: a tuple or a list with three integers (d, i, s)
               where d defines the costs for a deletion
                     i defines the costs for an insertion and
                     s defines the costs for a substitution
    """
    rows = len(s)+1
    cols = len(t)+1
    deletes, inserts, substitutes = costs
    
    dist = [[0 for x in range(cols)] for x in range(rows)]
    # source prefixes can be transformed into empty strings 
    # by deletions:
    for row in range(1, rows):
        dist[row][0] = row * deletes
    # target prefixes can be created from an empty source string
    # by inserting the characters
    for col in range(1, cols):
        dist[0][col] = col * inserts
        
    for col in range(1, cols):
        for row in range(1, rows):
            if s[row-1] == t[col-1]:
                cost = 0
            else:
                cost = substitutes
            dist[row][col] = min(dist[row-1][col] + deletes,
                                 dist[row][col-1] + inserts,
                                 dist[row-1][col-1] + cost) # substitution
    for r in range(rows):
        print(dist[r])
    
 
    return dist[row][col]


def verify_acuracy(distance,l1,l2):
	#check if the distance between the halves of the repeat is at most 0.1 times the minimum length of the two
	percentage = 0.1 * min([l1,l2])
	print ('percentage:',percentage)
	return (distance <= percentage) 





#main


def analyze(param1=2,param2=10):
	#hist= np.zeros(len(sequence),int)
	print('\n\n -----  starting analysis ----- \n\n')
	hist = [] 
	words=[]
	words=variat(words,param1,param2)
	print('sequence:\n\n',sequence,'\n\n')

	#print 'words:\n\n',words

	for i in words:
		locs=find_occurs(sequence,i)
		if(len(locs)>1):
			print (i,' locs:',locs)
			dist=distances(locs)
			hist=histogram(len(sequence),dist,hist)
			#hist=hist+h
			print (i,' distances of locs:',dist)
			if(dist[0]==len(i)):
				print ("-------->tandem repeat finded!")
				tandem = sequence[:locs[0]]+'-'+sequence[locs[0]:locs[1]]+'-'+sequence[locs[1]:] 
				print (sequence)
				print (tandem)
	show_h(hist)

def calcule_distance(s1='sitting',s2='sitting'):
	print ('word_a:',s1)
	print ('word_b:',s2)
	a = iterative_levenshtein(s1,s2)
	print ('levenshtein_distance:',a)
	#show_h(hist)
	print ('acuracy:',verify_acuracy(a,len(s1),len(s2)))


def get_seqs_repeats(seq,position,l1,l2):
	print(seq)
	seq1=seq[position:position+l1]
	seq2=seq[position+l1:position+l1+l2]
	recs=[SeqRecord(seq1,'repeat1', '', ''),SeqRecord(seq2,'repeat2', '', '')]
	print (recs)
	SeqIO.write(recs,"compare_sequences.fasta", "fasta")
	print('repeats saved!')
	


try:
	#gets minumum(a) and max(b) length of substrings to analyze
	file=str(sys.argv[1])
	n_seq=int(sys.argv[2])
	min_length=int(sys.argv[3])
	max_length=int(sys.argv[4])
	sequence=open_file(n_seq,file)
	analyze(min_length,max_length)
except:
	sequence = open_file(n_seq,file)
	analyze()
	#seq=open_file(1,'test3.fasta')
	#get_seqs_repeats(seq,4515,268,268)
#test_repeat()




#print 'words:',ex_words(sequence,2)
#print 'words:',ex_words(sequence,3)
#print 'words:',ex_words(sequence,4)
#print 'occurs:',find_occurs(sequence,'TA')

