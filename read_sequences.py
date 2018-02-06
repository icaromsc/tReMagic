from __future__ import print_function
from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment
from operator import itemgetter, attrgetter




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
    #for r in range(rows):
        #print(dist[r])
    
 
    return dist[row][col]


def verify_acuracy(distance,l1,l2):
	#check if the distance between the halves of the repeat is at most 0.1 times the minimum length of the two
	percentage = 0.1 * min([l1,l2])
	print ('percentage:',percentage)
	print("is accurate?")
	return (distance <= percentage) 



def extract_words(seq,wsize,start=0,cut=1):
	words = []
	L = len(seq)
	for i in range(start,L-wsize+1,cut):
		words.append(seq[i:i+wsize])
	words = set(words)
	return list(words)


def extract_words_positions(seq,wsize,start=0,cut=1):
	pos = []
	pos.append(cut)
	L = len(seq)
	for i in range(start,L-wsize+1,cut):
		temp=i+wsize
		if temp <= L:
			pos.append(temp)
		else:
			break
	return pos


def get_distances(seq,wsize,positions):
	#candidates=[]
	distances=[]
	for position in positions:
		seqs=get_seqs_repeats(seq,position,wsize,wsize)
		distance=iterative_levenshtein(seqs[0],seqs[1])
		#candidate=(position,position+wsize,distance)
		#candidates.append(candidate)
		distances.append(distance)

	#print candidates
	return candidates

def get_candidates(positions,wsize,distances):
	candidates=[]
	for position,distance in zip(positions,distances):
		candidate=(position,wsize,distance)
		print(candidate)
		candidates.append(candidate)

	print("\nSorting candidates...\n")
	candidates = sorted(candidates, key=itemgetter(2))
	print("\n getting 10 best candidates...\n")
	for i in range(10):
		print(candidates[i])
	




def get_best_candidate(distances):
	index_min = min(xrange(len(distances)), key=distances.__getitem__)
	return index_min 

# seqs=[]
# with open("test3.fasta", "rU") as handle:
#     for record in SeqIO.parse(handle, "fasta"):
#        seqs.append(record.seq)

# log = open("out_test.txt", "w")
# print('first sequence:\n',seqs[0], file = log)


def open_file(position=0,file="test3.fasta"):
	print('file opened...')
	seqs=[]
	with open(file, "rU") as handle:
	    for record in SeqIO.parse(handle, "fasta"):
	       seqs.append(record.seq)
	    return seqs[position]


def get_seqs_repeats(seq,position,l1,l2,saveFile=False):
	#print(seq)
	seq1=seq[position:position+l1]
	seq2=seq[position+l1:position+l1+l2]
	#print (recs)
	if saveFile:
		recs=[SeqRecord(seq1,'repeat1', '', ''),SeqRecord(seq2,'repeat2', '', '')]
		seq_name=str("repeats_pos_"+str(position)+"_len_"+str(l1)+".fasta")
		SeqIO.write(recs,seq_name, "fasta")
		print('repeats saved!')
	#align(seq1,seq2)
	return [seq1,seq2]
	
def align(seq1,seq2):
	alignments = pairwise2.align.globalxx(seq1, seq2)
	ali=format_alignment(*alignments[0])
	print (ali)
	return ali 

start=2361
l1,l2=480,480
sequence=open_file()
seq=get_seqs_repeats(sequence,start,l1,l2)
align(seq[0],seq[1])
d=iterative_levenshtein(seq[0],seq[1])
print('distance:',d)
print (verify_acuracy(d,l1,l2))
