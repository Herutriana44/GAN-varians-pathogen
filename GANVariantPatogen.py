import random as rd
import time as t
from decimal import Decimal

def randomProtein(maxLen):
	amino_acid_codes = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	res = ""
	length = rd.randint(2, maxLen)
	for i in range(length):
		randomOfIndex = rd.randint(0,19)
		res += amino_acid_codes[randomOfIndex]
	return res.lower()

def global_alignment(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
	# Create a matrix to store the scores
	matrix = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
	
	# Initialize the first row and column of the matrix
	for i in range(1, len(seq1)+1):
		matrix[i][0] = matrix[i-1][0] + gap_penalty
	for j in range(1, len(seq2)+1):
		matrix[0][j] = matrix[0][j-1] + gap_penalty
	
	# Fill in the rest of the matrix
	for i in range(1, len(seq1)+1):
		for j in range(1, len(seq2)+1):
			# Calculate the score for the three possible alignments
			match = matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
			delete = matrix[i-1][j] + gap_penalty
			insert = matrix[i][j-1] + gap_penalty
			
			# Choose the highest scoring alignment
			matrix[i][j] = max(match, delete, insert)
	
	# Trace back through the matrix to find the optimal alignment
	align1, align2 = '', ''
	i, j = len(seq1), len(seq2)
	while i > 0 and j > 0:
		score = matrix[i][j]
		diag_score = matrix[i-1][j-1]
		up_score = matrix[i][j-1]
		left_score = matrix[i-1][j]
		
		# Check which cell the score came from
		if score == diag_score + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
			align1 += seq1[i-1]
			align2 += seq2[j-1]
			i -= 1
			j -= 1
		elif score == up_score + gap_penalty:
			align1 += '-'
			align2 += seq2[j-1]
			j -= 1
		elif score == left_score + gap_penalty:
			align1 += seq1[i-1]
			align2 += '-'
			i -= 1
	
	# Finish tracing if necessary
	while i > 0:
		align1 += seq1[i-1]
		align2 += '-'
		i -= 1
	while j > 0:
		align1 += '-'
		align2 += seq2[j-1]
		j -= 1
	
	# Return the optimal alignment
	return align1[::-1], reverse_string(align2)

def reverse_string(s):
    reversed_string = ""
    for i in range(len(s) - 1, -1, -1):
        reversed_string += s[i]
    return reversed_string

def calculate_homology(align1, align2):
	# Count the number of matches
	matches = sum(1 for a, b in zip(align1, align2) if a == b)
	
	# Calculate the homology
	homology = matches / len(align1)
	
	return homology


def calculate_similarity(align1, align2):
	# Count the number of matches and gaps
	matches = sum(1 for a, b in zip(align1, align2) if a == b)
	gaps = sum(1 for a, b in zip(align1, align2) if a == '-' or b == '-')
	
	# Calculate the similarity
	similarity = (matches - gaps) / len(align1)
	
	return similarity
	
	
# GAN Section
	
def Generator(maxLength=30):
	prot = randomProtein(maxLength)
	return prot

def Discriminator(inp,seqTemplate, tHomology,tSimilar):
	align1, align2 = global_alignment(inp, seqTemplate)
	homology = calculate_homology(align1, align2)
	similarity = calculate_similarity(align1, align2)
	if((homology >= tHomology) and (similarity >= tSimilar)):
		return True
	else:
		return False
	
def run(inputSeq,tHomology,tSimilarity, maxSeq=10,maxLength=30):
	tBegin = t.time()
	res = []
	while(len(res) <= maxSeq):
		prot = Generator(maxLength)
		#print(f'protein generator : {prot}', end=" ")
		discriminator = Discriminator(prot, inputSeq, tHomology, tSimilarity)
		if(discriminator):
			#print(f'=> discriminator : lolos')
			res.append(prot)
		else:
			while(True):
				prot = revice(prot, inputSeq)
				discriminator = Discriminator(prot, inputSeq, tHomology, tSimilarity)
				if(discriminator):
					res.append(prot)
					break
		
		if(len(res) > maxSeq):
			break
	
	tEnd = t.time()
	waktu = tEnd - tBegin
	return res, Decimal(waktu)
	
	
def revice(inp, seqTemplate):
	align1, align2 = global_alignment(inp, seqTemplate)
	amino_acid_codes = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	
	res = ""
	for i in range(0,len(align1)):
		if(align1[i] == align2[i]):
			res += align1[i]
		if(align1[i] != align2[i]):
			indRand = rd.randint(0,19)
			while(align1[i] == amino_acid_codes[indRand]):
				indRand = rd.randint(0,19)
			res += amino_acid_codes[indRand]
	return res.lower()
	
	
	
input = "earrfnqyyssikrsgsiq"
conTemp =  "earlflqyygsjkmq"
"""
print(input)
print(conTemp)
print(revice(input,conTemp))
"""
#print(run(input.lower(), 0.7, 0.1, 1, len(input)))


#print(input)
#print(run(input.lower(), 0.7, 0.1, 1, len(input)))

#print(Decimal(0.0000000000001))

import pandas as pd

df = pd.read_csv('time.csv')
df['result'] = ""
print(df)

for i in range(0,len(df)):
    print(i)
    df['result'][i], df['time'][i] = run(df['protein'][i].lower(), df['tHomology'][i], df['tSimilarity'][i], 10, df['panjang'][i])

df.to_csv('timeRes.csv', index=False)

#print(randomProtein(10))
