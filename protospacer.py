#--------------------------------------------------------------------------------------------------------------------------------------------------
'''
 a function that identifies all candidate guide (protospacer) sequences from a given genomic region.
'''
#--------------------------------------------------------------------------------------------------------------------------------------------------
def candidate_sequence(dna_seq:str, ps_length:int, pam:str):
    
    length_of_dna_seq = len(dna_seq)
    pam_length = 0  #variable to hold the length of PAM depending on NCC or NGG
    pam_str = ""    #variable to hold the PAM string from the sequence
    candidates = [] #result list of the protospacers
    candidate = ""  #variable to hold a candidate
    #--------------------------------------------------------
    #check if the PAM is NCC or NGG. pam_length is 3 if pam is NGG, else it's 5
    if pam is "NGG":
        pam_length = 3
    elif pam is "NCC":
        pam_length = 5
    #--------------------------------------------------------
    #get the pam_str depending on the pam_length
    if pam_length == 3:
        pam_str = dna_seq[-3:]
    else:
        pam_str = dna_seq[-5:]
    rev_pam = pam_str[::-1] #get the reverse of the pam_str to search in the input DNA sequence
    #--------------------------------------------------------
    #check if the rev_pam is actually in the dna_seq. 
    #If not present then simply return empty list 
    if rev_pam in dna_seq:
        index = dna_seq.index(rev_pam)
    else:
        return candidates 
    #--------------------------------------------------------
    #find candidate from the dna_seq and add to the result candidates list
    for seq in dna_seq:
        diff = length_of_dna_seq - index
        if diff < ps_length:
            continue
        if pam_length > 0 :
            end = index + ps_length
            candidate = dna_seq[index:end]
            candidates.append(candidate)
            index += 1
            pam_length -= 1
    return candidates

#-------------------------------------------------
#       sample test
#-------------------------------------------------
s = "TGATCTACTAGAGACTACTAACGGGGATACATAG"
l = 20
p = "NGG"
print(candidate_sequence(s,l,p))
#--------------------------------------------------------------------------------------------------------------------------------------------------
'''
take a reference genome FASTA file 
'''
#--------------------------------------------------------------------------------------------------------------------------------------------------
def candidate_sequences_from_fasta(fasta_file):
    fasta = open("sample.fa")   #read sequences from the sample.fa file
    candidates = [] #result list of the protospacers

    for sequence in fasta:
        sequence = sequence.strip()

        if sequence[0] == ">":
            continue
        else:
            c = candidate_sequence(sequence, 20, "NGG")
            if c:
                candidates.extend(c)
    return candidates

print(len(candidate_sequences_from_fasta("sample.fa")))
#--------------------------------------------------------------------------------------------------------------------------------------------------
