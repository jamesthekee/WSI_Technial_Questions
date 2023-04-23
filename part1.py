# import pandas as pd

# df = pd.read_csv("seq.tsv", sep="\t")


# Question 1 & 2 (without pandas)
def loadData():
    data = []
    with open("seq.tsv", mode="r") as f:
        lines = f.read().splitlines()
        cols = lines[0].split("\t")
        for k in lines[1:]:
            x = k.split("\t")
            x[0] = int(x[0])
            data.append(x)
            
    print(cols)
            
    # TEST DATA
    # Check there is the expected number of cols,
    # And of the expected datatype
    
    for row in data:
        assert len(row) == 3
        assert type(row[0]) == int   
        assert type(row[1]) == str 
        assert type(row[2]) == str        
    
    return data
    
def count_differences(seqa, seqb):
    # Count number of deletions, insertions, mutations edits.
    # Since an edit could represent an infinite number of edits, e.g TAT->TACT could be 3 deletions and 4 additions, 
    # I will assume we are only interested in the most short sequence of edits
    
    # Implemented with an adapted needleman-wunsch,
    # so that is resembles calculating edit distance more.
    
    MUTATION = -1
    DELETION = -1
    INSERTION = -1
    
    dp = [[-(x+y) for x in range(len(seqa)+1)] for y in range(len(seqb)+1)]
    
    # Fill in dp matrix
    for j in range(1, len(seqb)+1):
        for i in range(1, len(seqa)+1):
            case1 = dp[j-1][i-1]
            if seqa[i-1] != seqb[j-1]:
                case1 += MUTATION
            case2 = dp[j][i-1]+DELETION
            case3 = dp[j-1][i]+INSERTION
            dp[j][i] = max(case1, case2, case3)
    
    # Backtrack from dp[len(b)][ len(a)]
    x = len(seqb)-1
    y = len(seqa)-1
    mutations = 0
    deletions = 0
    insertions = 0
    while x>0 or y>0:
        if dp[x][y] == dp[x-1][y-1]:
            x -= 1
            y -= 1
        elif dp[x][y] == dp[x-1][y-1] + MUTATION:
            x -= 1
            y -= 1
            mutations += 1
        elif dp[x][y] == dp[x][y-1] + DELETION:
            y -= 1
            deletions += 1
        elif dp[x][y] == dp[x-1][y] + INSERTION:  
            x -= 1
            insertions += 1
    
    return mutations, deletions, insertions
    
        

def question3(data):
    for _, seq, edit in data:
        result = count_differences(seq,edit)
        print(seq, edit)
        print(f"Mutations: {result[0]}\nDeletions: {result[1]}\nInsertions: {result[2]}")
        

question3(loadData())
