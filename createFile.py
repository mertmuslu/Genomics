from random import randint
import random

def createFile(fileName):
    lines = []
    k = 20  # Define k for k-mer length
    with open(fileName, 'w') as f:
        for i in range(10):
            line = ''.join(random.choice('ATGC') for _ in range(500))
            lines.append(line)
        lines = implant(lines, k)  
        
        with open(fileName, 'w') as f:  
            f.writelines(line + '\n' for line in lines)
            

# Generate a random k-mer of length k to implant in the file    
def random_k_mer(k):
    k_mer = []
    for _ in range(k):
        x = random.randint(1, 4)
        if x == 1:
            k_mer.append('A')
        elif x == 2:
            k_mer.append('T')
        elif x == 3:
            k_mer.append('G')
        elif x == 4:
            k_mer.append('C')
    return ''.join(k_mer)

# Implant the random k-mer into the file at random positions in the lines
def mutate(motifs, k):
    positions = random.sample(range(k), 2)  # Select 4 unique random positions
    for j in positions:
        current_base = motifs[j]
        if current_base == 'A':
            new_base = random.choice(['G', 'C', 'T'])
        elif current_base == 'T':
            new_base = random.choice(['A', 'G', 'C'])
        elif current_base == 'G':
            new_base = random.choice(['A', 'T', 'C'])
        elif current_base == 'C':
            new_base = random.choice(['A', 'T', 'G'])
        motifs = motifs[:j] + new_base + motifs[j + 1:]
    return motifs

def implant(lines, k):
    k_mer = random_k_mer(k)
    print("Inserted k-mer is: " + k_mer)
    for i in range(len(lines)):
        x = random.randint(0, 500 - k)
        temp_k_mer = mutate(k_mer, k)
        print("Implanting " + temp_k_mer + " in line " + str(i))
        lines[i] = lines[i][:x] + temp_k_mer + lines[i][x + k:]
    return lines  # Return the updated lines

createFile("input.txt")