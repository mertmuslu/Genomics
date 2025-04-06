import random
import time

class rms_and_gibbs:
    def __init__(self, file_path, k):
        self.file_path = file_path
        self.k = k
        with open(file_path, 'r') as f:
            self.lines = f.readlines()
      
    # Select 10 motifs with random positions in the lines
    def select_motifs(self):
        motifs = []
        for i in range(10):
            x = random.randint(0, 500 - self.k)
            motifs.append(self.lines[i][x:x + self.k])
        return motifs

    # Create a profile matrix of the motifs
    def profile(self, motifs, laplace_smoothing=False):
        profile_motifs = [[0, 0, 0, 0] for _ in range(self.k)]  # Initialize 2D array
        numLines = len(motifs)  
        for i in range(self.k):
            count_A = 0
            count_T = 0
            count_G = 0
            count_C = 0
            for j in range(numLines):
                if motifs[j][i] == 'A':
                    count_A += 1
                elif motifs[j][i] == 'T':
                    count_T += 1
                elif motifs[j][i] == 'G':
                    count_G += 1
                elif motifs[j][i] == 'C':
                    count_C += 1

            if laplace_smoothing:
                # Apply Laplace smoothing
                profile_motifs[i][0] = (count_A + 1)
                profile_motifs[i][1] = (count_T + 1)
                profile_motifs[i][2] = (count_G + 1)
                profile_motifs[i][3] = (count_C + 1)
            else:
                profile_motifs[i][0] = count_A / numLines
                profile_motifs[i][1] = count_T / numLines
                profile_motifs[i][2] = count_G / numLines
                profile_motifs[i][3] = count_C / numLines

        return profile_motifs

    # Calculate the score of the motifs
    def calc_score(self, motifs):
        score = 0
        consensus = self.consensus_string(motifs)
        for i in range(10):
            for j in range(self.k):
                if motifs[i][j] != consensus[j]:
                    score += 1
        return score   
    
    # Consensus string of the motifs
    def consensus_string(self, motifs):
        consensus = ''
        for i in range(self.k):
            count_A = 0
            count_T = 0
            count_G = 0
            count_C = 0
            for j in range(len(motifs)):
                if motifs[j][i] == 'A':
                    count_A += 1
                elif motifs[j][i] == 'T':
                    count_T += 1
                elif motifs[j][i] == 'G':
                    count_G += 1
                elif motifs[j][i] == 'C':
                    count_C += 1

            if count_A >= max(count_T, count_G, count_C):
                consensus += 'A'
            elif count_T >= max(count_A, count_G, count_C):
                consensus += 'T'
            elif count_G >= max(count_A, count_T, count_C):
                consensus += 'G'
            else:
                consensus += 'C'
        
        return consensus    
    
    # RMS algorithm
    def randomized_motif_search(self): 
        best_overall_motifs = None
        best_overall_score = float('inf')
        start_time = time.time()
        for run in range(150):  # Repeat the process 100 times
            
            
            # Randomly select 10 motifs from the lines
            best_motifs = self.select_motifs()
            best_score = self.calc_score(best_motifs)
            
            iterations = 0 
            
            while True:
                iterations += 1
                # Profile(motifs)
                profile_motifs = self.profile(best_motifs)
                for i in range(0, 10):
                    max_prob = 0
                    for j in range(500 - self.k):
                        prob = 1
                        compare_motif = self.lines[i][j:j + self.k]
                        # calculation of probabilities for k-mer motif in the lines
                        for z in range(self.k):
                            if compare_motif[z] == 'A':
                                prob *= profile_motifs[z][0]
                            elif compare_motif[z] == 'T':
                                prob *= profile_motifs[z][1]
                            elif compare_motif[z] == 'G':
                                prob *= profile_motifs[z][2]
                            elif compare_motif[z] == 'C':
                                prob *= profile_motifs[z][3]
                        # keep highest probability motif and assign it to the best_motif
                        if prob > max_prob:                    
                            max_prob = prob
                            best_motif = compare_motif
                    # keep best motifs for each line
                    best_motifs[i] = best_motif
                
                current_score = self.calc_score(best_motifs)
                self.printBestMotifsInsideTxt_rms(best_motifs)

                # Check score, if remains the same for 10 iterations, print and break
                if current_score >= best_score:

                    break
                else:
                    best_score = current_score
            
            

            # Update the best overall motifs if the current run has a better score
            if best_score < best_overall_score:
                best_overall_score = best_score
                best_overall_motifs = best_motifs[:]
        
        # Print the best motifs after 100 runs
        print("\nBest motifs after 100 runs:")
        for motif in best_overall_motifs:
            print(motif)
        print("Best score: " + str(best_overall_score))
        print("Consensus string: " + self.consensus_string(best_overall_motifs))
        self.printBestMotifsInsideTxt_rms(best_overall_motifs)
        end_time = time.time()  # End timing
        print(f"RMS Execution time: {end_time - start_time:.2f} seconds")

    def gibbs_sampler(self):
        best_overall_motifs = None
        best_overall_score = float('inf')
        start_time_gibbs = time.time()  # Start timing
        counter = 0  # Counter to track iterations without improvement
        
        for run in range(10000):  # Repeat the process up to 10000 times
            
            best_motifs = self.select_motifs()
            best_score = self.calc_score(best_motifs)
            iterations = 0
            
            while True:
                max_prob = 0
                iterations += 1
                # pop random motif
                i = random.randint(0, 9)
                removed_motif = best_motifs.pop(i)
                # profile motif with Laplace smoothing enabled
                profile_motifs = self.profile(best_motifs, laplace_smoothing=True)
                
                for j in range(500 - self.k):
                    # compare_motif selected from removed motif's string
                    compare_motif = self.lines[i][j:j + self.k]
                    prob = 1
                    for z in range(self.k):
                        if compare_motif[z] == 'A':
                            prob *= profile_motifs[z][0]
                        elif compare_motif[z] == 'T':
                            prob *= profile_motifs[z][1]
                        elif compare_motif[z] == 'G':
                            prob *= profile_motifs[z][2]
                        elif compare_motif[z] == 'C':
                            prob *= profile_motifs[z][3]
                    if prob > max_prob:                    
                        max_prob = prob
                        best_motif = compare_motif

                # Add updated motif back to the list
                best_motifs.insert(i, best_motif)
                
                # Calculate score
                current_score = self.calc_score(best_motifs)
                
                # Check score, if remains the same for 10 iterations, print and break
                if current_score >= best_score:
                    break
                else:
                    best_score = current_score
            
            # Update the best overall motifs
            if best_score < best_overall_score:
                best_overall_score = best_score
                best_overall_motifs = best_motifs[:]
                counter = 0 
            else:
                counter += 1
            
            # Halt the loop if no improvement for 10 iterations
            if counter >= 1000:
                print("No improvement")
                break
        
        # Print the best motifs after the loop
        print("\nBest motifs after Gibbs sampling:")
        for motif in best_overall_motifs:
            print(motif)
        print("Best score: " + str(best_overall_score))
        print("Consensus string: " + self.consensus_string(best_overall_motifs))
        self.printBestMotifsInsideTxt_gibbs(best_overall_motifs)
        end_time_gibbs = time.time()  # End timing
        print(f"Gibbs Execution time: {end_time_gibbs - start_time_gibbs:.2f} seconds")
    
    def printBestMotifsInsideTxt_rms(self, motifs):
        with open('output_rms.txt', 'a') as f:  # Open in append mode
            for i in range(10):
                f.write(motifs[i] + '\n')
            f.write('\n' + self.consensus_string(motifs) + '\n')
            f.write(str(self.calc_score(motifs)) + '\n')
    def printBestMotifsInsideTxt_gibbs(self, motifs):
        with open('output_gibbs.txt', 'a') as f:  # Open in append mode
            for i in range(10):
                f.write(motifs[i] + '\n')
            f.write('\n' + self.consensus_string(motifs) + '\n')
            f.write(str(self.calc_score(motifs)) + '\n')
                    

k = int(input("Enter the k: "))
tempk = k
rms = rms_and_gibbs('input.txt', k)
rms.randomized_motif_search()

print("\n--------------------------------------------------\n")

rms.gibbs_sampler()
