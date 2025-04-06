import random
import time
from collections import Counter

def read_dna_file(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f.readlines()]

def random_kmer(dna, k):
    start = random.randint(0, len(dna) - k)
    return dna[start:start + k]

def profile_with_pseudocounts(motifs, k):
    profile = {'A':[1]*k, 'C':[1]*k, 'G':[1]*k, 'T':[1]*k}
    for motif in motifs:
        for i, char in enumerate(motif):
            profile[char][i] += 1
    total = [sum(x) for x in zip(profile['A'], profile['C'], profile['G'], profile['T'])]
    for base in 'ACGT':
        profile[base] = [val / total[i] for i, val in enumerate(profile[base])]
    return profile

def most_probable_kmer(dna, k, profile):
    max_prob = -1
    most_prob = dna[0:k]
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        prob = 1
        for j, char in enumerate(kmer):
            prob *= profile[char][j]
        if prob > max_prob:
            max_prob = prob
            most_prob = kmer
    return most_prob

def score(motifs):
    consensus = ''
    total_score = 0
    for i in range(len(motifs[0])):
        counts = Counter(motif[i] for motif in motifs)
        max_freq = counts.most_common(1)[0][1]
        total_score += len(motifs) - max_freq
    return total_score

def consensus(motifs):
    consensus_str = ''
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        most_common = Counter(column).most_common(1)[0][0]
        consensus_str += most_common
    return consensus_str
def randomized_motif_search(dna_list, k, iterations=1000):
    best_motifs = [random_kmer(dna, k) for dna in dna_list]
    best_score = score(best_motifs)

    for _ in range(iterations):
        profile = profile_with_pseudocounts(best_motifs, k)
        motifs = [most_probable_kmer(dna, k, profile) for dna in dna_list]
        current_score = score(motifs)
        if current_score < best_score:
            best_motifs = motifs
            best_score = current_score
        else:
            break
    return best_motifs, best_score
def weighted_choice(kmer_probs):
    total = sum(prob for kmer, prob in kmer_probs)
    r = random.uniform(0, total)
    upto = 0
    for kmer, prob in kmer_probs:
        if upto + prob >= r:
            return kmer
        upto += prob
    return kmer_probs[-1][0]

def gibbs_sampler(dna_list, k, max_iter=1000):
    motifs = [random_kmer(dna, k) for dna in dna_list]
    best_motifs = motifs[:]
    best_score = score(motifs)
    stable_count = 0

    for it in range(max_iter):
        i = random.randint(0, len(dna_list)-1)
        excluded = motifs[:i] + motifs[i+1:]
        profile = profile_with_pseudocounts(excluded, k)

        kmers = []
        for j in range(len(dna_list[i]) - k + 1):
            kmer = dna_list[i][j:j+k]
            prob = 1
            for m, char in enumerate(kmer):
                prob *= profile[char][m]
            kmers.append((kmer, prob))
        motifs[i] = weighted_choice(kmers)

        current_score = score(motifs)
        if current_score < best_score:
            best_motifs = motifs[:]
            best_score = current_score
            stable_count = 0
        else:
            stable_count += 1
            if stable_count >= 10:
                break
    return best_motifs, best_score
def run_algorithms(file_path, k):
    dna_list = read_dna_file(file_path)

    print(f"\nRunning for k = {k}...\n")

    # Randomized Motif Search
    start_rms = time.time()
    rms_motifs, rms_score = randomized_motif_search(dna_list, k)
    end_rms = time.time()
    print("RMS Best Score:", rms_score)
    print("RMS Consensus:", consensus(rms_motifs))
    print("RMS Time:", round(end_rms - start_rms, 3), "seconds")

    # Gibbs Sampler
    start_gibbs = time.time()
    gibbs_motifs, gibbs_score = gibbs_sampler(dna_list, k)
    end_gibbs = time.time()
    print("\nGibbs Best Score:", gibbs_score)
    print("Gibbs Consensus:", consensus(gibbs_motifs))
    print("Gibbs Time:", round(end_gibbs - start_gibbs, 3), "seconds")
if __name__ == "__main__":
    for k in [20]:
        run_algorithms("input.txt", k)
