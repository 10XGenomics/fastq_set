'''
Compares the result of adapter trimming between the rust trimmer
and cutadapt for the configuration described. See `testlist`
in the main function
'''

import tenkit.fasta as tk_fasta
import cellranger.utils as cr_utils
import random
import json
import os
import subprocess
import time

SEQ_LEN = 150
# Random seq of SEQ_LEN bases
# Insert adapter at some position [0,SEQ_LEN) into 50% of reads
# Trim the read to SEQ_LEN bases
def insert_adapter_at_random(out_writer, adapter, num_seqs, err_rate, allow_indel):
    for i in xrange(num_seqs):
        hdr = str(i)
        seq = random_seq(SEQ_LEN)
        if random.random() > 0.5:
            pos = random.randint(0, SEQ_LEN-1)
            seq = seq[:pos] + mutate(adapter, err_rate, allow_indel) + seq[pos:]
            seq = seq[:SEQ_LEN]
        tk_fasta.write_read_fasta(out_writer, hdr, seq)

# Mutate the sequence with an error rate
# Optionally allow indels
def mutate(seq, err_rate, allow_indel):
    mis_seq = random_seq(len(seq))
    ins_seq = random_seq(len(seq))
    mutated = ""
    for (s, m, i) in zip(seq, mis_seq, ins_seq):
        val = random.random()
        if val <= err_rate:
            if allow_indel:
                if val <= err_rate/2.0: # mismatch
                    mutated += m
                elif val <= 3.0 * err_rate / 4.0: # insert
                    mutated += s
                    mutated += i
                else: # deletion
                    pass
            else:
                mutated += m
        else:
            mutated += s
    return mutated

def random_seq(n):
    alphabets = "ACGT"
    return ''.join([alphabets[random.randint(0,3)] for _ in xrange(n)])

def run_adapter_trimmer(test, input, output, config_json):
    config = {
        'adapter': {
            'name': test.get('name', 'my_primer'),
            'end': test['end'],
            'location': test['location'],
            'seq': test['adapter'],
        },
        'input_fasta': input,
        'output_fasta': output,
    }
    json.dump(config, open(config_json, 'w'), indent=4)

    cmd = ['cargo', 'r', '--release', '--example', 'trim', '--', config_json]
    output = subprocess.check_output(cmd, stderr=subprocess.PIPE)

    elapsed_time = None
    for line in output.strip().split('\n'):
        if line.startswith('Elapsed time'):
            elapsed_time = float(line.strip().split(' ')[2][2:-1])
            break

    assert elapsed_time is not None
    return elapsed_time

def run_cutadapt(test, input, output):
    version = subprocess.check_output(['cutadapt', '--version']).strip()
    if test['location'] == 'non_internal':
        assert version in ['1.17', '1.18'] # Only supported in these versions
    cmd = ['cutadapt']
    if test['end'] == 'three_prime':
        cmd.append('-a')
    else: # Not doing any error checking
        cmd.append('-g')
    if test['location']=='anywhere':
        cmd.append(test['adapter'])
    elif test['location']=='anchored':
        if test['end'] == 'three_prime':
            cmd.append(test['adapter'] + '$')
        else: # Five prime, Not doing any error checking
            cmd.append('^' + test['adapter'])
    else: # Non internal
        if test['end'] == 'three_prime':
            cmd.append(test['adapter'] + 'X')
        else: # Five prime, Not doing any error checking
            cmd.append('X' + test['adapter'])
    cmd.extend(['-o', output])
    cmd.extend(['--overlap', '5'])
    cmd.append(input)

    t0 = time.time()
    _ = subprocess.check_output(cmd)
    return float(time.time() - t0)

def check(input, trimmer_out, cutadapt_out):

    correct_untrimmed = 0 # Not trimmed by both
    correct_trimmed = 0 # Consistently trimmed by both
    rust_trimmed_only = 0 # Trimmed only by rust
    rust_untrimmed_only = 0 # Not trimmed only by rust
    inconsistent_trimmed = 0 # Trimmed inconsistently

    with open(cutadapt_out) as fcut, open(trimmer_out) as frust, open(input) as finput:
        cutiter = cr_utils.get_fasta_iter(fcut)
        rustiter = cr_utils.get_fasta_iter(frust)
        inputiter = cr_utils.get_fasta_iter(finput)
        for (c, r, u) in zip(cutiter, rustiter, inputiter):
            if c==u and r==u:
                correct_untrimmed += 1
            elif c==u:
                rust_trimmed_only += 1
            elif r==u:
                rust_untrimmed_only += 1
            elif c==r or abs(len(r[1])-len(c[1]))<=2: # Allow upto a slop of two
                correct_trimmed += 1
            else:
                inconsistent_trimmed += 1
        
    # Sensitivity = (Rust & Cutadapt) / Cutadapt
    sensitivity = float(correct_trimmed + inconsistent_trimmed) / float(correct_trimmed + inconsistent_trimmed + rust_untrimmed_only)
    # PPV = (Rust & Cutadapt)/Rust
    ppv = float(correct_trimmed + inconsistent_trimmed) / float(correct_trimmed + inconsistent_trimmed + rust_trimmed_only)
    # Concordance = (Rust & Cutadapt & Rust==Cutadpt) / (Rust & Cutadapt)
    concordance = float(correct_trimmed) / float(correct_trimmed + inconsistent_trimmed)

    return sensitivity, ppv, concordance

def print_result(test, sensitivity, ppv, concordance, speedup):
    if sensitivity > 0.99 and ppv > 0.99 and concordance > 0.99:
        result = 'SUCCESS'
    else:
        result = 'FAILED!'
    print("> {} -> PPV: {:.2%} Sensitivity: {:.2%} Concordance: {:.2%} \t SPEEDUP: {:.2} \t Adapter: {} end: {} location: {}".format(\
        result, ppv, sensitivity, concordance, speedup, test['adapter'], test['end'], test['location']))

if __name__=="__main__":

    random.seed(0)

    testlist = [
        {'adapter': 'AGATGTGTATAAGAGACAG', 'end': 'three_prime', 'location': 'anywhere', 'name': 'Agore-ME'},
        {'adapter': 'CTGTCTCTTATACACATCT', 'end': 'three_prime', 'location': 'anywhere', 'name': 'Agore-MErc'},
        {'adapter': 'TTTCTTATATGGG', 'end': 'five_prime', 'location': 'anchored', 'name': 'VDJ-spacer'},
        {'adapter': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 'end': 'three_prime', 'location': 'anywhere', 'name': 'VDJ-R2rc'},
        {'adapter': 'ATCTCGTATGCCGTCTTCTGCTTG', 'end': 'three_prime', 'location': 'anywhere', 'name': 'VDJ-P7rc'},
        {'adapter': 'AAAAAAAAAAAAAAAAAAAA', 'end': 'three_prime', 'location': 'anywhere', 'name': 'VDJ-polyA'},
        {'adapter': 'AAGCAGTGGTATCAACGCAGAGTACAT', 'end': 'three_prime', 'location': 'anywhere', 'name': 'VDJ-rtprimer'},
        {'adapter': 'CCCATATAAGAAA', 'end': 'three_prime', 'location': 'anywhere', 'name': 'VDJ-spacer-rc'},
        {'adapter': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', 'end': 'three_prime', 'location': 'anywhere', 'name': 'VDJ-R1rc'},
        {'adapter': 'AGATCTCGGTGGTCGCCGTATCATT', 'end': 'three_prime', 'location': 'anywhere', 'name': 'VDJ-P5rc'},
    ]

    err_rates = [0, 0.01, 0.02, 0.05, 0.1]

    num_seqs = 10000 # Number of sequences for each error rate
    
    abs_dir_path = os.path.dirname(os.path.realpath(__file__))
    trimmer_out = os.path.join(abs_dir_path, "trimmer_out.fa")
    cutadapt_out = os.path.join(abs_dir_path, "cutadapt_out.fa")
    input_fa = os.path.join(abs_dir_path, "input.fa")
    config_json = os.path.join(abs_dir_path, "config.json")

    for i, test in enumerate(testlist):
        print 'Generating data for test ', test
        adapter = test['adapter']
        with open(input_fa, 'w') as writer:
            for err_rate in err_rates:
                insert_adapter_at_random(writer, adapter, num_seqs, err_rate, False)
                insert_adapter_at_random(writer, adapter, num_seqs, err_rate, True)
        print 'Executing rust adapter trimmer...'
        trimmer_time = run_adapter_trimmer(test, input_fa, trimmer_out, config_json)
        print 'Done in {} seconds..'.format(trimmer_time)
        print 'Executing cutadapt...'
        cutadapt_time = run_cutadapt(test, input_fa, cutadapt_out)
        print 'Done in {} seconds..'.format(cutadapt_time)
        speedup = cutadapt_time/trimmer_time
        sensitivity, ppv, concordance = check(input_fa, trimmer_out, cutadapt_out)
        print '> TEST RESULT {}/{}'.format(i+1, len(testlist))
        print_result(test, sensitivity, ppv, concordance, speedup)
        
