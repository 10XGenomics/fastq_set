v=$(cutadapt --version | sed 's/\./_/g')
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o cutadapt_v${v}_anywhere_3p.fa --overlap 5 input.fa > /dev/null
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACX -o cutadapt_v${v}_non_internal_3p.fa --overlap 5 input.fa > /dev/null
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC$ -o cutadapt_v${v}_anchored_3p.fa --overlap 5 input.fa > /dev/null
cutadapt -g AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o cutadapt_v${v}_anywhere_5p.fa --overlap 5 input.fa > /dev/null
cutadapt -g XAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o cutadapt_v${v}_non_internal_5p.fa --overlap 5 input.fa > /dev/null
cutadapt -g ^AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o cutadapt_v${v}_anchored_5p.fa --overlap 5 input.fa > /dev/null

cutadapt -a AAAAAAAAAAAAAAAAAAAA -o cutadapt_v${v}_anywhere_3p_polyA.fa --overlap 5 input_polyA.fa > /dev/null
cutadapt -a AAAAAAAAAAAAAAAAAAAAX -o cutadapt_v${v}_non_internal_3p_polyA.fa --overlap 5 input_polyA.fa > /dev/null
cutadapt -a AAAAAAAAAAAAAAAAAAAA$ -o cutadapt_v${v}_anchored_3p_polyA.fa --overlap 5 input_polyA.fa > /dev/null
cutadapt -g AAAAAAAAAAAAAAAAAAAA -o cutadapt_v${v}_anywhere_5p_polyA.fa --overlap 5 input_polyA.fa > /dev/null
cutadapt -g XAAAAAAAAAAAAAAAAAAAA -o cutadapt_v${v}_non_internal_5p_polyA.fa --overlap 5 input_polyA.fa > /dev/null
cutadapt -g ^AAAAAAAAAAAAAAAAAAAA -o cutadapt_v${v}_anchored_5p_polyA.fa --overlap 5 input_polyA.fa > /dev/null