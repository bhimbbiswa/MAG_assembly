import os

configfile: "config.yaml"

# Accessions: single id or list file
if "accession" in config:
    accession_ids = [str(config["accession"]).strip()]
else:
    with open(config["accession_list"]) as f:
        accession_ids = [line.strip() for line in f if line.strip()]

rule all:
    input:
        # cleaned paired reads (post host-depletion)
        expand(os.path.join(config["output_dir"], "{accession}", "{accession}_clean_1.fastq.gz"), accession=accession_ids),
        expand(os.path.join(config["output_dir"], "{accession}", "{accession}_clean_2.fastq.gz"), accession=accession_ids),
        # cleaned singletons (post host-depletion)
        expand(os.path.join(config["output_dir"], "{accession}", "{accession}_clean_singletons.fastq.gz"), accession=accession_ids),
        # assembly outputs
        expand(os.path.join(config["output_dir"], "{accession}", "spades", "contigs.fasta"), accession=accession_ids),
        expand(os.path.join(config["output_dir"], "{accession}", "spades_done.flag"), accession=accession_ids)

# Download FASTQ (paired SRA)
rule download_fastq:
    output:
        forward = temp(os.path.join(config["output_dir"], "{accession}", "{accession}_1.fastq.gz")),
        reverse = temp(os.path.join(config["output_dir"], "{accession}", "{accession}_2.fastq.gz")),
        flag    = temp(os.path.join(config["output_dir"], "{accession}", "download_done.flag"))
    log:
        os.path.join(config["output_dir"], "{accession}", "download.log")
    params:
        out_dir = lambda wc: os.path.join(config["output_dir"], wc.accession),
        extra   = config.get("fastq_dl_extra", "")
    conda: "envs/fastq-dl.yaml"
    shell:
        r"""
        mkdir -p {params.out_dir}
        fastq-dl --accession {wildcards.accession} --outdir {params.out_dir} {params.extra} &> {log}
        # Sanity check
        test -s {output.forward} && test -s {output.reverse}
        touch {output.flag}
        """

# Trim & QC
rule trim_and_quality_filter:
    input:
        forward = os.path.join(config["output_dir"], "{accession}", "{accession}_1.fastq.gz"),
        reverse = os.path.join(config["output_dir"], "{accession}", "{accession}_2.fastq.gz"),
        flag    = os.path.join(config["output_dir"], "{accession}", "download_done.flag")
    output:
        trimmed_forward = temp(os.path.join(config["output_dir"], "{accession}", "{accession}_trimmed_1.fq.gz")),
        trimmed_reverse = temp(os.path.join(config["output_dir"], "{accession}", "{accession}_trimmed_2.fq.gz")),
        flag            = temp(os.path.join(config["output_dir"], "{accession}", "trim_done.flag"))
    log:
        os.path.join(config["output_dir"], "{accession}", "trim.log")
    params:
        threads = int(config.get("trim_threads", 8)),
        out_dir = lambda wc: os.path.join(config["output_dir"], wc.accession),
        extra   = config.get("trim_galore_extra", "--quality 20 --length 50 --fastqc")
    conda: "envs/metagen_qc.yaml"
    shell:
        r"""
        mkdir -p {params.out_dir}
        trim_galore --paired -o {params.out_dir} -j {params.threads} {params.extra} \
            {input.forward} {input.reverse} &> {log}
        mv {params.out_dir}/{wildcards.accession}_1_val_1.fq.gz {output.trimmed_forward}
        mv {params.out_dir}/{wildcards.accession}_2_val_2.fq.gz {output.trimmed_reverse}
        touch {output.flag}
        """

# Map to host genome and keep ALL unmapped (pairs + singletons)
rule map_to_host_genome:
    input:
        trimmed_forward = os.path.join(config["output_dir"], "{accession}", "{accession}_trimmed_1.fq.gz"),
        trimmed_reverse = os.path.join(config["output_dir"], "{accession}", "{accession}_trimmed_2.fq.gz"),
        flag            = os.path.join(config["output_dir"], "{accession}", "trim_done.flag")
    output:
        unmapped_bam = temp(os.path.join(config["output_dir"], "{accession}", "{accession}_unmapped_all.bam")),
        flag         = temp(os.path.join(config["output_dir"], "{accession}", "mapping_done.flag"))
    log:
        os.path.join(config["output_dir"], "{accession}", "mapping.log")
    params:
        ref     = config["reference_file"],   # BWA index prefix (.bwt etc. must exist)
        threads = int(config.get("map_threads", 8)),
        extra   = config.get("bwa_extra", "-M")
    conda: "envs/metagen_qc.yaml"
    shell:
        r"""
        # Ensure BWA index exists
        test -s {params.ref}.bwt || (echo "Missing BWA index for {params.ref}" >&2; exit 2)

        bwa mem {params.extra} -t {params.threads} {params.ref} {input.trimmed_forward} {input.trimmed_reverse} 2> {log} | \
        # Keep primary unmapped (-f 4), drop secondary (256) and supplementary (2048) => -F 2304
        samtools view -@ {params.threads} -b -f 4 -F 2304 -o {output.unmapped_bam} -
        touch {output.flag}
        """

# Sort & convert BAM → FASTQ (emit pairs AND singletons)
rule sort_and_convert:
    input:
        unmapped_bam = os.path.join(config["output_dir"], "{accession}", "{accession}_unmapped_all.bam"),
        flag         = os.path.join(config["output_dir"], "{accession}", "mapping_done.flag")
    output:
        sorted_bam  = temp(os.path.join(config["output_dir"], "{accession}", "{accession}_unmapped_all_sorted.bam")),
        clean_fwd   = temp(os.path.join(config["output_dir"], "{accession}", "{accession}_clean_1.fastq")),
        clean_rev   = temp(os.path.join(config["output_dir"], "{accession}", "{accession}_clean_2.fastq")),
        clean_single= temp(os.path.join(config["output_dir"], "{accession}", "{accession}_clean_singletons.fastq")),
        flag        = temp(os.path.join(config["output_dir"], "{accession}", "sort_done.flag"))
    log:
        os.path.join(config["output_dir"], "{accession}", "sort_convert.log")
    params:
        threads = int(config.get("samtools_threads", 8))
    conda: "envs/metagen_qc.yaml"
    shell:
        r"""
        samtools sort -@ {params.threads} -n {input.unmapped_bam} -o {output.sorted_bam} &>> {log}
        samtools fastq -@ {params.threads} -n \
            -1 {output.clean_fwd} \
            -2 {output.clean_rev} \
            -s {output.clean_single} \
            -0 /dev/null {output.sorted_bam} &>> {log}
        touch {output.flag}
        """

# Compress (pigz) & tidy
rule compress_and_cleanup:
    input:
        clean_fwd    = os.path.join(config["output_dir"], "{accession}", "{accession}_clean_1.fastq"),
        clean_rev    = os.path.join(config["output_dir"], "{accession}", "{accession}_clean_2.fastq"),
        clean_single = os.path.join(config["output_dir"], "{accession}", "{accession}_clean_singletons.fastq"),
        flag         = os.path.join(config["output_dir"], "{accession}", "sort_done.flag")
    output:
        cleaned_forward = protected(os.path.join(config["output_dir"], "{accession}", "{accession}_clean_1.fastq.gz")),
        cleaned_reverse = protected(os.path.join(config["output_dir"], "{accession}", "{accession}_clean_2.fastq.gz")),
        cleaned_single  = protected(os.path.join(config["output_dir"], "{accession}", "{accession}_clean_singletons.fastq.gz")),
        flag            = temp(os.path.join(config["output_dir"], "{accession}", "compress_done.flag"))
    log:
        os.path.join(config["output_dir"], "{accession}", "compress.log")
    params:
        out_dir = lambda wc: os.path.join(config["output_dir"], wc.accession)
    conda: "envs/metagen_qc.yaml"
    shell:
        r"""
        pigz -p {config.get("pigz_threads", 4)} {input.clean_fwd} {input.clean_rev} {input.clean_single} &> {log}
        rm -f {params.out_dir}/{wildcards.accession}_unmapped_all.bam \
              {params.out_dir}/{wildcards.accession}_unmapped_all_sorted.bam
        touch {output.flag}
        """

# metaSPAdes assembly (uses singletons via -s)
rule assemble_metaspades:
    input:
        cleaned_forward = os.path.join(config["output_dir"], "{accession}", "{accession}_clean_1.fastq.gz"),
        cleaned_reverse = os.path.join(config["output_dir"], "{accession}", "{accession}_clean_2.fastq.gz"),
        cleaned_single  = os.path.join(config["output_dir"], "{accession}", "{accession}_clean_singletons.fastq.gz"),
        flag            = os.path.join(config["output_dir"], "{accession}", "compress_done.flag")
    output:
        contigs = os.path.join(config["output_dir"], "{accession}", "spades", "contigs.fasta"),
        flag    = os.path.join(config["output_dir"], "{accession}", "spades_done.flag")
    log:
        os.path.join(config["output_dir"], "{accession}", "spades", "spades.log")
    params:
        out_dir      = lambda wc: os.path.join(config["output_dir"], wc.accession, "spades"),
        spades_bin   = config.get("spades_bin", "metaspades.py"),
        threads      = int(config.get("spades_threads", 16)),
        memory       = int(config.get("spades_memory_gb", 128)),
        kmers        = config.get("spades_kmers", "21,33,55"),
        phred_offset = str(config.get("spades_phred_offset", 33)),
        omp_threads  = int(config.get("spades_omp_threads", 16)),
        mkl_threads  = int(config.get("spades_mkl_threads", 1)),
        extra        = config.get("spades_extra", "")
    threads: int(config.get("spades_threads", 16))
    conda: "envs/spades.yaml"
    shell:
        r"""
        mkdir -p {params.out_dir}
        export OMP_NUM_THREADS={params.omp_threads}
        export MKL_NUM_THREADS={params.mkl_threads}

        {params.spades_bin} \
          --phred-offset {params.phred_offset} \
          -o {params.out_dir} \
          -t {params.threads} \
          -m {params.memory} \
          -k {params.kmers} \
          -1 {input.cleaned_forward} \
          -2 {input.cleaned_reverse} \
          -s {input.cleaned_single} \
          {params.extra} &> {log}

        test -s {output.contigs}
        touch {output.flag}
        """
