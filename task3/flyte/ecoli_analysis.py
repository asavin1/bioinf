from flytekit import task, workflow, Resources
from typing import Tuple
import subprocess
import os
import re

@task(limits=Resources(cpu="2", mem="2Gi"))
def run_fastqc(fastq_path: str) -> str:
    """Run FastQC quality control"""
    output_dir = "fastqc_results"
    os.makedirs(output_dir, exist_ok=True)
    subprocess.run(["fastqc", fastq_path, "-o", output_dir], check=True)
    return f"{output_dir}/{os.path.basename(fastq_path)}_fastqc.html"

@task(limits=Resources(cpu="4", mem="8Gi"))
def align_reads(fastq_path: str, ref_genome: str) -> str:
    """Align reads to reference"""
    sam_file = "aligned.sam"
    subprocess.run([
        "bwa", "mem",
        "-t", "4",
        ref_genome,
        fastq_path,
        "-o", sam_file
    ], check=True)
    return sam_file

@task
def calculate_stats(sam_file: str) -> float:
    """Calculate mapping statistics"""
    # Конвертируем SAM в BAM
    bam_file = "aligned.bam"
    subprocess.run(["samtools", "view", "-S", "-b", sam_file, "-o", bam_file], check=True)
    
    # Сортируем BAM файл
    sorted_bam = "aligned_sorted.bam"
    subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam], check=True)
    
    # Индексируем
    subprocess.run(["samtools", "index", sorted_bam], check=True)
    
    # Получаем статистику (правильный синтаксис)
    flagstat_output = subprocess.check_output(
        ["samtools", "flagstat", sorted_bam],
        universal_newlines=True
    )
    
    # Парсим процент картирования
    match = re.search(r"(\d+) \+ \d+ mapped \((\d+\.\d+)%", flagstat_output)
    if not match:
        raise ValueError("Could not parse mapping percentage from flagstat")
    
    return float(match.group(2))

@task
def evaluate_qc(mapped_percent: float) -> str:
    """Evaluate QC result"""
    return "PASS" if mapped_percent >= 75.0 else "FAIL"

@workflow
def ecoli_qc_pipeline(
    fastq_path: str = "data/SRR33692911.fastq",
    ref_genome: str = "data/ecoli_ref.fna"
) -> Tuple[float, str]:
    """Main analysis workflow"""
    # Quality control
    fastqc_report = run_fastqc(fastq_path=fastq_path)
    
    # Alignment
    sam_file = align_reads(fastq_path=fastq_path, ref_genome=ref_genome)
    
    # Statistics
    mapped_percent = calculate_stats(sam_file=sam_file)
    
    # Evaluation
    status = evaluate_qc(mapped_percent=mapped_percent)
    
    return (mapped_percent, status)

if __name__ == "__main__":
    # Test run
    print("Running pipeline...")
    percent, status = ecoli_qc_pipeline()
    print(f"Result: {percent:.2f}% mapped - {status}")
