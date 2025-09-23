import time
import os
import pandas as pd
from cyvcf2 import VCF

#Data path examples
#datadir = os.path.abspath("../../Data/")
#maindir   = datadir+"/GWAS/PRJNA1104541/"
#snpeff    = "/home/louie/Documents/Software/snpEff/snpEff.jar"

def IndexRefGenome(reference):
    start = time.time()
    # Index the reference genome with BWA
    os.system(f'bwa index {reference}')
    print(f'Reference index complete | {(time.time()-start)/60} minutes')

def GenerateVCF(maindir, snpeff):
    count = 0
    for i in os.listdir(maindir):
        start = time.time()
        if "fastq" in i and ".gz" not in i:
            count += 1
            stem = i.split(".fastq")[0]
            loc_short = f'{maindir}{stem}/{stem}'
            print(f'{"\n"*2}{"~"*10}\n({count}) - Main file - {i}')
            # Make directory
            if os.path.isdir(f'{maindir}{stem}') == False:
                os.system(f'mkdir {maindir}{stem}')
            os.system(f'cp {maindir}{i} {maindir}{stem}')
            # Align reads to reference (single-end)
            os.system(f'\nbwa mem -t 12 {reference} {loc_short}.fastq > {loc_short}_aligned.sam')
            # Convert SAM to BAM
            os.system(f'\nsamtools view -Sb {loc_short}_aligned.sam | samtools sort -o {loc_short}_aligned.bam')
            # Pileup generation
            os.system(f'\nbcftools mpileup -f {reference} {loc_short}_aligned.bam > {loc_short}.bcf')
            # Call variants
            os.system(f'\nbcftools call -mv -Ov -o {loc_short}.vcf {loc_short}.bcf')
            # Annotate VCF
            os.system(f'\njava -Xmx4g -jar {snpeff} GRCh38.p14 {loc_short}.vcf > {loc_short}_annotated.vcf')
            # Clean intermediary files
            os.system(f'rm {loc_short}.bcf {loc_short}_aligned.bam {loc_short}_aligned.sam {loc_short}.fastq {loc_short}.vcf')
            print(f'{i} analysis complete | {round((time.time()-start)/60, 3)} minutes')

def ParseVCF(maindir):
    #Set default df
    df        = pd.DataFrame(columns = ["Chr", "Pos","Ref", "Alt","Gene","Effect","Impact","Transcript","ProteinChange"])

    #Iterate through VCF files in directory, return df with "High" impact varients
    start = time.time()
    for path, subdirs, files in os.walk(maindir):
        for name in files:
            if "annotated.vcf" in name:
                vcf_file = os.path.join(path, name)
                print(vcf_file)
                records = []
                for variant in VCF(vcf_file):
                    ann_field = variant.INFO.get("ANN")
                    for ann in ann_field.split(","):
                        fields = ann.split("|")
                        records.append({
                            "Chr": variant.CHROM,
                            "Pos": variant.POS,
                            "Ref": variant.REF,
                            "Alt": ",".join(variant.ALT),
                            "Gene": fields[3],
                            "Effect": fields[1],
                            "Impact": fields[2],
                            "Transcript": fields[6],
                            "ProteinChange": fields[10]
                        })
                tmpdf = pd.DataFrame(records)
                df = pd.concat([df, tmpdf.loc[tmpdf["Impact"] == "HIGH"]])         
    print(f'DF generation complete | {round((time.time()-start)/60, 3)} minutes')
    
    return df