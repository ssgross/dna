#!/usr/bin/python

import string
import sys
import subprocess

bwa_bin = '/home/ssgross/bwa-0.7.4/bwa'
samtools_bin = '/home/ssgross/samtools-0.1.19/samtools'
java_bin = '/home/ssgross/jre1.7.0_25/bin/java'
gatk_jar = '/home/ssgross/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar'
picard_dir = '/home/ssgross/picard-tools-1.91'
snpeff_jar = '/home/ssgross/snpEff/snpEff.jar'
snpeff_config = '/home/ssgross/snpEff/snpEff.config'

reference = '/home/ssgross/dna/hs37d5.fa'
dbsnp = '/home/ssgross/dna/dbsnp_137.b37.vcf'
thousand_indels = '/home/ssgross/dna/1000G_phase1.indels.b37.vcf'
mills_indels = '/home/ssgross/dna/Mills_and_1000G_gold_standard.indels.b37.vcf'
hapmap = '/home/ssgross/dna/hapmap_3.3.b37.vcf'
omni = '/home/ssgross/dna/1000G_omni2.5.b37.vcf'

def RunCommand(command):
  print 'Running %s' % command
  subprocess.call(command, shell=True)


sample_name = sys.argv[1]
in_prefix = sys.argv[2]
num_lanes = int(sys.argv[3])
files_per_lane = int(sys.argv[4])


def Map():
  sorted_bams = []
  for lane_num in range(1, num_lanes + 1):
    for file_num in range(1, files_per_lane + 1):
      read_group_header = ('"@RG\\tID:%s_L_%d_%d\\tSM:%s\\tPL:ILLUMINA\\tLB:%s"'
                           % (sample_name, lane_num, file_num, sample_name, sample_name))
      mapped_prefix = 'L00%d_00%d' % (lane_num, file_num)
      bwa_command = ('%s mem -M -R %s -t 8 %s %s_L00%d_R1_00%d.fastq %s_L00%d_R2_00%d.fastq > %s.sam'
                     % (bwa_bin, read_group_header, reference, in_prefix, lane_num, file_num, in_prefix,
                        lane_num, file_num, mapped_prefix))
      RunCommand(bwa_command)
      compress_command = '%s view -bS %s.sam > %s.bam' % (samtools_bin, mapped_prefix, mapped_prefix)
      RunCommand(compress_command)
      sort_command = ('%s sort -m 2000000000 %s.bam %s_sorted'
                      % (samtools_bin, mapped_prefix, mapped_prefix))
      RunCommand(sort_command)
      sorted_bams.append('%s_sorted.bam' % mapped_prefix)
  
  merge_inputs = string.join(['INPUT=%s' % x for x in sorted_bams]) 
  merge_command = ('%s -Xmx2g -jar %s/MergeSamFiles.jar %s OUTPUT=mapped.bam'
                   % (java_bin, picard_dir, merge_inputs))
  RunCommand(merge_command)


def IndexBam(bam_file):
  command = '%s index %s' % (samtools_bin, bam_file)
  RunCommand(command)


def MarkDuplicates():
  mark_command = ('%s -Xmx2g -jar %s/MarkDuplicates.jar INPUT=mapped.bam OUTPUT=deduped.bam '
                  'METRICS_FILE=dedup_metrics.txt' % (java_bin, picard_dir))
  RunCommand(mark_command)
  IndexBam('deduped.bam')


def Realign():
  candidates_command = ('%s -Xmx2g -jar %s -T RealignerTargetCreator -nt 8 -I deduped.bam -R %s '
                        '-known %s -known %s -o realign_targets.intervals'
                        % (java_bin, gatk_jar, reference, thousand_indels, mills_indels))
  RunCommand(candidates_command)
  realign_command = ('%s -Xmx2g -jar %s -T IndelRealigner -R %s -I deduped.bam '
                     '-targetIntervals realign_targets.intervals -known %s -known %s '
                     '-o realigned.bam'
                     % (java_bin, gatk_jar, reference, thousand_indels, mills_indels))
  RunCommand(realign_command)


def RecalibrateBases():
  recal_data_command = ('%s -Xmx2g -jar %s -T BaseRecalibrator -nct 8 -I realigned.bam -R %s '
                        '-knownSites %s -o recal_data.txt'
                        % (java_bin, gatk_jar, reference, dbsnp))
  RunCommand(recal_data_command)
  print_reads_command = ('%s -Xmx2g -jar %s -T PrintReads -nct 8 -I realigned.bam -R %s '
                         '--BQSR recal_data.txt -o recalibrated.bam'
                         % (java_bin, gatk_jar, reference))
  RunCommand(print_reads_command)


def CallVariants():
  call_command = ('%s -Xmx2g -jar %s -T HaplotypeCaller -nct 8 -I recalibrated.bam -R %s '
                  '--dbsnp %s -stand_emit_conf 10.0 -o raw_variants.vcf'
                  % (java_bin, gatk_jar, reference, dbsnp))
  RunCommand(call_command)


def RecalibrateVariants():
  snp_recal_command = ('%s -Xmx2g -jar %s -T VariantRecalibrator -nt 8 -input raw_variants.vcf -R %s '
                       '--maxGaussians 4 '
                       '-resource:hapmap,known=false,training=true,truth=true,prior=15.0 %s '
                       '-resource:omni,known=false,training=true,truth=false,prior=12.0 %s '
                       '-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 %s '
                       '-an QD -an MQRankSum -an ReadPosRankSum -an FS '
                       '-mode SNP -recalFile snp.recal -tranchesFile snp.tranches'
                       % (java_bin, gatk_jar, reference, hapmap, omni, dbsnp))
  RunCommand(snp_recal_command)
  indel_recal_command = ('%s -Xmx2g -jar %s -T VariantRecalibrator -nt 8 -input raw_variants.vcf -R %s '
                         '--maxGaussians 2 '
                         '-resource:mills,known=false,training=true,truth=true,prior=12.0 %s '
                         '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %s '
                         '-an FS -an ReadPosRankSum -an MQRankSum '
                         '-mode INDEL -recalFile indel.recal -tranchesFile indel.tranches'
                         % (java_bin, gatk_jar, reference, mills_indels, dbsnp))
  RunCommand(indel_recal_command)
  apply_snp_recal_command = ('%s -Xmx2g -jar %s -T ApplyRecalibration -nt 8 -input raw_variants.vcf -R %s '
                             '-mode SNP --ts_filter_level 99.9 -recalFile snp.recal -tranchesFile snp.tranches '
                             '-o snp_recal_variants.vcf'
                             % (java_bin, gatk_jar, reference))
  RunCommand(apply_snp_recal_command)
  apply_indel_recal_command = ('%s -Xmx2g -jar %s -T ApplyRecalibration -nt 8 -input snp_recal_variants.vcf -R %s '
                             '-mode INDEL --ts_filter_level 99.9 -recalFile indel.recal -tranchesFile indel.tranches '
                             '-o recal_variants.vcf'
                             % (java_bin, gatk_jar, reference))
  RunCommand(apply_indel_recal_command)


def AnnotateVariants():
  snpeff_command = ('%s -Xmx2g -jar %s -c %s GRCh37.72 -canon recal_variants.vcf > annotated_variants.vcf'
                    % (java_bin, snpeff_jar, snpeff_config))
  RunCommand(snpeff_command)


Map()
MarkDuplicates()
Realign()
RecalibrateBases()
CallVariants()
RecalibrateVariants()
AnnotateVariants()
