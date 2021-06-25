# ancient_LoF
Data source files: 
  https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz 
  https://gnomad-public-us-east-1.s3.amazonaws.com/papers/2019-flagship-lof/v1.0/gnomad.v2.1.1.all_lofs.txt.bgz
  sgdp_data (5 samples)
  Archaics_data (3 samples)
  
Commands to produce files for lof variant counts and corresponding synonymous variant counts:
(the synonymous variants are from the same gene as the LoF variant and are chosen based on similar allele frequency to LoF variant)
  
python lof_count.py gnomad.v2.1.1.all_lofs.txt.bgz /project/mathilab/data/ancient/snpAD/VCF/Archaics/Archaics_chr*_mq25_mapab100_nonref_only.vcf.gz > outputs/arch_lofCount.txt

python lof_list.py gnomad.v2.1.1.all_lofs.txt.bgz /project/mathilab/data/ancient/snpAD/VCF/Archaics/Archaics_chr*_mq25_mapab100_nonref_only.vcf.gz > outputs/arch_lofList.txt

python exomes_command.py outputs/arch_lofList.txt /project/mathilab/data/gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz | gzip > outputs/arch_syn_freqs.gz

zgrep LoF outputs/arch_syn_freqs.gz | gzip > outputs/arch_lof_freqs.gz

python syn_list.py outputs/arch_lof_freqs.gz outputs/arch_syn_freqs.gz --out outputs/arch_synList.txt

python syn_count.py outputs/arch_synList.txt /project/mathilab/data/ancient/snpAD/VCF/Archaics/Archaics_chr*_mq25_mapab100_nonref_only.vcf.gz > outputs/arch_synCount.txt
