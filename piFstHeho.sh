#!/bin/bash
file=/mnt/ext35/ws/vcfsh/gea/gea_109_qc.recode.vcf
Pk=/mnt/ext35/ws/vcfsh/gea/21Pk_modified.txt
Pm=/mnt/ext35/ws/vcfsh/gea/27Pm_modified.txt
Pd=/mnt/ext35/ws/vcfsh/gea/19Pd_modified.txt
Py=/mnt/ext35/ws/vcfsh/gea/36Py_modified.txt
Pt=/mnt/ext35/ws/vcfsh/gea/3Pt_modified.txt
Pyv=/mnt/ext35/ws/vcfsh/gea/3Pyv_modified.txt


#1.pi
vcftools --vcf $file --keep $Pk --site-pi --out 21Pk_pi
vcftools --vcf $file --keep $Pm --site-pi --out 27Pm_pi
vcftools --vcf $file --keep $Pd --site-pi --out 19Pd_pi
vcftools --vcf $file --keep $Py --site-pi --out 36Py_pi
vcftools --vcf $file --keep $Pt --site-pi --out 3Pt_pi
vcftools --vcf $file --keep $Pyv --site-pi --out 3Pyv_p
#3.Ho
vcftools --vcf $file --keep $Pk --het --out 21Pk_het
vcftools --vcf $file --keep $Pm --het --out 27Pm_het
vcftools --vcf $file --keep $Pd --het --out 19Pd_het
vcftools --vcf $file --keep $Py --het --out 36Py_het
vcftools --vcf $file --keep $Pt --het --out 3Pt_het
vcftools --vcf $file --keep $Pyv --het --out 3Pyv_het

##需要基于结果手动计算Ho
awk '{if (NR>1) print $1, 1 - $2/$4}' 21Pk_het.het > 21Pk_Ho.txt
awk '{if (NR>1) print $1, 1 - $2/$4}' 27Pm_het.het > 27Pm_Ho.txt
awk '{if (NR>1) print $1, 1 - $2/$4}' 19Pd_het.het > 19Pd_Ho.txt
awk '{if (NR>1) print $1, 1 - $2/$4}' 36Py_het.het > 36Py_Ho.txt
awk '{if (NR>1) print $1, 1 - $2/$4}' 3Pt_het.het > 3Pt_Ho.txt
awk '{if (NR>1) print $1, 1 - $2/$4}' 3Pyv_het.het > 3Pyv_Ho.txt


#4.He
vcftools --vcf $file --keep $Pk --hardy --out 21Pk_hardy
vcftools --vcf $file --keep $Pm --hardy --out 27Pm_hardy
vcftools --vcf $file --keep $Pd --hardy --out 19Pd_hardy
vcftools --vcf $file --keep $Py --hardy --out 36Py_hardy
vcftools --vcf $file --keep $Pt --hardy --out 3Pt_hardy
vcftools --vcf $file --keep $Pyv --hardy --out 3Pyv_hardy

##基于结果手动计算期望杂合度（He=E_HET/总样本数）
awk 'NR>1 {split($4,a,"/"); He=a[2]/(a[1]+a[2]+a[3]); print $1,$2,He}' 21Pk_hardy.hwe > 21Pk_He.txt
awk 'NR>1 {split($4,a,"/"); He=a[2]/(a[1]+a[2]+a[3]); print $1,$2,He}' 27Pm_hardy.hwe > 27Pm_He.txt
awk 'NR>1 {split($4,a,"/"); He=a[2]/(a[1]+a[2]+a[3]); print $1,$2,He}' 19Pd_hardy.hwe > 19Pd_He.txt
awk 'NR>1 {split($4,a,"/"); He=a[2]/(a[1]+a[2]+a[3]); print $1,$2,He}' 36Py_hardy.hwe > 36Py_He.txt
awk 'NR>1 {split($4,a,"/"); He=a[2]/(a[1]+a[2]+a[3]); print $1,$2,He}' 3Pt_hardy.hwe > 3Pt_He.txt
awk 'NR>1 {split($4,a,"/"); He=a[2]/(a[1]+a[2]+a[3]); print $1,$2,He}' 3Pyv_hardy.hwe > 3Pyv_He.txt

##顺便也算一个Ho，并且跳过分母为0的行
awk 'NR>1 {split($3, a, "/"); total = a[1] + a[2] + a[3]; if (total > 0) {Ho = a[2] / total; print $1, $2, Ho} else {print "Skipping line " NR ": zero denominator" > "/dev/stderr"}}' 21Pk_hardy.hwe > 21Pk_hardyHo.txt
awk 'NR>1 {split($3, a, "/"); total = a[1] + a[2] + a[3]; if (total > 0) {Ho = a[2] / total; print $1, $2, Ho} else {print "Skipping line " NR ": zero denominator" > "/dev/stderr"}}' 27Pm_hardy.hwe > 27Pm_hardyHo.txt
awk 'NR>1 {split($3, a, "/"); total = a[1] + a[2] + a[3]; if (total > 0) {Ho = a[2] / total; print $1, $2, Ho} else {print "Skipping line " NR ": zero denominator" > "/dev/stderr"}}' 19Pd_hardy.hwe > 19Pd_hardyHo.txt
awk 'NR>1 {split($3, a, "/"); total = a[1] + a[2] + a[3]; if (total > 0) {Ho = a[2] / total; print $1, $2, Ho} else {print "Skipping line " NR ": zero denominator" > "/dev/stderr"}}' 36Py_hardy.hwe > 36Py_hardyHo.txt
awk 'NR>1 {split($3, a, "/"); total = a[1] + a[2] + a[3]; if (total > 0) {Ho = a[2] / total; print $1, $2, Ho} else {print "Skipping line " NR ": zero denominator" > "/dev/stderr"}}' 3Pt_hardy.hwe > 3Pt_hardyHo.txt
awk 'NR>1 {split($3, a, "/"); total = a[1] + a[2] + a[3]; if (total > 0) {Ho = a[2] / total; print $1, $2, Ho} else {print "Skipping line " NR ": zero denominator" > "/dev/stderr"}}' 3Pyv_hardy.hwe > 3Pyv_hardyHo.txt
