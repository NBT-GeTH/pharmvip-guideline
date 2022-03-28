#!/bin/bash
#SBATCH -p compute                  # ระบุ partition หรือประภทเครื่องที่ใช้งาน [Compute/Memory/GPU]
#SBATCH -N 1 --ntasks-per-node=40   # ระบุจำนวนเครื่อง (nodes) และ จำนวน core ต่อ node
#SBATCH -t 04:00:00                 # ระบุเวลาที่ต้องการจองหรือใช้งาน (time limit) สูงสุด โดยมีรูปแบบคือ  ชั่วโมง:นาที:วินาที
#SBATCH -J vcf_inx                   # ระบุชื่อของ Job 

# module purge                        #unload module ทั้งหมด เพราะว่าอาจจะมีการ Load module ไว้ก่อนหน้านั้น
# module load intel                   #load module ที่ต้องการใช้งาน ตัวอย่างนี้คือ intel
# module load python
# python  pass_finder_rs72549303.py                 #สั่งรัน program/executable code ของท่าน
# path="/tarafs/biobank/data/home/ktraipar/sequence/sample/SRR769545_rs72549303.vcf"
path="/tarafs/biobank/data/home/ktraipar/1000genome/HG02012_hg38.vcf"
bash vcf_handle.sh -i $path
echo "done"