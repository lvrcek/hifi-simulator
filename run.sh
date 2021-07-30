python -u distributions.py ../hifi-data/escherichia/SRR10971019_subreads.fastq --length --qscore --accuracy &> ecoli.log &
python -u distributions.py ../hifi-data/drosophila/SRR12473480_subreads.fastq.gz --length --qscore --accuracy &> drosophila.log &
python -u distributions.py ../hifi-data/mosquito/SRR12121585_subreads.fastq --length --qscore --accuracy &> mosquito.log &
python -u distributions.py ../hifi-data/human/SRX5633451.fastq.gz --length --qscore --accuracy &> human_X.log &
python -u distributions.py ../hifi-data/human/20200511.C004C.ccs.fastq.gz --length --qscore --accuracy &> human_C.log &
python -u distributions.py ../hifi-data/human/20200513.I002C.ccs.fastq.gz --length --qscore --accuracy &> human_I.log &
python -u distributions.py ../hifi-data/human/SRR11292120_3_subreads.fastq.gz --length --qscore --accuracy &> human_R.log &

