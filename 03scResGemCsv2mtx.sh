cd ./data4test
cat img4test_cp_outlines_scResGem.csv |cut -d, -f2 | sed '1d' | sort -n|uniq > gene.tsv 
cat img4test_cp_outlines_scResGem.csv |cut -d, -f8| sed '1d' | sort -n|uniq > barcodes.tsv 

cat img4test_cp_outlines_scResGem.csv |cut -d, -f2,5,8|sed '1d'|cut -d, -f3|while read barcodesName; do grep -w -n ${barcodesName} barcodes.tsv ; done|cut -d: -f1 > barcodes_index.txt

cat gene.tsv|sed 's/^/bioin&/' |sed 's/$/&plant/' > geneLabeled.tsv 
cat img4test_cp_outlines_scResGem.csv |awk -F,  '{OFS = ","} {print $1,"bioin",$2,"plant",$3,$4,$5,$6,$7,$8}'|sed 's/bioin,geneID,plant/geneID/g'|sed 's/bioin,/bioin/g'|sed 's/,plant/plant/g' > img4test_cp_outlines_scResGemLabeled.csv
cat img4test_cp_outlines_scResGemLabeled.csv |cut -d, -f2,5,8|sed '1d'|cut -d, -f1|while read geneName; do grep -w -n ${geneName} geneLabeled.tsv ; done|cut -d: -f1 > gene_index.txt

cat img4test_cp_outlines_scResGem.csv |cut -d, -f5|sed '1d' > MIDcount.txt; 
paste  gene_index.txt  barcodes_index.txt  MIDcount.txt >matrixPre.mtx

cat gene.tsv |wc -l > headPre.txt 
cat barcodes.tsv |wc -l >> headPre.txt
cat img4test_cp_outlines_scResGem.csv |wc -l >> headPre.txt
sed  ':a;N;$!ba;s/\n/ /g' headPre.txt |sed '1i %'|sed '1i %%MatrixMarket matrix coordinate real general'> head.txt
cat head.txt  matrixPre.mtx > matrix.mtx

paste gene.tsv gene.tsv |awk ' {print $0"\tGene Expression"}' > genes.tsv
cat img4test_cp_outlines_scResGem.csv |cut -d , -f 6,7,8 > img4test_cp_outlines_scResGem_mini.csv

mkdir matrix4test
cp genes.tsv matrix4test
cp barcodes.tsv matrix4test
cp matrix.mtx matrix4test

rm geneLabeled.tsv
rm headPre.txt
rm head.txt
rm img4test_cp_outlines_scResGemLabeled.csv
img4test_cp_outlines_scResGem_mini.csv
rm matrixPre.mtx
rm MIDcount.txt