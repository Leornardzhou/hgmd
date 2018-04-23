# mkdir log out gene
# split -l 100 -a 3 -d HGMD.genelist gene/HGMD.genelist. 

for i in {0..8};do

  for g in `ls gene| grep "0$i[0-9]"`;do

    echo "python craw_gene.py -gl gene/$g -o out/$g 1>log/$g.o.log 2>log/$g.e.log" >> down.$i.sh

  done

done

