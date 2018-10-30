 #!/bin/bash

##Pt 50
hed /msu/data/t3work9/voetberg/vari_dist/mc_dis_50.root\
  -o ../results/mc_dis_lt50.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 100 
hed /msu/data/t3work9/voetberg/vari_dist/mc_dis_100.root\
  -o ../results/mc_dis_50To100.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 150 
hed /msu/data/t3work9/voetberg/vari_dist/mc_dis_150.root\
  -o ../results/mc_dis_100To150.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 200 
hed /msu/data/t3work9/voetberg/vari_dist/mc_dis_200.root\
  -o ../results/mc_dis_150To200.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 250 
hed /msu/data/t3work9/voetberg/vari_dist/mc_dis_250.root\
  -o ../results/mc_dis_200To250.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 300 
hed /msu/data/t3work9/voetberg/vari_dist/mc_dis_300.root\
  -o ../results/mc_dis_250To300.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 350 
hed /msu/data/t3work9/voetberg/vari_dist/mc_dis_350.root\
  -o ../results/mc_dis_300To350.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 400 
hed /msu/data/t3work9/voetberg/vari_dist/mc_dis_400.root\
  -o ../results/mc_dis_350To400.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 400p 
hed /msu/data/t3work9/voetberg/vari_dist/mc_dis_400p.root\
  -o ../results/mc_dis_gt400.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

