 #!/bin/bash
##MC
##Pt 300 
hed /msu/data/t3work9/voetberg/vari_dist/ptratio_dis_300.root\
  -o ../results/ptratio_mc_dis_250To300.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 350 
hed /msu/data/t3work9/voetberg/vari_dist/ptratio_dis_350.root\
  -o ../results/ptratio_mc_dis_300To350.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 400 
hed /msu/data/t3work9/voetberg/vari_dist/ptratio_dis_400.root\
  -o ../results/ptratio_mc_dis_350To400.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Pt 400p 
hed /msu/data/t3work9/voetberg/vari_dist/ptratio_dis_400p.root\
  -o ../results/ptratio_mc_dis_gt400.pdf \
  -e 's/^(Signal|Background)__//' \
     'nl/^(Signal|Background).*/\1/' \
     'n/^Signal__/norm'\
     'n/^Background__/norm'\
  -g 'leg tr'\
  --colors=600 635

##Data
##Pt 300 
hed /msu/data/t3work9/voetberg/vari_dist/ptratio_dis_300.root\
  -o ../results/ptratio_data_dis_250To300.pdf \
  -e 's/^Data_(Signal|Right|Left)__//' \
     'nl/^Data_(Signal|Right|left).*/\1/' \
     'n/^Data_Signal__/norm'\
     'n/^Data_Left__/norm'\
     'n/^Data_Right__/norm'\
  -g 'leg tr'\
  --colors=600 635 65 

##Pt 350 
hed /msu/data/t3work9/voetberg/vari_dist/ptratio_dis_350.root\
  -o ../results/ptratio_data_dis_300To350.pdf \
  -e 's/^Data_(Signal|Right|Left)__//' \
     'nl/^Data_(Signal|Right|left).*/\1/' \
     'n/^Data_Signal__/norm'\
     'n/^Data_Left__/norm'\
     'n/^Data_Right__/norm'\
  -g 'leg tr'\
  --colors=600 635 65


##Pt 400 
hed /msu/data/t3work9/voetberg/vari_dist/ptratio_dis_400.root\
  -o ../results/ptratio_data_dis_350To400.pdf \
  -e 's/^Data_(Signal|Right|Left)__//' \
     'nl/^Data_(Signal|Right|left).*/\1/' \
     'n/^Data_Signal__/norm'\
     'n/^Data_Left__/norm'\
     'n/^Data_Right__/norm'\
  -g 'leg tr'\
  --colors=600 635 65

##Pt 400p 
hed /msu/data/t3work9/voetberg/vari_dist/ptratio_dis_400p.root\
  -o ../results/ptratio_data_dis_gt400.pdf \
  -e 's/^Data_(Signal|Right|Left)__//' \
     'nl/^Data_(Signal|Right|left).*/\1/' \
     'n/^Data_Signal__/norm'\
     'n/^Data_Left__/norm'\
     'n/^Data_Right__/norm'\
  -g 'leg tr'\
  --colors=600 635 65
