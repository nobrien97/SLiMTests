cd /mnt/c/GitHub/SLiMTests/tests/standingVar/analysis/PBS
DATAPATH=/mnt/d/SLiMTests/tests/standingVar

#./addheadertohet.sed $DATAPATH/slim_locusHo.csv > $DATAPATH/slim_locusHo_hdr.csv

HEADERFILE=./addheadertohet.txt

cat $HEADERFILE $DATAPATH/slim_locusHo.csv > $DATAPATH/slim_locusHo_hdr.csv

