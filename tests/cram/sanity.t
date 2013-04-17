  $ export INH5=`python -c "from pbcore import data ; print data.getCmpH5()['cmph5']"`
  $ export INBH51=`python -c "from pbcore import data ; print data.getCmpH5()['bash5s'][0]"`
  $ export INBH52=`python -c "from pbcore import data ; print data.getCmpH5()['bash5s'][1]"`
  $ export BARCODE_FASTA=$TESTDIR/../../etc/barcode.fasta
  $ echo $INBH51 > bas.fofn
  $ echo $INBH52 >> bas.fofn
  $ pbbarcode.py labelZMWs $BARCODE_FASTA bas.fofn
  $ echo $?
  0
  $ pbbarcode.py emitFastqs bas.fofn barcode.fofn
  $ echo $?
  0
  $ cp $INH5 ./aligned_reads.cmp.h5         
  $ pbbarcode.py labelAlignments barcode.fofn aligned_reads.cmp.h5
  $ echo $?
  0
  $ pbbarcode.py emitFastqs --subreads bas.fofn barcode.fofn
  $ echo $?
  0
