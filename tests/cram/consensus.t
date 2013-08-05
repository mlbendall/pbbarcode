  $ export INH5=`python -c "from pbcore import data ; print data.getCmpH5()['cmph5']"`
  $ export INBH51=`python -c "from pbcore import data ; print data.getCmpH5()['bash5s'][0]"`
  $ export INBH52=`python -c "from pbcore import data ; print data.getCmpH5()['bash5s'][1]"`
  $ export BARCODE_FASTA=$TESTDIR/../../etc/barcode.fasta
  $ echo $INBH51 > bas.fofn
  $ echo $INBH52 >> bas.fofn
  $ pbbarcode.py labelZmws $BARCODE_FASTA bas.fofn
  $ pbbarcode.py consensus bas.fofn barcode.fofn
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46[INFO]  [blasr] started.
  2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] [INFO] 2013-08-02T00:28:462013-08-02T00:28:46 [blasr] started. [blasr] started.
  
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
  [INFO] 2013-08-02T00:28:46 [blasr] started.
  [INFO] 2013-08-02T00:28:46 [blasr] ended.
