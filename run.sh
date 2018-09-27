rm x_record.txt
rm *.plt
make
#./ntuplace2 -aux ibm01.aux -util 0.8 -step 0.3 -loop 10 -isparallel 0 -debug 0 -gd 0
# 0.58 0.405
./ntuplace2 -aux ibm01.aux -util 0.8 -step 0.505 -loop 1 -isparallel 0 -debug 0 -checkstep 5 
gprof ntuplace2 > record_ibm01.log
#0.000006507   0.003467539
