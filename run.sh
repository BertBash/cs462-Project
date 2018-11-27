rm simple4.txt
rm simple16.txt
rm simple64.txt
rm simple256.txt


for run in {1..20}
do
	mpirun -n 4 ./mmsimple >> simple4.txt
	mpirun -n 16 ./mmsimple >> simple16.txt
	mpirun -n 64 ./mmsimple >> simple64.txt
	mpirun -n 256 ./mmsimple >> simple256.txt
done
