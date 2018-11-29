rm simple4.txt cannon4.txt
rm simple16.txt cannon16.txt
rm simple64.txt cannon64.txt
rm simple256.txt cannon256.txt
rm serial.txt

for run in {1..20}
do
	mpirun -n 1 ./mmserial >> serial.txt

	mpirun -n 4 ./mmsimple >> simple4.txt
	mpirun -n 16 ./mmsimple >> simple16.txt
	mpirun -n 64 ./mmsimple >> simple64.txt
#	mpirun -n 256 ./mmsimple >> simple256.txt

	mpirun -n 4 ./mmcannon >> cannon4.txt
	mpirun -n 16 ./mmcannon >> cannon16.txt
	mpirun -n 64 ./mmcannon >> cannon64.txt
#	mpirun -n 256 ./mmcannon >> cannon256.txt
done
