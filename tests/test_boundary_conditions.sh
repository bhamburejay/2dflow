CMD="mpiexec-mpich-clang17 -n 4"
ehco $CMD
mpiexec-mpich-clang17 -n 4 ./test_boundary_conditions.exe
echo "rank 0" 
cat rank_0_testbc.txt
echo "rank 1" 
cat rank_1_testbc.txt
echo "rank 2" 
cat rank_2_testbc.txt
echo "rank 3" 
cat rank_3_testbc.txt
