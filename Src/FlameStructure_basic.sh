#case_folder="/scratch/b/bsavard/zisen347/scopingRuns/Burke9_Re6000_2J6_nonreactive"
#case_folder="/scratch/b/bsavard/zisen347/scopingRuns/Burke9_Re6000_2J6_nonreactive_1"
case_folder="/scratch/b/bsavard/zisen347/scopingRuns/MicroMix"
plts=plt_04*
plt_names=$(eval echo "${case_folder}/${plts}")

echo "Going to process the following files: "
for i in /scratch/b/bsavard/zisen347/scopingRuns/MicroMix/${plts}; do
	echo $i
done
printf "\n"

for i in /scratch/b/bsavard/zisen347/scopingRuns/MicroMix/${plts}; do
	echo "Procesing" $i
  srun -N 1 -n 40 FlameStructure_basic3d.gnu.MPI.ex FlameStructure_basic.input infile=$i
  printf "\n"
done
