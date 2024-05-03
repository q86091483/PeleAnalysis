plt_folder="/scratch/b/bsavard/zisen347/scopingRuns/NUIG_Re4000_2J6_4atm/Level_3"
plts=plt_10580
plt_names=$(eval echo "${case_folder}/${plts}")

echo "Going to process the following files: "
for i in ${plt_folder}/${plts}; do
	echo $i
done
printf "\n"

for i in ${plt_folder}/${plts}; do
	echo "Processing" $i
  srun -N 1 -n 40 RPA3d.gnu.MPI.ex  RPA.input infile=$i
  printf "\n"
done

#plt_folder="/scratch/b/bsavard/zisen347/scopingRuns/NUIG_Re4000_2J6_4atm/Level_2"
#plts=plt_*
#plt_names=$(eval echo "${case_folder}/${plts}")

#echo "Going to process the following files: "
#for i in ${plt_folder}/${plts}; do
#	echo $i
#done
#printf "\n"

#for i in ${plt_folder}/${plts}; do
#	echo "Processing" $i
#  srun -N 2 -n 80 MicroMix3d.gnu.MPI.ex MicroMix.input infile=$i
#  printf "\n"
#done