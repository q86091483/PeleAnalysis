case_folder="/scratch/b/bsavard/zisen347/scopingRuns/Burke9_Re6000_2J6_nonreactive"
plts="plt_0{8,9}*"
plt_names=$(eval echo "${case_folder}/${plts}")
for i in ${plt_names}; do
	echo $i
done
srun -N 1 -n 40 KolmogorovScale3d.gnu.MPI.ex input.Kolmogorov infile=${plt_names}
