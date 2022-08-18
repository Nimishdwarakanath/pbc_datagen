rm *.f90_summary summary.f90
for j in MODULE "TYPE (" FUNCTION SUBROUTINE CALL
do
 for i in *.f90
 do
 grep -in "\<${j}" ${i} >> ${i}_summary
 done
done

for i in *.f90
do
sort -n ${i}_summary > t
rm ${i}_summary
echo "!----------------------------------------------------------------------------------" >> summary.f90
echo "! ${i} " >> summary.f90
echo "!----------------------------------------------------------------------------------"  >> summary.f90
cat t >> summary.f90
rm t
done
