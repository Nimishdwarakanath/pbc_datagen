directoryList="bonded non-bonded dummy"
counter=0
for directory in ${directoryList}
do
molID=${counter}
data="../${directory}/data.example"
bondLine=`grep -n 'Bonds' $data | awk -F':' '{print $1}'`
angleLine=`grep -n 'Angles' $data | awk -F':' '{print $1}'`
beginLine=$((${bondLine}+2))
endLine=$((${angleLine}-2))
echo "pbc set [exec cat parm] -molid $molID"
echo "pbc box -molid $molID"
echo "topo clearbonds -molid $molID"
sed -n "${beginLine},${endLine}p" $data | awk -v molID=$molID '{printf "topo addbond %i %i -molid %i\n",$3-1,$4-1,molID}'
counter=$((${counter}+1))
done
