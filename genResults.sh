#! /bin/bash
# Usage: genResults.sh OMPrun 51
for step in "0.05" "0.005" "0.0005" "0.00005" "0.000005"
do
    # e.g. OMPrun 51 0.05
    { time ./$1 $2 $step > /dev/null ; } 2> tmp.txt
    echo "./$1 $2 $step" >> res.txt
    echo `cat tmp.txt | grep "real"` >> res.txt
done
echo "" >> res.txt
