#! /bin/bash

for n in $(ls)
  do mv $n "$(echo $n | sed 's/0\./0_/')"; 
done

exit
