#!/bin/bash
RANDOM=$$$(date +%s)

remove() {
  local -n _arr=$1      # underscore-prefixed name to reduce collision likelihood
  local idx=$2
  unset _arr[$idx]      # remove the undesired item
  _arr=( "${_arr[@]}" ) # renumber the indexes
}

containsElement () {
  local e match="$1"
  shift
  for e; do [[ "$e" == "$match" ]] && return 0; done
  return 1
}

if [ $# -lt 3 ]; then
	echo "usage: $0 total_question number_assigned nsets"
	exit 1
fi

if [ $(echo $1 | egrep -c '^[0-9]+$') -eq 0 ]; then
	echo "ERROR: Expected integer for \$1, Got '$1'. Aborting..."
	exit 1
fi
if [ $(echo $2 | egrep -c '^[0-9]+$') -eq 0 ]; then
	echo "ERROR: Expected integer for \$2, Got '$2'. Aborting..."
	exit 1
fi
if [ $(echo $3 | egrep -c '^[0-9]+$') -eq 0 ]; then
	echo "ERROR: Expected integer for \$3, Got '$3'. Aborting..."
	exit 1
fi
if [ $1 -lt $2 ]; then
	echo "ERROR: total_questions($1) less than number_assigned($2). Aborting..."
	exit 1
fi
#populate the tmp.dat file to ensure near uniform sampling
n=$(echo $1 $2 $3 | awk '{printf "%.0f",($2*$3/$1)}')
q=()
for i in $(seq 1 $n)
do
	for j in $(seq 1 $1)
	do
		q=(${q[*]} $j)
	done
done

j=$(($3-1))
for i in $(seq 1 $j)
do
	p=()
	echo "##set $i##"
	k=1
	while [ $k -le $2 ]
	do
		idx=$(($RANDOM % ${#q[@]}))
		num=${q[$idx]}
		if ! containsElement "$num" "${p[@]}"; then
			echo $num
			p=(${p[@]} $num)
			remove q $idx
			k=$(($k + 1))
		fi
	done
done
echo "##set ${3}##"
for i in ${q[@]}
do
	echo $i
done
