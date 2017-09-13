#!/bin/bash

for i in *.fasta
do

	name=`basename $i .fasta`
	echo $name
	transeq -sequence $i -outseq $name.1.faa -frame 1
	
	stars[1]=`grep -c \* $name.1.faa`
	echo "Frame 1 *s: ${stars[1]}"
	
	if [[ ${stars[1]} -gt 0 ]]
	then
		transeq -sequence $i -outseq $name.2.faa -frame 2
		transeq -sequence $i -outseq $name.3.faa -frame 3
		
		stars[2]=`grep -c \* $name.2.faa`
		stars[3]=`grep -c \* $name.3.faa`
		
		echo "Trying Frame 2: ${stars[2]} and frame 3: ${stars[3]}"

		best=`echo "${stars[@]}" | tr -s ' ' '\n' | awk '{print($0" "NR)}' | sort -g -k1,1 | head -1 | cut -f2 -d' '`
	
		echo "The best translation is $best with ${stars[$best]} *s"
	
		mv $name.$best.faa ../amino_acid/$name.faa
		rm $name.*.faa
	else
		echo "Frame 1 had 0 *s, going with that."
		
		mv $name.1.faa ../amino_acid/$name.faa
	fi
	
done
	
	
	