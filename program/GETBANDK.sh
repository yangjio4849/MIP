#!/bin/bash
# need file: POSCAR
# need lib: phonopy-1.12.0
#
#    It is based on Setyawan, W., & Curtarolo, S. (2010).
#    High-throughput electronic band structure calculations:
#    Challenges and tools. Computational Materials Science,
#    49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010


pai=3.14159265358979323846
abs(){ echo ${1#-};}
#cp ../../SCF/test_scan/CONTCAR POSCAR
phonopy --symmetry --tolerance=0.1 >symmetry.txt
spacegroup=`grep "space_group_number:" symmetry.txt | awk '{print $2}'`   ###########find spacegroup number
symbol=`grep "space_group_type:" symmetry.txt |awk '{print $2}' | cut -c2-2`  ######find symbol

######start to divide seven crystal system by the spacegroup number
if (( $spacegroup >=1 & $spacegroup <=2 ));then
  latticetype="triclinic"
  elif (( $spacegroup >=3 & $spacegroup <=15 ));then
  latticetype="monoclinic"
  elif (( $spacegroup >=16 & $spacegroup <=74 ));then
  latticetype="orthorhombic"
  elif (( $spacegroup >=75 & $spacegroup <=142 ));then
  latticetype="tetragonal"
  # elif (( $spacegroup >=143 & $spacegroup <=167 ));then    # recode 20188/11/05
  elif (( $spacegroup ==146 || $spacegroup == 148 || $spacegroup == 155 || $spacegroup == 160 || $spacegroup == 161 || $spacegroup == 166 || $spacegroup == 167 ));then
  latticetype="rhombohedral"
  # elif (( $spacegroup >=168 & $spacegroup <=194 ));then
  elif (( $spacegroup >=143 & $spacegroup <=194 ));then
  latticetype="hexagonal"
  elif (( $spacegroup >=195 & $spacegroup <=230 ));then
  latticetype="cubic"
  else
  echo "lattice type error"
  exit
fi
##########finished the divide

######get the lattice parameter using the lattice_parameter program
cp POSCAR POSCAR_tmp   #
cp BPOSCAR POSCAR      #  

for ((i=1;i<=6;i++))
do
  lattice[$i]=`lattice_parameters | awk -v j=$i '{print $j}'`
done 

mv POSCAR_tmp POSCAR   #
##########start to classify different KPOINTS type
if [ $latticetype == "cubic" ];then
  if [ $symbol == P ];then
  typemode="CUB"
  elif [ $symbol == F ];then
  typemode="FCC"
  elif [ $symbol == I ];then
  typemode="BCC"
  else
  echo "(Unexpected value for symbol: $symbol)"
  exit
  fi
elif [ $latticetype == "tetragonal" ];then
  if [ $symbol == P ];then
  typemode="TET"
  elif [ $symbol == I ];then
    # conventional lattice  a a c , a > c BCT1 ;a<c BCT2 . a=${lattice[1]} b=${lattice[2]} c=${lattice[3]}
    cp POSCAR POSCAR_tmp   #
    cp BPOSCAR POSCAR      #  

    for ((i=1;i<=6;i++))
    do
       lattice[$i]=`lattice_parameters | awk -v j=$i '{print $j}'`
    done

    mv POSCAR_tmp POSCAR   #
## add 18/11/15


    if (( `echo ${lattice[1]} == ${lattice[2]}  | bc` )); then      
        if (( `echo "${lattice[1]} > ${lattice[3]}" | bc ` ));then
        typemode="BCT1"
        else
        typemode="BCT2"
        fi
    elif (( `echo ${lattice[1]} == ${lattice[3]} | bc ` )); then      
        if (( `echo "${lattice[1]} > ${lattice[2]}" | bc ` ));then
        typemode="BCT1"
        else
        typemode="BCT2"
        fi
    elif (( `echo ${lattice[2]} == ${lattice[3]}  | bc ` )); then      
        if (( `echo "${lattice[1]} < ${lattice[3]}" | bc ` ));then
        typemode="BCT1"
        else
        typemode="BCT2"
        fi
    fi
  else
  echo "(Unexpected value for symbol: $symbol)"
  exit
  fi
elif [ $latticetype == "orthorhombic" ];then

# conventional lattice  a < b < c ,  a=${lattice[1]} b=${lattice[2]} c=${lattice[3]}
#
after_sort=` echo -e "${lattice[1]}\n${lattice[2]}\n${lattice[3]}" | sort -n `
for ((i=1;i<=3;i++))
do
  lattice[$i]=`echo $after_sort | awk -v j=$i '{print $j}'`
done

#echo   ${lattice[*]}  # todo

  if [ $symbol == P ];then
  typemode="ORC"
  elif [ $symbol == F ];then
    if (( ` echo " scale=7; 1.0000000 / ${lattice[1]} / ${lattice[1]} > 1.0000000 / ${lattice[2]} / ${lattice[2]} + 1.0000000 / ${lattice[3]} / ${lattice[3]} " | bc  ` ));then
    typemode="ORCF1"
    elif (( `echo " scale=7; 1.0000000 / ${lattice[1]} / ${lattice[1]} < 1.0000000 / ${lattice[2]} / ${lattice[2]} + 1.0000000 / ${lattice[3]} / ${lattice[3]} " | bc ` ));then
    typemode="ORCF2"
    else
    typemode="ORCF3"
    fi
  elif [ $symbol == I ];then
  typemode="ORCI"
#echo ${lattice[*]}  # todo
  elif [ $symbol == C ] || [  $symbol ==A ] ;then
  typemode="ORCC"

else
  echo "(Unexpected value for symbol: $symbol)"
  exit
  fi
elif [ $latticetype == "hexagonal" ];then
  typemode="HEX"
elif [ $latticetype == "rhombohedral" ];then

# alpha = self._prim.lattice.lengths_and_angles[1][0] # PPOSCAR_lattice
cp POSCAR POSCAR_tmp   #
cp PPOSCAR POSCAR      #  

for ((i=1;i<=6;i++))
do
  lattice[$i]=`lattice_parameters | awk -v j=$i '{print $j}'`
done 

mv POSCAR_tmp POSCAR   #

  if (( `echo "${lattice[4]} < 90" | bc ` ));then
  typemode="RHL1"
  else
  typemode="RHL2"
  fi
elif [ $latticetype == "monoclinic" ];then

cp POSCAR POSCAR_tmp   #
cp BPOSCAR POSCAR      #

for ((i=1;i<=6;i++))
do
  lattice[$i]=`lattice_parameters | awk -v j=$i '{print $j}'`
done
mv POSCAR_tmp POSCAR #

#echo 'lattice BP' ${lattice[*]} # todo

# swapped a b
tmp=${lattice[1]}
lattice[1]=${lattice[2]}
lattice[2]=${tmp}

#after_sort=` echo -e "${lattice[1]}\n${lattice[2]}\n${lattice[3]}" | sort -n `
# # a b <= c
#for ((i=1;i<=3;i++))
#do
#  lattice[$i]=`echo $after_sort | awk -v j=$i '{print $j}'`
#done

# alpha < 90, beta=gama=90
for ((i=3;i<=6;i++))
do
    if [[ $(echo " ${lattice[$i]} != 90 " | bc) == 1 ]]
    then
        if [[ $( echo " ${lattice[$i]} > 90 " | bc ) == 1 ]]
        then
        alpha=$( echo " 180 - ${lattice[$i]}  " | bc )
        else
        alpha=${lattice[$i]}
        fi
    fi

done

#echo 'after sort '
#echo ${lattice[*]}
#echo 'alpha' $alpha # todo
# Pass value
lattice[4]=$alpha
lattice[6]=$gamma

  if [ $symbol == P ];then
  typemode="MCL"
  elif [ $symbol == C ];then
################################################# add by MJ

# # get a b c
# a1=` awk 'NR==3{print $1}' BPOSCAR `
# b1=` awk 'NR==3{print $2}' BPOSCAR `
# c1=` awk 'NR==3{print $3}' BPOSCAR `
# a2=` awk 'NR==4{print $1}' BPOSCAR `
# b2=` awk 'NR==4{print $2}' BPOSCAR `
# c2=` awk 'NR==4{print $3}' BPOSCAR `
# a3=` awk 'NR==5{print $1}' BPOSCAR `
# b3=` awk 'NR==5{print $2}' BPOSCAR `
# c3=` awk 'NR==5{print $3}' BPOSCAR `

a=${lattice[1]}
b=${lattice[2]}
c=${lattice[3]}

#Conventional lattice
#a1 = ( a, 0, 0)
#a2 = (0, b, 0)
#a3 = (0, c cos alpha, c sin alpha)
# ----------------------------------------------
#Primitive lattice
#a1 = ( a/2, b/2, 0)                        # a11, a12, a13
#a2 = (-a/2, b/2, 0)                        # a21, a22, a23
#a3 = (0, c cos alpha, c sin alpha)   # a31, a32, a33
# a11, a12, a13
# a21, a22, a23
# a31, a32, a33

a11=$( echo  "scale=7; $a / 2 " | bc )
a12=$( echo  "scale=7; $b/ 2 " | bc )
a13=0
a21=$( echo  "scale=7; -$a / 2 " | bc )
a22=$( echo  "scale=7; $b / 2 " | bc )
a23=0
a31=0
a32=$( echo  "scale=7; $c * c ( $alpha * $pai / 180) " | bc -l )
a33=$( echo  "scale=7; $c * s ( $alpha * $pai / 180) " | bc -l )

#echo 'prim' # todo
#echo $a11, $a12, $a13
#echo $a21, $a22, $a23
#echo $a31, $a32, $a33
#echo '-------'

a1=$(echo $a11 $a12 $a13 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
a2=$(echo $a21 $a22 $a23 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
a3=$(echo $a31 $a32 $a33 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
#b1 * b2
a2_a3=$( echo  "scale=7; $a21*$a31 + $a22 *$a32 + $a23 *$a33" | bc )
a1_a3=$( echo  "scale=7; $a11*$a31 + $a12 *$a32 + $a13 *$a33" | bc )
a1_a2=$( echo  "scale=7; $a11*$a21 + $a12 *$a22 + $a13 *$a23" | bc )


x=$( echo  "scale=7; $a2_a3/($a2*$a3)" | bc )
y=$( echo  "scale=7; $a1_a3/($a1*$a3)" | bc )
z=$( echo  "scale=7; $a1_a2/($a1*$a2)" | bc )
# arccos()
#alpha=$( awk "BEGIN{x=$x ;print atan2(sqrt(1-(x ^2 )), x)/3.1415927*180}" )
#beta=$( awk "BEGIN{x=$y ;print atan2(sqrt(1-(x ^2 )), x)/3.1415927*180}" )
#gamma=$( awk "BEGIN{x=$z ;print atan2(sqrt(1-(x ^2 )), x)/3.1415927*180}" )

#echo 'a1,a2,a3,al,be,ga'# todo
#echo $a1, $a2, $a3, $alpha, $beta, $gamma
#


#  A ^(-1) = 1/|A|  A*
# |A|
detA=$( echo "$a11 * $a22 * $a33 + $a12 * $a23 * $a31 + $a21 * $a32 * $a13 - $a13 * $a22 * $a31 - $a12 * $a21 * $a33 - $a11 * $a23 * $a32 "  | bc )
#echo 'detA' $detA #todo

#A*
# c11, c12, c13
# c21, c22, c23
# c31, c32, c33
c11=$( echo  "scale=7; $a22 * $a33 - $a32 * $a23 " | bc )
c12=$( echo  "scale=7; -($a12 * $a33 - $a32 * $a13) " | bc )
c13=$( echo  "scale=7; $a12 * $a23 - $a22 * $a13 " | bc )
c21=$( echo  "scale=7; -($a21 * $a33 - $a31 * $a23) " | bc )
c22=$( echo  "scale=7; $a11 * $a33 - $a31 * $a13 " | bc )
c23=$( echo  "scale=7; -($a11 * $a23 - $a21 * $a13) " | bc )
c31=$( echo  "scale=7; $a21 * $a32 - $a31 * $a22 " | bc )
c32=$( echo  "scale=7; -($a11 * $a32 - $a31 * $a12) " | bc )
c33=$( echo  "scale=7; $a11 * $a22 - $a21 * $a12 " | bc )
#echo 'c: A*' # todo
#echo $c11, $c12, $c13
#echo $c21, $c22, $c23
#echo $c31, $c32, $c33
#echo '-------'


# A^-1
# d11, d12, d13
# d21, d22, d23
# d31, d32, d33
d11=$( echo  "scale=7; $c11 / $detA " | bc )
d12=$( echo  "scale=7; $c12 / $detA " | bc )
d13=$( echo  "scale=7; $c13 / $detA " | bc )
d21=$( echo  "scale=7; $c21 / $detA " | bc )
d22=$( echo  "scale=7; $c22 / $detA " | bc )
d23=$( echo  "scale=7; $c23 / $detA " | bc )
d31=$( echo  "scale=7; $c31 / $detA " | bc )
d32=$( echo  "scale=7; $c32 / $detA " | bc )
d33=$( echo  "scale=7; $c33 / $detA " | bc )

#echo 'd : A^-1 ' # todo
#echo $d11, $d12, $d13
#echo $d21, $d22, $d23
#echo $d31, $d32, $d33
#echo '-------'

# A^-1 .T
#  b11, b12, b13 = { d11, d21, d31
#  b21, b22, b23 = { d12, d22, d32
#  b31, b32, b33 = { d13, d23, d33
b11=$d11
b12=$d21
b13=$d31
b21=$d12
b22=$d22
b23=$d32
b31=$d13
b32=$d23
b33=$d33

#echo 'b : A^-1 .T ' # todo
#echo $b11, $b12, $b13
#echo $b21, $b22, $b23
#echo $b31, $b32, $b33
#echo '-------'




# kgamma : between b1 and b2
#|b1|
b1=$(echo $b11 $b12 $b13 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
b2=$(echo $b21 $b22 $b23 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
b3=$(echo $b31 $b32 $b33 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
#b1 * b2
b2_b3=$( echo  "scale=7; $b21*$b31 + $b22 *$b32 + $b23 *$b33" | bc )
b1_b3=$( echo  "scale=7; $b11*$b31 + $b12 *$b32 + $b13 *$b33" | bc )
b1_b2=$( echo  "scale=7; $b11*$b21 + $b12 *$b22 + $b13 *$b23" | bc )


x=$( echo  "scale=7; $b2_b3/($b2*$b3)" | bc )
y=$( echo  "scale=7; $b1_b3/($b1*$b3)" | bc )
z=$( echo  "scale=7; $b1_b2/($b1*$b2)" | bc )
# arccos()
kalpha=$( awk "BEGIN{x=$x ;print atan2(sqrt(1-(x ^2 )), x)/$pai*180}" )
kbeta=$( awk "BEGIN{x=$y ;print atan2(sqrt(1-(x ^2 )), x)/$pai*180}" )
kgamma=$( awk "BEGIN{x=$z ;print atan2(sqrt(1-(x ^2 )), x)/$pai*180}" )

#echo 'b1,b2,b3,kal,kbe,kga'# todo
#echo $b1, $b2, $b3, $kalpha, $kbeta, $kgamma
################################################## end
# Pass value
lattice[4]=$alpha
lattice[6]=$kgamma

#echo 'if' ${lattice[*]} # todo :del

    if (( `echo "${lattice[6]} > 90" | bc ` ));then
    typemode="MCLC1"
    elif (( `echo "${lattice[6]} < 90" | bc ` ));then
      if (( `echo " scale=7; ${lattice[2]} * c ( ${lattice[4]} * $pai / 180) / ${lattice[3]} + ${lattice[2]} ^ 2 * (s ( ${lattice[4]} * $pai / 180) ^ 2) / ${lattice[1]} ^ 2 < 1" | bc -l` ));then
      typemode="MCLC3"
      elif (( `echo "scale=7; ${lattice[2]} * c ( ${lattice[4]} * $pai / 180) / ${lattice[3]} + ${lattice[2]} ^ 2 * (s ( ${lattice[4]} * $pai / 180) ^ 2) / ${lattice[1]} ^ 2 > 1" | bc -l`));then
      typemode="MCLC5"
      else
      typemode="MCLC4"
      fi
    else
    typemode="MCLC2"
    fi
  else
  echo "(Unexpected value for symbol: $symbol)"
  exit
  fi
elif [ $latticetype == "triclinic" ];then
#      kalpha = self._prim_rec.lengths_and_angles[1][0]
#      kbeta = self._prim_rec.lengths_and_angles[1][1]
#      kgamma = self._prim_rec.lengths_and_angles[1][2]

##########
# lattice
#a1 = ( a, 0, 0)
#a2 = (b cos beta, b sin beta, 0)
#a3 = (-------------------------------------)
# # get a b c
a11=` awk 'NR==3{print $1}' BPOSCAR `
a12=` awk 'NR==3{print $2}' BPOSCAR `
a13=` awk 'NR==3{print $3}' BPOSCAR `
a21=` awk 'NR==4{print $1}' BPOSCAR `
a22=` awk 'NR==4{print $2}' BPOSCAR `
a23=` awk 'NR==4{print $3}' BPOSCAR `
a31=` awk 'NR==5{print $1}' BPOSCAR `
a32=` awk 'NR==5{print $2}' BPOSCAR `
a33=` awk 'NR==5{print $3}' BPOSCAR `

#head -6 PPOSCAR
#echo 'prim' # todo
#echo $a11, $a12, $a13
#echo $a21, $a22, $a23
#echo $a31, $a32, $a33
#echo '-------'

a1=$(echo $a11 $a12 $a13 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
a2=$(echo $a21 $a22 $a23 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
a3=$(echo $a31 $a32 $a33 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
#b1 * b2
a2_a3=$( echo  "scale=7; $a21*$a31 + $a22 *$a32 + $a23 *$a33" | bc )
a1_a3=$( echo  "scale=7; $a11*$a31 + $a12 *$a32 + $a13 *$a33" | bc )
a1_a2=$( echo  "scale=7; $a11*$a21 + $a12 *$a22 + $a13 *$a23" | bc )


x=$( echo  "scale=7; $a2_a3/($a2*$a3)" | bc )
y=$( echo  "scale=7; $a1_a3/($a1*$a3)" | bc )
z=$( echo  "scale=7; $a1_a2/($a1*$a2)" | bc )
# arccos()
alpha=$( awk "BEGIN{x=$x ;print atan2(sqrt(1-(x ^2 )), x)/$pai*180}" )
beta=$( awk "BEGIN{x=$y ;print atan2(sqrt(1-(x ^2 )), x)/$pai*180}" )
gamma=$( awk "BEGIN{x=$z ;print atan2(sqrt(1-(x ^2 )), x)/$pai*180}" )

#echo 'a1,a2,a3,al,be,ga'# todo
#echo $a1, $a2, $a3, $alpha, $beta, $gamma



#  A ^(-1) = 1/|A|  A*
# |A|
detA=$( echo "$a11 * $a22 * $a33 + $a12 * $a23 * $a31 + $a21 * $a32 * $a13 - $a13 * $a22 * $a31 - $a12 * $a21 * $a33 - $a11 * $a23 * $a32 "  | bc )
#echo 'detA' $detA #todo

#A*
# c11, c12, c13
# c21, c22, c23
# c31, c32, c33
c11=$( echo  "scale=7; $a22 * $a33 - $a32 * $a23 " | bc )
c12=$( echo  "scale=7; -($a12 * $a33 - $a32 * $a13) " | bc )
c13=$( echo  "scale=7; $a12 * $a23 - $a22 * $a13 " | bc )
c21=$( echo  "scale=7; -($a21 * $a33 - $a31 * $a23) " | bc )
c22=$( echo  "scale=7; $a11 * $a33 - $a31 * $a13 " | bc )
c23=$( echo  "scale=7; -($a11 * $a23 - $a21 * $a13) " | bc )
c31=$( echo  "scale=7; $a21 * $a32 - $a31 * $a22 " | bc )
c32=$( echo  "scale=7; -($a11 * $a32 - $a31 * $a12) " | bc )
c33=$( echo  "scale=7; $a11 * $a22 - $a21 * $a12 " | bc )
#echo 'c: A*' # todo
#echo $c11, $c12, $c13
#echo $c21, $c22, $c23
#echo $c31, $c32, $c33
#echo '-------'


# A^-1
# d11, d12, d13
# d21, d22, d23
# d31, d32, d33
d11=$( echo  "scale=7; $c11 / $detA " | bc )
d12=$( echo  "scale=7; $c12 / $detA " | bc )
d13=$( echo  "scale=7; $c13 / $detA " | bc )
d21=$( echo  "scale=7; $c21 / $detA " | bc )
d22=$( echo  "scale=7; $c22 / $detA " | bc )
d23=$( echo  "scale=7; $c23 / $detA " | bc )
d31=$( echo  "scale=7; $c31 / $detA " | bc )
d32=$( echo  "scale=7; $c32 / $detA " | bc )
d33=$( echo  "scale=7; $c33 / $detA " | bc )

#echo 'd : A^-1 ' # todo
#echo $d11, $d12, $d13
#echo $d21, $d22, $d23
#echo $d31, $d32, $d33
#echo '-------'

# A^-1 .T
#  b11, b12, b13 = { d11, d21, d31
#  b21, b22, b23 = { d12, d22, d32
#  b31, b32, b33 = { d13, d23, d33
b11=$d11
b12=$d21
b13=$d31
b21=$d12
b22=$d22
b23=$d32
b31=$d13
b32=$d23
b33=$d33

#echo 'b : A^-1 .T ' # todo
#echo $b11, $b12, $b13
#echo $b21, $b22, $b23
#echo $b31, $b32, $b33
#echo '-------'




# kgamma : between b1 and b2
#|b1|
b1=$(echo $b11 $b12 $b13 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
b2=$(echo $b21 $b22 $b23 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
b3=$(echo $b31 $b32 $b33 | awk '{printf("%.7f",  sqrt($1^2 + $2^2 + $3^2 ) )}')
#b1 * b2
b2_b3=$( echo  "scale=7; $b21*$b31 + $b22 *$b32 + $b23 *$b33" | bc )
b1_b3=$( echo  "scale=7; $b11*$b31 + $b12 *$b32 + $b13 *$b33" | bc )
b1_b2=$( echo  "scale=7; $b11*$b21 + $b12 *$b22 + $b13 *$b23" | bc )


x=$( echo  "scale=7; $b2_b3/($b2*$b3)" | bc )
y=$( echo  "scale=7; $b1_b3/($b1*$b3)" | bc )
z=$( echo  "scale=7; $b1_b2/($b1*$b2)" | bc )
# arccos()
kalpha=$( awk "BEGIN{x=$x ;print atan2(sqrt(1-(x ^2 )), x)/$pai*180}" )
kbeta=$( awk "BEGIN{x=$y ;print atan2(sqrt(1-(x ^2 )), x)/$pai*180}" )
kgamma=$( awk "BEGIN{x=$z ;print atan2(sqrt(1-(x ^2 )), x)/$pai*180}" )

#echo 'b1,b2,b3,kal,kbe,kga'# todo
#echo $b1, $b2, $b3, $kalpha, $kbeta, $kgamma
################################################## end
## Pass value
lattice[4]=$kalpha
lattice[5]=$kbeta
lattice[6]=$kgamma

##########

  if (( `echo " ${lattice[4]} > 90" | bc` & `echo "${lattice[5]} > 90" | bc` & `echo "${lattice[6]} > 90" | bc` ));then
  typemode="TRIA"
  elif (( `echo "${lattice[4]} < 90" | bc` & `echo "${lattice[5]} < 90" | bc` & `echo "${lattice[6]} < 90" | bc` ));then
  typemode="TRIB"
  elif (( `echo "${lattice[4]} > 90" | bc` & `echo "${lattice[5]} > 90" | bc` & `echo "${lattice[6]} >= 89.99999995" | bc` & `echo "${lattice[6]} < 90.00000005" | bc` ));then
  typemode="TRIA"
  elif (( `echo "${lattice[4]} < 90" | bc` & `echo "${lattice[5]} < 90" | bc` & `echo "${lattice[6]} >= 89.99999995" | bc` & `echo "${lattice[6]} < 90.00000005" | bc` ));then
  typemode="TRIB"
### add 18/11/14
  else
################## if error
  typemode="TRIB"

  fi
else
echo "Uniknown lattice type $latticetype"
exit
fi
rm -rf BPOSCAR  PPOSCAR  symmetry.txt
#
#echo $typemode     # todo
#echo ${lattice[*]}
##################END by lixin
##################${lattice[1-6]} means a b c alpha beta gamma
##################typemode means different KPOINTS mode


################## cat KPOINTS.band part 
################## read kmesh density from command line 

name=$typemode
# read -p 'please input the kmesh of band: '  mesh_in_put

if [[ $# -gt 0 ]]
then 
    # kmesh=$mesh_in_put
    kmesh=$1
else
    kmesh=40
fi
##### pi 
a=${lattice[1]}
b=${lattice[2]}
c=${lattice[3]}
alpha=$(echo ${lattice[4]} $pai | awk '{printf("%.11f", ( $1 / 180 * $2 ))}')    #pi
beta=$(echo ${lattice[5]} $pai | awk '{printf("%.11f", ( $1 / 180 * $2 ))}')

#echo 'jisuancanshu:' $a $b $c $alpha $beta # todo

case $name in
  BCT1) 

eta=$(echo $c $a | awk '{printf("%.11f", ((1 + $1 ** 2 / $2 ** 2 )/ 4.0))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.0 0.0 0.5 ! X

0.0 0.0 0.5 ! X
-0.5 0.5 0.5 ! M

-0.5 0.5 0.5 ! M
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
$eta $eta -$eta ! Z

$eta $eta -$eta ! Z
0.25 0.25 0.25 ! P

0.25 0.25 0.25 ! P
0.0 0.5 0.0 ! N

0.0 0.5 0.0 ! N
-$eta $eta_1 $eta ! Z_1

-$eta $eta_1 $eta ! Z_1
-0.5 0.5 0.5 ! M

0.0 0.0 0.5 ! X
0.25 0.25 0.25 ! P

EOF
;;
  BCT2) 

eta=$(echo $a $c | awk '{printf("%.11f", ((1 + $1 ** 2 / $2 ** 2 )/ 4.0))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1))}')
zeta=$(echo $a $c | awk '{printf("%.11f", ($1 ** 2/(2 * $2 ** 2)))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.0 0.0 0.5 ! X

0.0 0.0 0.5 ! X
-$zeta $zeta 0.5 ! Y

-$zeta $zeta 0.5 ! Y
-$eta $eta $eta ! \Sigma

-$eta $eta $eta ! \Sigma
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 0.5 -0.5 ! Z

0.5 0.5 -0.5 ! Z
$eta $eta_1 -$eta ! \Sigma_1

$eta $eta_1 -$eta ! \Sigma_1
0.0 0.5 0.0 ! N

0.0 0.5 0.0 ! N
0.25 0.25 0.25 ! P

0.25 0.25 0.25 ! P
0.5 0.5 -$zeta ! Y_1

0.5 0.5 -$zeta ! Y_1
0.5 0.5 -0.5 ! Z

0.0 0.0 0.5 ! X
0.25 0.25 0.25 ! P

EOF
;;
  MCL) 

eta=$(echo $b $c $alpha | awk '{printf("%.11f", ((1 - $1 * cos($3) / $2) / (2 * sin($3) ** 2)))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1))}')
nu=$(echo $b $c $alpha $eta | awk '{printf("%.11f", (0.5 - $4 * $2 * cos($3) / $1))}')
nu_1=$(echo $nu | awk '{printf("%.11f", (1 - $1))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.0 0.0 0.5 ! Y

0.0 0.0 0.5 ! Y
0.0 $eta $nu_1 ! H

0.0 $eta $nu_1 ! H
0.0 0.5 0.5 ! C

0.0 0.5 0.5 ! C
0.5 0.5 0.5 ! E

0.5 0.5 0.5 ! E
0.5 $eta_1 $nu ! M_1

0.5 $eta_1 $nu ! M_1
0.5 0.5 0.0 ! A

0.5 0.5 0.0 ! A
0.0 0.5 0.0 ! X

0.0 0.5 0.0 ! X
0.0 $eta_1 $nu ! H_1

0.5 $eta $nu_1 ! M
0.5 0.0 0.5 ! D

0.5 0.0 0.5 ! D
0.5 0.0 0.0 ! Z

0.0 0.0 0.5 ! Y
0.5 0.0 0.5 ! D

EOF
;;
  MCLC1)

zeta=$(echo $b $c $alpha | awk '{printf("%.11f", ((2 - $1 * cos($3) / $2) / (4 * sin($3) ** 2)))}')
zeta_1=$(echo $zeta | awk '{printf("%.11f", (1 - $1))}')
eta=$(echo $b $c $alpha $zeta | awk '{printf("%.11f", (0.5 + 2 * $4 * $2 * cos($3) / $1))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1))}')
psi=$(echo $a $b $alpha | awk '{printf("%.11f", (0.75 - $1 ** 2 / (4 * $2 ** 2 * sin($3) ** 2)))}')
psi_1=$(echo $psi | awk '{printf("%.11f", (1 - $1))}')
psi_2=$(echo $psi | awk '{printf("%.11f", ($1 - 1))}')
phi=$(echo $b $c $alpha $psi | awk '{printf("%.11f", ($4 + (0.75 - $4) * $1 * cos($3) / $2))}')
phi_1=$(echo $phi | awk '{printf("%.11f", (1 - $1))}')
phi_2=$(echo $phi | awk '{printf("%.11f", ($1 - 1))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.5 0.5 0.0 ! Y

0.5 0.5 0.0 ! Y
$zeta_1 $zeta_1 $eta_1 ! F

$zeta_1 $zeta_1 $eta_1 ! F
0.5 0.5 0.5 ! L

0.5 0.5 0.5 ! L
$phi $phi_1 0.5 ! I

$phi_1 $phi_2 0.5 ! I_1
0.0 0.0 0.5 ! Z

0.0 0.0 0.5 ! Z
$zeta $zeta $eta ! F_1

0.5 0.5 0.0 ! Y
$psi $psi_1 0.0 ! X_1

$psi_1 $psi_2 0.0 ! X
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.0 ! N

0.5 0.0 0.5 ! M
0.0 0.0 0.0 ! \Gamma

EOF
;;
  MCLC2)

zeta=$(echo $b $c $alpha | awk '{printf("%.11f", ((2 - $1 * cos($3) / $2) / (4 * sin($3) ** 2)))}')
zeta_1=$(echo $zeta | awk '{printf("%.11f", (1 - $1))}')
eta=$(echo $b $c $alpha $zeta | awk '{printf("%.11f", (0.5 + 2 * $4 * $2 * cos($3) / $1))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1))}')
psi=$(echo $a $b $alpha | awk '{printf("%.11f", (0.75 - $1 ** 2 / (4 * $2 ** 2 * sin($3) ** 2)))}')
psi_1=$(echo $psi | awk '{printf("%.11f", (1 - $1))}')
psi_2=$(echo $psi | awk '{printf("%.11f", ($1 - 1))}')
phi=$(echo $b $c $alpha $psi | awk '{printf("%.11f", ($4 + (0.75 - $4) * $1 * cos($3) / $2))}')
phi_1=$(echo $phi | awk '{printf("%.11f", (1 - $1))}')
phi_2=$(echo $phi | awk '{printf("%.11f", ($1 - 1))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.5 0.5 0.0 ! Y

0.5 0.5 0.0 ! Y
$zeta_1 $zeta_1 $eta_1 ! F

$zeta_1 $zeta_1 $eta_1 ! F
0.5 0.5 0.5 ! L

0.5 0.5 0.5 ! L
$phi $phi_1 0.5 ! I

$phi_1 $phi_2 0.5 ! I_1
0.0 0.0 0.5 ! Z

0.0 0.0 0.5 ! Z
$zeta $zeta $eta ! F_1

0.5 0.0 0.0 ! N
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.5 ! M

EOF
;;
  MCLC3)

mu=$(echo $a $b | awk '{printf("%.11f", ((1 + $2 ** 2 / $1 ** 2) / 4.0))}')
mu_1=$(echo $mu | awk '{printf("%.11f", (1 - $1))}')
delta=$(echo $a $b $c $alpha | awk '{printf("%.11f", ( $2 * $3 * cos($4) / (2 * $1 ** 2)))}')
zeta=$(echo $b $c $alpha $mu | awk '{printf("%.11f", ($4 - 0.25 + (1 - $1 * cos($3) / $2)/ (4 * sin($3) ** 2)))}')
zeta_1=$(echo $zeta | awk '{printf("%.11f", (1 - $1))}')
eta=$(echo $b $c $alpha $zeta | awk '{printf("%.11f", (0.5 + 2 * $4 * $2 * cos($3) / $1))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1))}')
psi=$(echo $eta $delta | awk '{printf("%.11f", ($1 - 2 * $2))}')
psi_1=$(echo $psi | awk '{printf("%.11f", (1 - $1))}')
psi_2=$(echo $psi | awk '{printf("%.11f", ($1 - 1))}')
phi=$(echo $zeta $mu | awk '{printf("%.11f", (1 + $1 - 2 * $2))}')
phi_1=$(echo $phi | awk '{printf("%.11f", (1 - $1))}')
phi_2=$(echo $phi | awk '{printf("%.11f", ($1 - 1))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
$mu $mu $delta ! Y

$mu $mu $delta ! Y
$phi_1 $phi_1 $psi_1 ! F

$phi_1 $phi_1 $psi_1 ! F
$zeta $zeta $eta ! H

$zeta $zeta $eta ! H
0.0 0.0 0.5 ! Z

0.0 0.0 0.5 ! Z
0.5 -0.5 0.5 ! I

0.5 -0.5 0.5 ! I
$phi $phi_2 $psi ! F_1

$zeta_1 -$zeta $eta_1 ! H_1
$mu_1 -$mu -$delta ! Y_1

$mu_1 -$mu -$delta ! Y_1
0.5 -0.5 0.0 ! X

0.5 -0.5 0.0 ! X
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.0 ! N

0.5 0.0 0.5 ! M
0.0 0.0 0.0 ! \Gamma

EOF
;;
  MCLC4)

mu=$(echo $a $b | awk '{printf("%.11f", ((1 + $2 ** 2 / $1 ** 2) / 4.0))}')
mu_1=$(echo $mu | awk '{printf("%.11f", (1 - $1))}')
delta=$(echo $a $b $c $alpha | awk '{printf("%.11f", ( $2 * $3 * cos($4) / (2 * $1 ** 2)))}')
zeta=$(echo $b $c $alpha $mu | awk '{printf("%.11f", ($4 - 0.25 + (1 - $1 * cos($3) / $2)/ (4 * sin($3) ** 2)))}')
zeta_1=$(echo $zeta | awk '{printf("%.11f", (1 - $1))}')
eta=$(echo $b $c $alpha $zeta | awk '{printf("%.11f", (0.5 + 2 * $4 * $2 * cos($3) / $1))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1))}')
psi=$(echo $eta $delta | awk '{printf("%.11f", ($1 - 2 * $2))}')
psi_1=$(echo $psi | awk '{printf("%.11f", (1 - $1))}')
psi_2=$(echo $psi | awk '{printf("%.11f", ($1 - 1))}')
phi=$(echo $zeta $mu | awk '{printf("%.11f", (1 + $1 - 2 * $2))}')
phi_1=$(echo $phi | awk '{printf("%.11f", (1 - $1))}')
phi_2=$(echo $phi | awk '{printf("%.11f", ($1 - 1))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
$mu $mu $delta ! Y

$mu $mu $delta ! Y
$phi_1 $phi_1 $psi_1 ! F

$phi_1 $phi_1 $psi_1 ! F
$zeta $zeta $eta ! H

$zeta $zeta $eta ! H
0.0 0.0 0.5 ! Z

0.0 0.0 0.5 ! Z
0.5 -0.5 0.5 ! I

$zeta_1 -$zeta $eta_1 ! H_1
$mu_1 -$mu -$delta ! Y_1

$mu_1 -$mu -$delta ! Y_1
0.5 -0.5 0.0 ! X

0.5 -0.5 0.0 ! X
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.0 ! N

0.5 0.0 0.5 ! M
0.0 0.0 0.0 ! \Gamma

EOF
;;
  MCLC5)

zeta=$(echo $a $b $c $alpha | awk '{printf("%.11f", (($2 ** 2 / $1 ** 2 + (1 - $2 * cos($4) / $3) / sin($4) ** 2) / 4))}')
zeta_1=$(echo $zeta | awk '{printf("%.11f", (1 - $1))}')
eta=$(echo $b $c $alpha $zeta | awk '{printf("%.11f", (0.5 + 2 * $4 * $2 * cos($3) / $1))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1))}')
mu=$(echo $a $b $c $alpha $eta | awk '{printf("%.11f", ($5 / 2 + $2 ** 2 / (4 * $1 ** 2) - $2 * $3 * cos($4) / (2 * $1 ** 2)))}')
mu_1=$(echo $mu | awk '{printf("%.11f", (1 - $1))}')
nu=$(echo $mu $zeta | awk '{printf("%.11f", ( 2 * $1 - $2))}')
nu_1=$(echo $nu | awk '{printf("%.11f", (1 - $1))}')
rho=$(echo $a $b $zeta | awk '{printf("%.11f", (1 - $3 * $1 ** 2 / $2 ** 2))}')
rho_1=$(echo $rho | awk '{printf("%.11f", (1 - $1))}')
rho_2=$(echo $rho | awk '{printf("%.11f", ($1 - 1))}')
omega=$(echo $a $b $c $alpha $nu | awk '{printf("%.11f", ((4 * $5 - 1 - $2 ** 2 * sin($4) ** 2 / $1 ** 2) * $3 / (2 * $2 * cos($4))))}')
omega_1=$(echo $omega | awk '{printf("%.11f", (1 - $1))}')
delta=$(echo $b $c $alpha $zeta $omega | awk '{printf("%.11f", ( $4 * $2 * cos($3) / $1 + $5 / 2 - 0.25 ))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
$mu $mu $delta ! Y

$mu $mu $delta ! Y
$nu $nu $omega ! F

$nu $nu $omega ! F
0.5 0.5 0.5 ! L

0.5 0.5 0.5 ! L
$rho $rho_1 0.5 ! I

$rho_1 $rho_2 0.5 ! I_1
0.0 0.0 0.5 ! Z

0.0 0.0 0.5 ! Z
$zeta $zeta $eta ! H

$zeta $zeta $eta ! H
$nu_1 $nu_1 $omega_1 ! F_1

$zeta_1 $zeta $eta_1 ! H_1
$mu_1 -$mu -$delta ! Y_1

$mu_1 -$mu -$delta ! Y_1
0.5 -0.5 0.0 ! X

0.5 -0.5 0.0 ! X
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.0 ! N

0.5 0.0 0.5 ! M
0.0 0.0 0.0 ! \Gamma

EOF
;;
  ORCC)

zeta=$(echo $a $b | awk '{printf("%.11f", ( (1 + $1 ** 2 / $2 ** 2) / 4 ))}')
zeta_1=$(echo $zeta | awk '{printf("%.11f", (1 - $1 ))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
$zeta $zeta 0.0 ! X

$zeta $zeta 0.0 ! X
0.0 0.5 0.0 ! S

0.0 0.5 0.0 ! S
0.0 0.5 0.5 ! R

0.0 0.5 0.5 ! R
$zeta $zeta 0.5 ! A

$zeta $zeta 0.5 ! A
0.0 0.0 0.5 ! Z

0.0 0.0 0.5 ! Z
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
-0.5 0.5 0.0 ! Y

-0.5 0.5 0.0 ! Y
-$zeta $zeta_1 0.0 ! X_1

-$zeta $zeta_1 0.0 ! X_1
-$zeta $zeta_1 0.5 ! A_1

-$zeta $zeta_1 0.5 ! A_1
-0.5 0.5 0.5 ! T

-0.5 0.5 0.5 ! T
-0.5 0.5 0.0 ! Y

0.0 0.0 0.5 ! Z
-0.5 0.5 0.5 ! T

EOF
;;
  ORCF1)

zeta=$(echo $a $b $c | awk '{printf("%.11f", ( (1 + $1 ** 2 / $2 ** 2 - $1 ** 2 / $3 ** 2) / 4 ))}')
zeta_1=$(echo $zeta | awk '{printf("%.11f", (1 - $1 ))}')
zeta_5=$(echo $zeta | awk '{printf("%.11f", (0.5 - $1 ))}')
zeta_52=$(echo $zeta | awk '{printf("%.11f", (0.5 + $1 ))}')
eta=$(echo $a $b $c | awk '{printf("%.11f", ( (1 + $1 ** 2 / $2 ** 2 + $1 ** 2 / $3 ** 2) / 4  ))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1 ))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.5 ! Y

0.5 0.0 0.5 ! Y
1.0 0.5 0.5 ! T

1.0 0.5 0.5 ! T
0.5 0.5 0.0 ! Z

0.5 0.5 0.0 ! Z
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.0 $eta $eta ! X

0.0 $eta $eta ! X
0.5 $zeta_5 $zeta_1 ! A_1

0.5 $zeta_5 $zeta_1 ! A_1
0.5 0.0 0.5 ! Y

1.0 0.5 0.5 ! T
1.0 $eta_1 $eta_1 ! X_1

0.0 $eta $eta ! X
0.5 $zeta_52 $zeta ! A

0.5 $zeta_52 $zeta ! A
0.5 0.5 0.0 ! Z

0.5 0.5 0.5 ! L
0.0 0.0 0.0 ! \Gamma

EOF
;;
  ORCF2)

phi=$(echo $a $b $c | awk '{printf("%.11f", ( (1 + $3 ** 2 / $2 ** 2 - $3 ** 2 / $1 ** 2) / 4 ))}')
phi_1=$(echo $phi | awk '{printf("%.11f", (1 - $1 ))}')
phi_5=$(echo $phi | awk '{printf("%.11f", (0.5 - $1 ))}')
phi_52=$(echo $phi | awk '{printf("%.11f", (0.5 + $1 ))}')
eta=$(echo $a $b $c | awk '{printf("%.11f", ( (1 + $1 ** 2 / $2 ** 2 - $1 ** 2 / $3 ** 2) / 4  ))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1 ))}')
eta_5=$(echo $eta | awk '{printf("%.11f", (0.5 - $1 ))}')
eta_52=$(echo $eta | awk '{printf("%.11f", (0.5 + $1 ))}')
delta=$(echo $a $b $c | awk '{printf("%.11f", ( (1 + $2 ** 2 / $1 ** 2 - $2 ** 2 / $3 ** 2) / 4  ))}')
delta_1=$(echo $delta | awk '{printf("%.11f", (1 - $1 ))}')
delta_5=$(echo $delta | awk '{printf("%.11f", (0.5 - $1 ))}')
delta_52=$(echo $delta | awk '{printf("%.11f", (0.5 + $1 ))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.5 ! Y

0.5 0.0 0.5 ! Y
0.5 $eta_5 $eta_1 ! C

0.5 $eta_5 $eta_1 ! C
$delta_5 0.5 $delta_1 ! D

$delta_5 0.5 $delta_1 ! D
0.0 0.5 0.5 ! X

0.0 0.5 0.5 ! X
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 0.5 0.0 ! Z

0.5 0.5 0.0 ! Z
$delta_52 0.5 $delta ! D_1

$delta_52 0.5 $delta ! D_1
$phi_1 $phi_5 0.5 ! H

$phi_1 $phi_5 0.5 ! H
0.5 $eta_5 $eta_1 ! C

0.5 $eta_52 $eta ! C_1
0.5 0.5 0.0 ! Z

0.0 0.5 0.5 ! X
$phi $phi_52 0.5 ! H_1

$phi_1 $phi_5 0.5 ! H
0.5 0.0 0.5 ! Y

0.5 0.5 0.5 ! L
0.0 0.0 0.0 ! \Gamma

EOF
;;
  ORCF3)

zeta=$(echo $a $b $c | awk '{printf("%.11f", ( (1 + $1 ** 2 / $2 ** 2 - $1 ** 2 / $3 ** 2) / 4 ))}')
zeta_1=$(echo $zeta | awk '{printf("%.11f", (1 - $1 ))}')
zeta_5=$(echo $zeta | awk '{printf("%.11f", (0.5 - $1 ))}')
zeta_52=$(echo $zeta | awk '{printf("%.11f", (0.5 + $1 ))}')
eta=$(echo $a $b $c | awk '{printf("%.11f", ( (1 + $1 ** 2 / $2 ** 2 + $1 ** 2 / $3 ** 2) / 4  ))}')
eta_1=$(echo $eta | awk '{printf("%.11f", (1 - $1 ))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.5 ! Y

0.5 0.0 0.5 ! Y
1.0 0.5 0.5 ! T

1.0 0.5 0.5 ! T
0.5 0.5 0.0 ! Z

0.5 0.5 0.0 ! Z
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.0 $eta $eta ! X

0.0 $eta $eta ! X
0.5 $zeta_5 $zeta_1 ! A_1

0.5 $zeta_5 $zeta_1 ! A_1
0.5 0.0 0.5 ! Y

0.0 $eta $eta ! X
0.5 $zeta_52 $zeta ! A

0.5 $zeta_52 $zeta ! A
0.5 0.5 0.0 ! Z

0.5 0.5 0.5 ! L
0.0 0.0 0.0 ! \Gamma

EOF
;;
  ORCI)

zeta=$(echo $a $c | awk '{printf("%.11f", ( (1 + $1 ** 2 / $2 ** 2) / 4 ))}')
zeta_1=$(echo $zeta | awk '{printf("%.11f", (1 - $1 ))}')
eta=$(echo $b $c | awk '{printf("%.11f", ( (1 + $1 ** 2 / $2 ** 2) / 4 ))}')
eta_1=$(echo $eta | awk '{printf("%.11f", ( 1 - $1 ))}')
delta=$(echo $a $b $c | awk '{printf("%.11f", ( ($2 ** 2 - $1 ** 2) / (4 * $3 ** 2) ))}')
delta_5=$(echo $delta | awk '{printf("%.11f", (0.5 - $1 ))}')
delta_52=$(echo $delta | awk '{printf("%.11f", (0.5 + $1 ))}')
mu=$(echo $a $b $c | awk '{printf("%.11f", ( ($1 ** 2 + $2 ** 2) / (4 * $3 ** 2) ))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
-$zeta $zeta $zeta ! X

-$zeta $zeta $zeta ! X
-$mu $mu $delta_5 ! L

-$mu $mu $delta_5 ! L
0.0 0.0 0.5 ! T

0.0 0.0 0.5 ! T
0.25 0.25 0.25 ! W

0.25 0.25 0.25 ! W
0.0 0.5 0.0 ! R

0.0 0.5 0.0 ! R
$zeta $zeta_1 -$zeta ! X_1

$zeta $zeta_1 -$zeta ! X_1
0.5 0.5 -0.5 ! Z

0.5 0.5 -0.5 ! Z
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
$eta -$eta $eta ! Y

$eta -$eta $eta ! Y
0.5 0.0 0.0 ! S

0.5 0.0 0.0 ! S
0.25 0.25 0.25 ! W

$mu -$mu $delta_52 ! L_1
$eta -$eta $eta ! Y

$eta_1 $eta -$eta ! Y_1
0.5 0.5 -0.5 ! Z

EOF
;;
  RHL1)

eta=$(echo $alpha | awk '{printf("%.11f", ( (1 + 4 * cos($1)) / (2 + 4 * cos($1)) ))}')
eta_1=$(echo $eta | awk '{printf("%.11f", ( 1 - $1 ))}')
eta_2=$(echo $eta | awk '{printf("%.11f", ( $1 - 1 ))}')
nu=$(echo $eta | awk '{printf("%.11f", ( 3.0 / 4.0 - $1 / 2.0 ))}')
nu_1=$(echo $nu | awk '{printf("%.11f", ( 1 - $1 ))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.0 ! L

0.5 0.0 0.0 ! L
0.5 $eta_1 $eta_2 ! B_1

$eta 0.5 $eta_1 ! B
0.5 0.5 0.5 ! Z

0.5 0.5 0.5 ! Z
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
$nu 0.0 -$nu ! X

$nu_1 $nu 0.0 ! Q
0.5 0.5 0.0 ! F

0.5 0.5 0.0 ! F
$nu_1 $nu_1 $eta_1 ! P_1

$nu_1 $nu_1 $eta_1 ! P_1
0.5 0.5 0.5 ! Z

0.5 0.0 0.0 ! L
$eta $nu $nu ! P

EOF
;;
  RHL2)

eta=$(echo $alpha | awk '{printf("%.11f", ( 1 / (2 * (sin($1 / 2.0) / cos($1 / 2.0) ) ** 2) ))}')
eta_1=$(echo $eta | awk '{printf("%.11f", ( 1 - $1 ))}')
nu=$(echo $eta | awk '{printf("%.11f", ( 3.0 / 4.0 - $1 / 2.0 ))}')
nu_1=$(echo $nu | awk '{printf("%.11f", ( 1 - $1 ))}')
nu_2=$(echo $nu | awk '{printf("%.11f", ( $1 - 1 ))}')

cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
$nu_1 -$nu $nu_1 ! P

$nu_1 -$nu $nu_1 ! P
0.5 -0.5 0.5 ! Z

0.5 -0.5 0.5 ! Z
$eta $eta $eta ! Q

$eta $eta $eta ! Q
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 -0.5 0.0 ! F

0.5 -0.5 0.0 ! F
$nu $nu_2 $nu_2 ! P_1

$nu $nu_2 $nu_2 ! P_1
$eta_1 -$eta -$eta ! Q_1

$eta_1 -$eta -$eta ! Q_1
0.5 0.0 0.0 ! L

0.5 0.0 0.0 ! L
0.5 -0.5 0.5 ! Z

EOF
;;
############# fine part
  BCC) 
cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.5 -0.5 0.5 ! H

0.5 -0.5 0.5 ! H
0.0 0.0 0.5 ! N

0.0 0.0 0.5 ! N
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.25 0.25 0.25 ! P

0.25 0.25 0.25 ! P
0.5 -0.5 0.5 ! H

0.25 0.25 0.25 ! P
0.0 0.0 0.5 ! N

EOF
;;
  CUB) 
cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.0 0.5 0.0 ! X

0.0 0.5 0.0 ! X
0.5 0.5 0.0 ! M

0.5 0.5 0.0 ! M
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 0.5 0.5 ! R

0.5 0.5 0.5 ! R
0.0 0.5 0.0 ! X

0.5 0.5 0.0 ! M
0.5 0.5 0.5 ! R

EOF
;;
  FCC) 
cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.5 ! X

0.5 0.0 0.5 ! X
0.5 0.25 0.75 ! W

0.5 0.25 0.75 ! W
0.375 0.375 0.75 ! K

0.375 0.375 0.75 ! K
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 0.5 0.5 ! L

0.5 0.5 0.5 ! L
0.625 0.25 0.625 ! U

0.625 0.25 0.625 ! U
0.5 0.25 0.75 ! W

0.5 0.25 0.75 ! W
0.5 0.5 0.5 ! L

0.5 0.5 0.5 ! L
0.375 0.375 0.75 ! K

0.625 0.25 0.625 ! U
0.5 0.0 0.5 ! X

EOF
;;
  HEX) 
cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.0 ! M

0.5 0.0 0.0 ! M
0.333333333333 0.333333333333 0.0 ! K

0.333333333333 0.333333333333 0.0 ! K
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.0 0.0 0.5 ! A

0.0 0.0 0.5 ! A
0.5 0.0 0.5 ! L

0.5 0.0 0.5 ! L
0.333333333333 0.333333333333 0.5 ! H

0.333333333333 0.333333333333 0.5 ! H
0.0 0.0 0.5 ! A

0.5 0.0 0.5 ! L
0.5 0.0 0.0 ! M

0.333333333333 0.333333333333 0.0 ! K
0.333333333333 0.333333333333 0.5 ! H

EOF
;;
  ORC) 
cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.0 ! X

0.5 0.0 0.0 ! X
0.5 0.5 0.0 ! S

0.5 0.5 0.0 ! S
0.0 0.5 0.0 ! Y

0.0 0.5 0.0 ! Y
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.0 0.0 0.5 ! Z

0.0 0.0 0.5 ! Z
0.5 0.0 0.5 ! U

0.5 0.0 0.5 ! U
0.5 0.5 0.5 ! R

0.5 0.5 0.5 ! R
0.0 0.5 0.5 ! T

0.0 0.5 0.5 ! T
0.0 0.0 0.5 ! Z

0.0 0.5 0.0 ! Y
0.0 0.5 0.5 ! T

0.5 0.0 0.5 ! U
0.5 0.0 0.0 ! X

0.5 0.5 0.0 ! S
0.5 0.5 0.5 ! R

EOF
;;
  TET) 
cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 0.0 0.0 ! \Gamma
0.0 0.5 0.0 ! X

0.0 0.5 0.0 ! X
0.5 0.5 0.0 ! M

0.5 0.5 0.0 ! M
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.0 0.0 0.5 ! Z

0.0 0.0 0.5 ! Z
0.0 0.5 0.5 ! R

0.0 0.5 0.5 ! R
0.5 0.5 0.5 ! A

0.5 0.5 0.5 ! A
0.0 0.0 0.5 ! Z

0.0 0.5 0.0 ! X
0.0 0.5 0.5 ! R

0.5 0.5 0.0 ! M
0.5 0.5 0.5 ! A

EOF
;;
  TRIA)
cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.5 0.0 0.0 ! X
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.0 0.5 0.0 ! Y

0.5 0.5 0.0 ! L
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.0 0.0 0.5 ! Z

0.5 0.0 0.5 ! N
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.0 0.5 0.5 ! M

0.5 0.5 0.5 ! R
0.0 0.0 0.0 ! \Gamma

EOF
;;
  TRIB)
cat >KPOINTS.band <<EOF
k-points along high symmetry lines
$kmesh
Line-mode
rec
0.0 -0.5 0.0 ! X
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.0 ! Y

0.5 -0.5 0.0 ! L
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
-0.5 0.0 0.5 ! Z

-0.5 -0.5 0.5 ! N
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.0 0.0 0.5 ! M

0.0 -0.5 0.5 ! R
0.0 0.0 0.0 ! \Gamma

EOF
;;
 *)

echo '!!!!!!!!!!!!!!!!!!!!!!!'
echo '!!!!!! NOT FOUND !!!!!!'
echo 'CHECK YOUR STRUCTURE!!!'
echo '!!!!!!!!!!!!!!!!!!!!!!!'
esac
