#!/bin/bash
star=1000000
acyclic=1000000
livarh=1000000
preacc=1000000


p1=100000
s1=100000
c1=100000
r1=100000

p0=100000
s0=100000
c0=100000
r0=100000

sh -c "make"

for i in {1..10}
do
    echo $i
    sh -c "./col_hessian > temp.out"
    starold=$(grep "STAR" temp.out | cut -d' ' -f5)
    ret=$(echo "$starold<$star" | bc)
    if [[ $ret -eq 1 ]]; then
        star=$starold
    fi
#breakdown star
    p1old=$(grep "P1" temp.out | cut -d' ' -f2)
    ret=$(echo "$p1old<$p1" | bc)
    if [[ $ret -eq 1 ]]; then
        p1=$p1old
    fi
    s1old=$(grep "S1" temp.out | cut -d' ' -f2)
    ret=$(echo "$s1old<$s1" | bc)
    if [[ $ret -eq 1 ]]; then
        s1=$s1old
    fi
    c1old=$(grep "C1" temp.out | cut -d' ' -f2)
    ret=$(echo "$c1old<$c1" | bc)
    if [[ $ret -eq 1 ]]; then
        c1=$c1old
    fi
    r1old=$(grep "R1" temp.out | cut -d' ' -f2)
    ret=$(echo "$r1old<$r1" | bc)
    if [[ $ret -eq 1 ]]; then
        r1=$r1old
    fi



    acyclicold=$(grep "ACYCLIC" temp.out | cut -d' ' -f5)
    ret=$(echo "$acyclicold<$acyclic" | bc)
    if [[ $ret -eq 1 ]]; then
        acyclic=$acyclicold
    fi
#breakdown acyclic
    p0old=$(grep "P0" temp.out | cut -d' ' -f2)
    ret=$(echo "$p0old<$p0" | bc)
    if [[ $ret -eq 1 ]]; then
        p0=$p0old
    fi
    s0old=$(grep "S0" temp.out | cut -d' ' -f2)
    ret=$(echo "$s0old<$s0" | bc)
    if [[ $ret -eq 1 ]]; then
        s0=$s0old
    fi
    c0old=$(grep "C0" temp.out | cut -d' ' -f2)
    ret=$(echo "$c0old<$c0" | bc)
    if [[ $ret -eq 1 ]]; then
        c0=$c0old
    fi
    r0old=$(grep "R0" temp.out | cut -d' ' -f2)
    ret=$(echo "$r0old<$r0" | bc)
    if [[ $ret -eq 1 ]]; then
        r0=$r0old
    fi

    livarhold=$(grep "LIVARHNO" temp.out | cut -d' ' -f5)
    ret=$(echo "$livarhold<$livarh" | bc)
    if [[ $ret -eq 1 ]]; then
        livarh=$livarhold
    fi
    preaccold=$(grep "LIVARHACC" temp.out | cut -d' ' -f5)
    ret=$(echo "$preaccold<$preacc" | bc -q)
    if [[ $ret -eq 1 ]]; then
        preacc=$preaccold
    fi
done
    
echo "direct: "$star
echo "P: " $p1
echo "S: " $s1
echo "C: " $c1
echo "R: " $r1
dcheck=$(echo "$p1+$s1+$c1+$r1" | bc)
echo "dcheck: "$dcheck
echo

echo "indirect: "$acyclic
echo "P: " $p0
echo "S: " $s0
echo "C: " $c0
echo "R: " $r0
icheck=$(echo "$p0+$s0+$c0+$r0" | bc)
echo "icheck: "$icheck
echo

echo "livarh: "$livarh
echo "preacc: "$preacc

sh -c "rm temp.out"
