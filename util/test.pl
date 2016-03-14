#!/usr/bin/perl -w

@array1=( "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9" );
@array2=( 6, 2, 4, 5 );
@array4=( 0.0006, 0.2, 0.004, 0.05 );
@array3=@array1[@array2];

print join(":",@array3),"\n";
@qq = 100*@array4;
@qq = map {1000000*$_} @array4;
print join(" ",@qq),"\n";
print STDOUT (sprintf("%6.6d " x ($#array4+1),@qq),"\n");
