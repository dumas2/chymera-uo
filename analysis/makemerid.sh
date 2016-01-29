#!/bin/tcsh
# $1 start file number
# $2 end file number
# $3 advance file
# $4 time file (output of sigplot)
# $5 plot title
# $6 input file prefix contain dierctory
# $7 output file prefix contain directory

  if ( $#argv != 7 ) then
       printf "\n\n\n A CHYMTOOL SPONSORED EVENT\n\n" ;
       echo makeplot v 1.  11.10.2007. Nora Bolig
       printf "\n\n";
       echo THERE SHOULD BE SEVEN ARGUMENTS
       printf "\n";
       echo Usage: $0 '{start file number} {end file number} {file number skip}'
       echo '{sigplot output file} {"PLOT TITLE"} {dir/input_prefix} {dir/output_prefix} '
       printf "\n";
       echo "Calling Cthulu...(run)"
       exit 127;
  endif

 printf "\n\n\n A CHYMTOOL SPONSORED EVENT\n\n" ;
 printf " makeplot v 1.  11.10.2007. Nora Bolig\n\n\n";
 printf " You typed: $0 $1 $2 $3 $4 '$5' $6 $7\n";
 printf "\n Check image for correct contours.\n\n";
 

 set dum1 = ` gnuplot --version | awk '{print substr($2,1,1)}' `
 set dum2 = ` gnuplot --version | awk '{print substr($2,3,3)}' `

 if ( $dum1 < 4 )  then
   echo ' Please update your version of gnuplot to 4.2 or higher '
   echo ' You are using' `gnuplot --version`
   exit
 else if ( $dum2 < 2 ) then
   echo ' Please update your version of gnuplot to 4.2 or higher '
   echo ' You are using' `gnuplot --version`
   exit
 endif

 @ i = $1 # starting file number

 while ( $i < $2 )  # endfile file number

 @ j = $i
 set j = `printf "%06d" $j`
   
 set l = `grep $j $4 | grep FILE | awk '{print $4}'`

gnuplot << LEAVE

unset key
set pm3d map
set size ratio 1
set view 0,0
set title "$5"
set tics out
set xlabel 'r(AU)'
set ylabel 'z(AU)'
set cblabel 'Log {/Symbol S}(g cm^{-2})'
set palette color
set colorbox
#set cbrange [-1:3]
set terminal postscript eps color enhanced  size 8in,8in font "Helvetica" 36
set output "$7$j.eps" 
plot  '$6$j' with image

LEAVE
 
 echo file $i
 @ i = $i + $3

 end

exit 0;
   
