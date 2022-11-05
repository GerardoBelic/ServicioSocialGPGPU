cat nativacg.pdb | awk '{printf "%s  %5i  %-3s %3s %1s %3i    %8.3f%8.3f%8.3f\n",$1,$2,$3,$4,"B",$6,$7+10.,$8,$9}' > copia.pdb
