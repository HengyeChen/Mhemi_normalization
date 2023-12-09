#!bin/bash
prefix="GM12878_Mhemi_rep123_CpG_dyad" #prefix of file name
E_CGNA=0.6 #cutting efficiency at CGNA
E_CGNG=0.95 #cutting efficiency at CGNG
E_PU=4 #purification efficiency at CGNG
#detemine the YNCGNR motif at each CpG
bedtools intersect -a $prefix.bed -b hg38_YNCGNR.fa.bed -wo > $prefix.YNCGNR.bed

#normalize the count of unme, hemiW, hemiC, and me
cutoff1=1
fd1=3
motif="CNCGNG"
Ew=$E_CGNG
Ec=$E_CGNG
grep -i 'C[A,T,C,G]CG[A,T,C,G]G' $prefix.YNCGNR.bed | awk -v Ew=$Ew -v Ec=$Ec -v Ep=$E_PU '{
    if(($4+$5+$6+$7)>2 && $2==($9+2))
    {
        unme=$4;hw=$5;hc=$6;me=$7
        if($7!=0)
        {
            me=($7/Ew/Ec*Ep);unme=0;hw=0;hc=0;
            if($4>me*(1-Ew)*(1-Ec)) unme=($4-me*(1-Ew)*(1-Ec));
            if($5>me*Ew*(1-Ec)) hw=($5-me*Ew*(1-Ec));
            if($6>me*Ec*(1-Ew)) hc=($6-me*Ec*(1-Ew))
        }
        else
        {
            me=0;unme=$4;hw=$5;hc=$6
        }
        print $1,$2,$3,unme,hw,hc,me,$8,$9,$10,$11
    }
    else if($2==($9+2))
    {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
    }         
}' OFS='\t' > $prefix.${motif}.v4.3.bed


motif='TNCGNG'
Ew=$E_CGNG
Ec=$E_CGNA
grep -i 'T[A,T,C,G]CG[A,T,C,G]G' $prefix.YNCGNR.bed | awk -v Ew=$Ew -v Ec=$Ec -v cutoff1=$cutoff1 -v fd1=$fd1 -v Ep=$E_PU '{
    if(($4+$5+$6+$7)>2 && $2==($9+2))
    {
        unme=$4;hw=$5;hc=$6;me=$7
        cutoff2=1*Ep/Ew/Ec
        if($7!=0)
        {
            unme=0;hw=0;hc=0;
            me=($7/Ew/Ec*Ep);
            if($4>unme*(1-Ew)*(1-Ec)) unme=($4-me*(1-Ew)*(1-Ec)); if($5>me*Ew*(1-Ec)) hw=($5-me*Ew*(1-Ec)); if($6>me*Ec*(1-Ew)) hc=($6-me*Ec*(1-Ew))
        }
        else if($7==0 && $5!=0 && $5>=$6)
        {
            unme=$4;hw=$5;hc=$6;
            if($5>cutoff2*Ew*(1-Ec))
            {
                me=cutoff2
                hc=$6-cutoff2*Ec*(1-Ew)
                hw=$5-cutoff2*Ew*(1-Ec)
                if(hc<0)
                {
                    hc=0
                }
                unme=$4-me*(1-Ew)*(1-Ec)
                if(unme<0)
                {
                    unme=0
                }
            }
            else if($5<=cutoff2*Ew*(1-Ec) && $5>=cutoff1)
            {
                me=($5-cutoff1)/Ew/(1-Ec)*Ep
                hc=$6-cutoff2*Ec*(1-Ew)
                unme=$4-me*(1-Ew)*(1-Ec)
                hw=$5-cutoff2*Ew*(1-Ec)
                if(hw<0)
                {
                    hw=0
                }
                if(unme<0)
                {
                    unme=0
                }
            }
        }
        else if($7==0 && $5!=0 && $6>$5)
        {
            me=$5/Ew/(1-Ec)
            hc=$6-me*Ec*(1-Ew)
            hw=0
            if(hc<0)
            {
                hc=0
            }
        }
        if(unme<0) unme=0;if(hw<0) hw=0;if(hc<0) hc=0;if(me<0) me=0
        print $1,$2,$3,unme,hw,hc,me,$8,$9,$10,$11
    }
    else if($2==($9+2))
    {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
    }     
}' OFS='\t' > $prefix.${motif}.v4.3.bed

motif='CNCGNA'
Ew=$E_CGNA
Ec=$E_CGNG
grep -i 'C[A,T,C,G]CG[A,T,C,G]A' $prefix.YNCGNR.bed | awk -v Ew=$Ew -v Ec=$Ec -v cutoff1=$cutoff1 -v fd1=$fd1 -v Ep=$E_PU '{
    if(($4+$5+$6+$7)>2 && $2==($9+2))
    {
        unme=$4;hw=$5;hc=$6;me=$7
        cutoff2=1*Ep/Ew/Ec
        if($7!=0)
        {
            unme=0;hw=0;hc=0;
            me=($7/Ew/Ec*Ep);
            if($4>unme*(1-Ew)*(1-Ec)) unme=($4-me*(1-Ew)*(1-Ec)); if($5>me*Ew*(1-Ec)) hw=($5-me*Ew*(1-Ec)); if($6>me*Ec*(1-Ew)) hc=($6-me*Ec*(1-Ew))
        }
        else if($7==0 && $6!=0 && $6>=$5)
        {
            unme=$4;hw=$5;hc=$6;
            if($6>cutoff2*Ec*(1-Ew))
            {
                me=cutoff2
                hc=$6-cutoff2*Ec*(1-Ew)
                hw=$5-cutoff2*Ew*(1-Ec)
                if(hw<0)
                {
                    hw=0
                }
                unme=$4-me*(1-Ew)*(1-Ec)
                if(unme<0)
                {
                    unme=0
                }
            }
            else if($6<=cutoff2*Ec*(1-Ew) && $6>=cutoff1)
            {
                me=($6-cutoff1)/Ec/(1-Ew)*Ep
                hc=cutoff1
                unme=$4-me*(1-Ew)*(1-Ec)
                hw=$5-cutoff2*Ew*(1-Ec)
                if(hw<0)
                {
                    hw=0
                }
                if(unme<0)
                {
                    unme=0
                }
            }
        }
        else if($7==0 && $5!=0 && $6<$5)
        {
            me=$6/Ew/(1-Ec)
            hc=0
            hw=$5-me*Ec*(1-Ew)
            if(hw<0)
            {
                hw=0
            }
        }
        if(unme<0) unme=0;if(hw<0) hw=0;if(hc<0) hc=0;if(me<0) me=0
        print $1,$2,$3,unme,hw,hc,me,$8,$9,$10,$11
    }
    else if($2==($9+2))
    {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
    }    
}' OFS='\t' > $prefix.${motif}.v4.3.bed


motif='TNCGNA'
Ew=$E_CGNA
Ec=$E_CGNA
grep -i 'T[A,T,C,G]CG[A,T,C,G]A' $prefix.YNCGNR.bed | awk -v Ew=$Ew -v Ec=$Ec -v cutoff1=$cutoff1 -v fd1=$fd1 -v Ep=$E_PU '{
    if(($4+$5+$6+$7)>2 && $2==($9+2))
    {
        unme=$4;hw=$5;hc=$6;me=$7
        if($7!=0)
        {
            unme=0;hw=0;hc=0;
            me=($7/Ew/Ec*Ep);
            if($4>unme*(1-Ew)*(1-Ec)) unme=($4-me*(1-Ew)*(1-Ec)); if($5>me*Ew*(1-Ec)) hw=($5-me*Ew*(1-Ec)); if($6>me*Ec*(1-Ew)) hc=($6-me*Ec*(1-Ew))
        }
        else if($5!=0 || $6!=0)
        {
            cutoff2=1/Ew/Ec*Ep
            unme=$4;hw=$5;hc=$6;
            if(unme>3*1/(1-Ew)/Ec*(1-Ew)*(1-Ec))
            {
                unmeadd=1/(1-Ew)/Ec
            }
            else
            {
                unmeadd=0
            }
            if($5>=$6)
            {
                if($6<=cutoff2*Ec*(1-Ew))
                {
                    hw=$5-$6
                    hc=0
                    me=1/Ec/(1-Ew)
                }
                else
                {
                    hw=$5-cutoff2*Ec*(1-Ew)
                    hc=$6-cutoff2*Ec*(1-Ew)
                    me=cutoff2*Ec*(1-Ew)
                    unme=$4-cutoff2*(1-Ec)*(1-Ew)
                }
            }
            else if($5<$6)
            {
                if($5<=cutoff2*Ew*(1-Ec))
                {
                    hw=0
                    hc=$6-$5
                    me=$5/(1-Ew)/Ec
                }
                else
                {
                    hw=$5-cutoff2*Ew*(1-Ec)
                    hc=$6-cutoff2*Ew*(1-Ec)
                    me=cutoff2*Ew*(1-Ec)
                    unme=$4-cutoff2*(1-Ew)*(1-Ec)
                    if(unme<0) unme=0                    
                }
            }
            if(me<unmeadd)
            {
                me=unmeadd
            }
        }
        else
        {
            cutoff2=1/Ew/Ec*Ep
            unme=$4;hw=$5;hc=$6;
            if(unme>3*1/Ew/(1-Ec)*(1-Ew)*(1-Ec))
            {
                me=1/Ew/(1-Ec)
                unme=$4-1/Ew/(1-Ec)*(1-Ew)*(1-Ec)
            }
            else
            {
                me=0
            }
        }
        if(unme<0) unme=0;if(hw<0) hw=0;if(hc<0) hc=0;if(me<0) me=0
        print $1,$2,$3,unme,hw,hc,me,$8,$9,$10,$11
    }
    else if($2==($9+2))
    {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
    }
}' OFS='\t' > $prefix.${motif}.v4.3.bed

cat $prefix.*.v4.3.bed | awk '{if(($4+$5+$6+$7)>2 || ($7==($4+$5+$6+$7) && $7>0)) print $0}' | sort -k1,1 -k2,2n > $prefix.YNCGNR.v4.3.bed