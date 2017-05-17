#!usr/bin/perl  -w
use strict;
#---------------------------------------------------------------------------
#脚本说明：对annnovar注释后的csv文件进行格式调整，生成下游要求的文件格式
#作者：谢剑邦 2016-10-21
#使用方法：perl  <Sample.vcf> <Sample.csv> <n>（vcf中样本数）
#注意项：
#1.从vcf文件中提取基因型和深度信息，注释到csv文件中。
#2.注释omim编号
#3.注释hgmd的ID和variant class
#4.提取variant的核苷酸和氨基酸信息
#5.去除繁杂的注释信息
#6.插入自动生报告的9列信息
#7.通过list过滤
#--------------------------------------------------------------------------
my $usage="perl $0 <Sample.vcf> <Sample.csv> <n> \n";
die $usage if @ARGV!=3;
#----------------------------------------
###从vcf文件中提取基因型和深度信息
#----------------------------------------

#提取vcf中深度和基因型等信息
sub extract{
	my $data =shift @_;
	my @gtarray =split /:/,$data;
	my $l = scalar@gtarray;
        my $value;
        if ($l!=5){####可能存在DP=0
                $value = "$gtarray[0],.,.";
        }else{
                pop @gtarray ;
                pop @gtarray ;
                my @AD = split /,/,$gtarray[1];
		$gtarray[1] = join ";",@AD;
                $value = join ",",@gtarray;
        }
	return $value;
}

my $header_info;
print "Loading Databases...\n";
my %genetype;
open VCF,"$ARGV[0]";
while(<VCF>){
        chomp;
        s/\r//g;
	next if /^##/;
	my $N = $ARGV[2];
	my @title_set;
	my @array = split;
        if($N==1){
                if(/^#CHROM/){
                        my @title = split ;
                        my $Allele_GT = "Geno_Type($title[-1])";
                        my $Allele_AD = "Allelic_depths($title[-1])";
                        my $Allele_DP = "Seq_depths($title[-1])";
                        my $title_info = join (",",$Allele_GT,$Allele_AD,$Allele_DP);
                        unshift @title_set,$title_info;
                        $header_info = join ",",@title_set;next;
                }
                my @value;
                if($array[-1]=~/^\.\/\./){ ##信息不全
                        my $key = $array[0]."\t".$array[1];
                        $genetype{$key}="./.,.,.";
                }
		else{
                        my $Vinfo = &extract($array[-1]);
                        unshift @value,$Vinfo;
                        my $value = join ",",@value;
                        my $key = $array[0]."\t".$array[1];
                        $genetype{$key}=$value;
                }

        }else{
		if(/^#CHROM/){
                my @title = split ;
                        for(1..$N){
                                my $Allele_GT = "Geno_Type($title[-$_])";
                                my $Allele_AD = "Allelic_depths($title[-$_])";
                                my $Allele_DP = "Seq_depths($title[-$_])";
                                my $title_info = join (",",$Allele_GT,$Allele_AD,$Allele_DP);
                                unshift @title_set,$title_info;
                        }
                        $header_info = join ",",@title_set;next;
                }
                my @valueN;
                for(1..$N){
                        if($array[-$_]=~/^\.\/\./){ ##信息不全
                                my $Vinfo = "./.,.,.";
                                unshift @valueN,$Vinfo;
                        }else{
                                my $Vinfo = &extract($array[-$_]);
                                unshift @valueN,$Vinfo;
                        }
                }
                my $value = join ",",@valueN;
                my $key = $array[0]."\t".$array[1];
                $genetype{$key}=$value;
        }

}
close VCF;

#---------------------------------------------
###提取omim的ID
#---------------------------------------------
my %omim;
open OMIM ,"/data/work/xiejb/database/mim2gene.txt";
while(<OMIM>){
	chomp;
	next if /^#/;
	my @array = split;
	next if !$array[3];
	$omim{$array[3]}=$array[0];
}
close OMIM;

#----------------------------------------------------------------------
#提取hgmd的数据信息，建立hash数据库
#1.snp 以绝对位置、转录本号、ref和alt作为key；hgmd的ID和variant class做values
#2.del 以绝对位置、转录本号和del+ref作为key；hgmd的ID和variant class做values
#3.indel 以绝对位置范围作为key、转录本号和indel信息；hgmd的ID和variant class做values
#----------------------------------------------------------------------
##读SNP
my %snp;
open SNP,"/data/work/xiejb/database/hgmd/hgmd-snp.txt";
while(<SNP>){
	chomp;
	my @array = split;
	$array[3] = "." if !$array[3];
	$array[5] = "." if !$array[5];
	my $key;
	if($array[6] = ~/(NM_\d+)\.\d*:c\.\d*(\w)-(\w)/){
		my $mut_info = $1.":".$2.":".$3;
		$key = join (":",$array[0],$array[1],$mut_info);
		$snp{$key} = $array[3].",".$array[5];
	}
}
close SNP;

my %del;
open DEL,"/data/work/xiejb/database/hgmd/hgmd-indel-dian.txt";
while(<DEL>){
	chomp;
	my @array = split;
	$array[3] = "." if !$array[3];
    $array[5] = "." if !$array[5];
	my $key;
	if($array[6] = ~/(NM_\d+)\.\d*:c\.\d*(del\w)/){
		$key = join (":",$array[0],$array[1],$1,$2);
		$del{$key} =  $array[3].",".$array[5];	
	}
}
close DEL;


my %indel;
open INDEL,"/data/work/xiejb/database/hgmd/hgmd-indel.txt";
while(<INDEL>){
	chomp;
	my $key;
	my @array = split;
	$array[3] = "." if !$array[3];
	$array[5] = "." if !$array[5] ;
	if($array[6] = ~/(NM_\d+)\.\d*:c\.\d*(\w+)/){
		$key = join (":",$array[0],$array[1],$1,$2);
		$indel{$key} =  $array[3].",".$array[5];
		
	}
		
}
close INDEL;

#------------------------------------------------
##测试数据
#foreach my $key (keys%indel){
#	my $value = $indel{$key};
#	print "$key\t$value\n";
#}
#--------------------------------------------------


#-------------------------------------------------
#方便注释hgmd，使用两个子函数进行注释。
###hgmd注释，snp和del输入键找出hgmd信息，没有就输出"."。
#1.hgmdspot对snp和单个核苷酸缺失进行注释
#2.hgmdindel对indel进行注释
#-------------------------------------------------
sub hgmdspot{
	my $tmp =shift @_;
	my $spot_info;
	if(exists $snp{$tmp}){
		$spot_info = $snp{$tmp};
	}elsif(exists $del{$tmp}){
		$spot_info = $del{$tmp};
	}else{
		$spot_info = ".,.";
	}
	return $spot_info;
}

####hgmd注释，只对indel。
sub hgmdindel{
	my $tmp =shift @_;
	my $indel_info;
	if(exists $indel{$tmp}){
		$indel_info = $indel{$tmp};
	}else{
		my @sample_key = split /:/,$tmp;
		my @pos = split /-/,$sample_key[1];
		foreach my $key (keys%indel){
			my @key_array = split /:/,$key;
			if($key_array[0] eq $sample_key[0]){
				my @key_pos = split /-/,$key_array[1];
				if ((($pos[0]<=$key_pos[0]) && ($pos[1]>=$key_pos[0])) ||( ($pos[0]>=$key_pos[0]) && ($pos[0]<=$key_pos[1]) )  ){
					my @value = split /,/,$indel{$key};
					$indel_info = $value[0].",?".$value[1];
				}
			}
		}
		 $indel_info = ".,." if !$indel_info;
	}
	return $indel_info;
	
}



#----------------------------------------------------------------
###处理注释文件（csv文件）
#如果注释CSV文件结构改变，只用修改文件中对应的参数
#----------------------------------------------------------------
print "Process is running...\n";
#修改输出文件名
open IN,"$ARGV[1]";
my @name  = split /\./,$ARGV[1];
my $SampleID = $name[0];
my $postfix = pop @name;
my $format = "standard";
#my $filterZF = "filtrGL";
my $new_name = join ("\.",@name,$format,$postfix);
#my $filter_name = join ("\.",@name,$format,$filterZF,$postfix);
open OUT,"> $new_name";
#open FGL,"> $filter_name";
###打印标题
print OUT "Chr,Start,End,Ref,Alt,$header_info,MutationType,ExonicFunc.refGene,Refer_Gene,dbSNP,Nucleotide_Change,AA_Change,Gene_Annotation,Variant_Assessment,Phyological_System,Variant_Risk,System_Assesment,System_Risk,CLINSIG,Hgmd_Variant_class,ExAC_EAS,1000g2015aug_eas,ExAC_ALL,1000g2015aug_all,MCAP,REVEL,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,omimID,hgmdID,AAChange.refGene\n";
while(<IN>){
	chomp;
	s/\r//g;
	
	####csv文件去双引号和字符中的逗号
	my @semicolon_array =  split /\"/;##以双引号为分割符
	my $l = scalar @semicolon_array;
	#print "$semicolon_array[3]\n";
	###分割后的数组中奇数项为字符串项，在这里替换所有的逗号。
        for (my $i=0;$i<=$l;$i++){
                ####找出奇数项处理
				my $yu = $i%2;
                if($yu == 1 ){
                        #print "$semicolon_array[$i]\n";
			if($semicolon_array[$i] &&  $semicolon_array[$i]=~/,/){####是奇数项并不为空，含有逗号的项进行处理
			#print "$semicolon_array[$i]\n";
			#$semicolon_array[$i] = s/,/\|/g;
			my @convert = split /,/,$semicolon_array[$i];
			$semicolon_array[$i] = join "\|",@convert;
			#print "$semicolon_array[$i]\n";
			}
                }
        }

	my $data = join "",@semicolon_array;####替换逗号后，重新合并进行后面的分析
#---------------------------------------------------------
####标题名核对，如果不对给出提示
#---------------------------------------------------------
	if (/^Chr/){
                my @title = split /,/,$data;
               
                print "CLINSIG don't match\n" if $title[10] ne "CLINSIG";print "ExAC_ALL don't match\n" if $title[16] ne "ExAC_ALL";
                print "ExAC_EAS don't match\n" if $title[19] ne "ExAC_EAS";print "1000g2015aug_eas don't match\n" if $title[24] ne "1000g2015aug_eas";
                print "1000g2015aug_all don't match\n" if $title[25] ne "1000g2015aug_all";print "SIFT_score don't match\n" if $title[30] ne "SIFT_score";
                print "SIFT_pred don't match\n" if $title[31] ne "SIFT_pred";print "Polyphen2_HDIV_score don't match\n" if $title[32] ne "Polyphen2_HDIV_score";
                print "Polyphen2_HDIV_pred don't match\n" if $title[33] ne "Polyphen2_HDIV_pred";print "Polyphen2_HVAR_score don't match\n" if $title[34] ne "Polyphen2_HVAR_score";
                print "Polyphen2_HVAR_pred don't match\n" if $title[35] ne "Polyphen2_HVAR_pred";print "LRT_score don't match\n" if $title[36] ne "LRT_score";
                print "LRT_pred don't match\n" if $title[37] ne "LRT_pred";print "MutationTaster_score don't match\n" if $title[38] ne "MutationTaster_score";
                print "MutationTaster_pred don't match\n" if $title[39] ne "MutationTaster_pred";print "MutationAssessor_score don't match\n" if $title[40] ne "MutationAssessor_score";
                print "MutationAssessor_pred don't match\n" if $title[41] ne "MutationAssessor_pred";print "FATHMM_score don't match\n" if $title[42] ne "FATHMM_score";
                print "FATHMM_pred don't match\n" if $title[43] ne "FATHMM_pred";

                next;
        }

#----------------------------------------------------------
#按行对数据进行分析
#----------------------------------------------------------

	my @array = split /,/,$data;
	
#	print "$array[6]\t$array[9]\n";
	my $genename = $array[6];####基因名
	###加omim编号
	my $omimNo;
	my @glist = $array[6];	
	for(@glist){
		if(exists $omim{$_}){
			$omimNo = $omim{$_};
		}else{
			$omimNo = ".";
		}
	}


	###将要输出文件的，最后14列固化。固化前5列。
	my $backward = join (",",$array[28],$array[29],$array[30],$array[31],$array[32],$array[33],$array[34],$array[35],$array[36],$array[37],$array[38],$array[39],$array[40],$array[41],$array[42],$array[43]);
	my $forward = join (",",$array[0],$array[1],$array[2],$array[3],$array[4]);
	my $key  = $array[0]."\t".$array[1];
	
	#-----------------------------------------------------
	#注释基因型
	#-----------------------------------------------------
	my $genetp;
	if (exists $genetype{$key}){
		$genetp = $genetype{$key};
	}else{
		my $length = $array[2]-$array[1];
		if ($length != 0){
			my $newstatpos = $array[1] - 1;
			my $key1  = $array[0]."\t".$newstatpos;
			if (exists $genetype{$key1}){
               			 $genetp = $genetype{$key1};
				#print "$genetp\n";
       			 }
		}else{
			my $newstatpos = $array[1] - 1;
                        my $key1  = $array[0]."\t".$newstatpos;
                        if (exists $genetype{$key1}){
                                 $genetp = $genetype{$key1};
			}
		}
	}
	$genetp=".,.,." if !$genetp;
	#-----------------------------------------------------
	###从第九列中分割核苷酸和氨基酸的突变
	#-----------------------------------------------------
	my ($Nucleotide_Change,$AA_Change,$hgmd);
	my @mut_info = split /\|/,$array[9];
	my (@nucleotide,@amino_acid);
	my $temp_Nuc= $array[0]."-".$array[1].":".$array[3].">".$array[4];
	for(@mut_info){
		if(/exon/){
			my @divide = split /:/;
			my $amino_acid_info = pop @divide ;
			my $nucleotide_info = pop @divide ;
			my @aa_tmp = split /\./,$amino_acid_info;##提取这转录本下的氨基酸突变
			my @n_tmp = split /\./,$nucleotide_info;
			my $amino_acid = pop @aa_tmp;
			my $nucleotide = pop @n_tmp;
			
			if($nucleotide=~/^\d+_\d*(ins.+)$/ ||$nucleotide=~/^\d*(dup\w+)$/ ){
				$nucleotide=$1;
				if($nucleotide = ~/dup/){
					my $dup_end = $array[2] + 1; ##dup 在sample中结束位置没有加1.
					$hgmd = &hgmdindel($array[0].":".$array[1]."-".$dup_end.":".$divide[1].":".$nucleotide);
				}else{
					$hgmd = &hgmdindel($array[0].":".$array[1]."-".$array[2].":".$divide[1].":".$nucleotide);
				}	
			}elsif($nucleotide=~/del/){
				if(length($array[3]) == 1){
				$nucleotide=$array[3].">.";
					$hgmd = &hgmdspot($array[0].":".$array[1].":".$divide[1].":del".$array[3]);
				}else{
					$nucleotide=$array[4].">.";
					$hgmd = &hgmdindel($array[0].":".$array[1]."-".$array[2].":".$divide[1].":del".$array[4]);
				}
			}else{
				$nucleotide=~/^(\w)\d+(\w)$/;
				#if(!$1 || !$2){print $nucleotide ,"\n";}
				my $ref = $1; my $alt = $2;
				$nucleotide = $ref.">".$alt;
				$hgmd = &hgmdspot($array[0].":".$array[1].":".$divide[1].":".$ref.":".$alt);
			}

			push @nucleotide,$nucleotide;
			push @amino_acid,$amino_acid;

		}else{
			#$Nucleotide_Change = ".";
			$AA_Change =  ".";
		}
	}

	my %nucl;@nucleotide= grep { ++$nucl{$_} < 2 } @nucleotide;
	
	my $Nucleotide = join ";",@nucleotide;
	if ($Nucleotide){
		#$Nucleotide=~//;
		$Nucleotide_Change = $array[0]."-".$array[1].":".$Nucleotide;
		my %hash;@amino_acid = grep { ++$hash{$_} < 2 } @amino_acid;
        	$AA_Change =  join ";",@amino_acid;
	}else{
		#$Nucleotide_Change = ".";
		$AA_Change = ".";
	}
	
	$hgmd = ".,." if !$hgmd;
	$genetp = "." if !$genetp;
	my @ahgmd = split /,/,$hgmd;
	my $EAAF = $array[19].",".$array[24].",".$array[16].",".$array[25];
	print OUT "$forward,$genetp,$array[5],$array[8],$genename,$array[27],$temp_Nuc,$AA_Change,-,-,-,1,1,1,$array[10],$ahgmd[1],$EAAF,$backward,$omimNo,$ahgmd[0],$array[9],$SampleID\n";
}
close IN;
close OUT;
print "Process End!!!\n";
