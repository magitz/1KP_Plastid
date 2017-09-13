#/usr/bin/perl

open PAUP_FILE, ">export_genes.nex" or die "Can't open paup file: $!";
my @gene_list=(accD, atpA, atpB, atpE, atpF, atpH, atpI, ccsA, cemA, clpP, infA, matK, ndhA, ndhB, ndhC, ndhD, ndhE, ndhF, ndhG, ndhH, ndhI, ndhJ, ndhK, petA, petB, petD, petG, petL, petN, psaA, psaB, psaC, psaI, psaJ, psbA, psbB, psbC, psbD, psbE, psbF, psbH, psbI, psbJ, psbK, psbL, psbM, psbN, psbT, psbZ, rbcL, rpl14, rpl16, rpl2, rpl20, rpl22, rpl23, rpl32, rpl33, rpl36, rpoA, rpoB, rpoC1, rpoC2, rps11, rps12, rps14, rps15, rps16, rps18, rps19, rps2, rps3, rps4, rps7, rps8, ycf2, ycf3, ycf4);

for my $gene(@gene_list){

	print PAUP_FILE "exclude all;\n";
    print PAUP_FILE "include $gene;\n";
	print PAUP_FILE "export file=individual_genes/$gene.data format=text;\n";
}

