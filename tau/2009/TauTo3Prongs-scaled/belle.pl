#!/usr/local/bin/perl

# print $ARGV[0]."\n";

$rescale = 1;
$rescale = 0.930930534;
$sfac1 = 1.61983 * $rescale;
$sfac2 = 2.40957 * $rescale;
$sfac3 = 2.32292 * $rescale;
$sfac4 = 5.43807 * $rescale;
$sfac5 = 1.42995 * $rescale;

printf "* B(tau- -> pi- pi- pi+ nu) [ex K0]\n";
printf "BEGIN    Belle PimPimPipNu submitted arXiv:1001.0083 (2010)\n";
printf "\n";
printf "MEASUREMENT m_PimPimPipNu statistical systematic\n";
printf "DATA        m_PimPimPipNu statistical systematic\n";
printf "            8.42000e-02   %.5e 0.00000e+00\n", 2.58811e-03*$sfac1;
printf "\n";
printf "STAT_CORR_WITH Belle PimKmPipNu submitted 0.1749885\n";
printf "STAT_CORR_WITH Belle PimKmKpNu submitted 0.04948276\n"; 
printf "STAT_CORR_WITH Belle KmKmKpNu submitted -0.05346557\n";
printf "\n";
printf "END\n";
printf "\n";
printf "* B(tau- -> pi- K- pi+ nu) [ex. K0]\n";
printf "BEGIN    Belle PimKmPipNu submitted arXiv:1001.0083 (2010)\n";
printf "\n";
printf "MEASUREMENT m_PimKmPipNu statistical systematic\n";
printf "DATA        m_PimKmPipNu statistical systematic\n";
printf "            3.30000e-03  %.5e 0.00000e+00\n", 1.66737e-04*$sfac2;
printf "\n";
printf "STAT_CORR_WITH Belle PimPimPipNu submitted 0.1749885\n";
printf "STAT_CORR_WITH Belle PimKmKpNu submitted 0.08026913\n";
printf "STAT_CORR_WITH Belle KmKmKpNu submitted 0.03505142\n";
printf "\n";
printf "END\n";
printf "\n";
printf "* B(tau- -> pi- K- K+ nu)\n";
printf "BEGIN    Belle PimKmKpNu submitted arXiv:1001.0083 (2010)\n";
printf "\n";
printf "MEASUREMENT m_PimKmKpNu statistical systematic\n";
printf "DATA        m_PimKmKpNu statistical systematic\n";
printf "            1.55000e-03 %.5e 0.00000e+00\n", 5.59666e-05*$sfac3;
printf "\n";
printf "STAT_CORR_WITH Belle PimPimPipNu submitted 0.04948276\n";
printf "STAT_CORR_WITH Belle PimKmPipNu submitted 0.08026913\n";
printf "STAT_CORR_WITH Belle KmKmKpNu submitted -0.008312885\n";
printf "\n";
printf "END\n";
printf "\n";
printf "* B(tau- -> K- K- K+ nu)\n";
printf "BEGIN    Belle KmKmKpNu submitted arXiv:1001.0083 (2010)\n";
printf "\n";
printf "MEASUREMENT m_KmKmKpNu   statistical  systematic\n";
printf "DATA        m_KmKmKpNu   statistical  systematic\n";
printf "            3.29000e-05  %.5e  0.00000e+00\n",2.59226e-06*$sfac4;
printf "\n";
printf "STAT_CORR_WITH Belle PimPimPipNu submitted -0.05346557\n";
printf "STAT_CORR_WITH Belle PimKmPipNu submitted 0.03505142\n";
printf "STAT_CORR_WITH Belle PimKmKpNu submitted -0.008312885\n";
printf "\n";
printf "END\n";
