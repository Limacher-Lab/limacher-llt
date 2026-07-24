# Ensure LaTeX can find shared files in ../_common/
$ENV{'TEXINPUTS'} = '../_common//:' . ($ENV{'TEXINPUTS'} || '');
$ENV{'BSTINPUTS'} = '../_common//:' . ($ENV{'BSTINPUTS'} || '');