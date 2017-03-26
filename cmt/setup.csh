# echo "Setting BOMKI_analysis v999 in /home/bomki/ICEDUST"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /home/bomki/ICEDUST/CMT/v1r20p20081118
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=BOMKI_analysis -version=v999 -path=/home/bomki/ICEDUST  -no_cleanup $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

