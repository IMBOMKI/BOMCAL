# echo "Setting BOMKI_analysis v999 in /home/bomki/ICEDUST"

if test "${CMTROOT}" = ""; then
  CMTROOT=/home/bomki/ICEDUST/CMT/v1r20p20081118; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh

tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=BOMKI_analysis -version=v999 -path=/home/bomki/ICEDUST  -no_cleanup $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

