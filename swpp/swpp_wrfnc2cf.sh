#! /bin/bash -x

use python2.7.2

function usage(){
  cat << End_of_help
    Usage: swpp_wrfnc2cf [options]

    --variables VAR1,VAR2,...,VARn
             -v VAR1,VAR2,...,VARn
      Comma-separatd list of variables to process
      
    --help, -h
      Print this help message
End_of_help
}

#####################################################################
#
#  Option parser
#
#  Try to get the configuration from the environment and the command
#  line.
#
#####################################################################
filetype="${SWPP_FILETYPE:-out}"
geofile=${SWPP_GEOFILE}
fullfile=${SWPP_FULLFILE}
wrfncxnj=${SWPP_WRFNCXNJ_PY}
xnjattr=${SWPP_WRFNCXNJ_GLOBALATTR}
xnjtable=${SWPP_WRFNCXNJ_TABLE}
outpattern=${SWPP_OUTPUT_PATTERN}
tmpdirbase=${SWPP_TEMPORARY_STORAGE}
checknofim=${SWPP_CHECK_NOFIM:-1} # Check Number Of Files In Month

while test -n "$1"; do
  case "$1" in
    --variables|-v) varlist=$2; shift ;;
    --levels|-l) levlist=$2; shift ;;
    --slevels|-s) slevlist=$2; shift ;;
    --framesperblock) echo "Option not implemented, yet: $1" ;;
    --maxblocks|-m) echo "Option not implemented, yet: $1" ;;
    --spinupframes) echo "Option not implemented, yet: $1" ;;
    --outpattern|-o) outpattern=$2; shift ;;
    --inputpath|-i) inputpath=$2; shift ;;
    --temporary-storage|-t) tmpdirbase=$2; shift ;;
    --geofile|-g) geofile=$2; shift ;;
    --fullfile|-f) fullfile=$2; shift ;;
    --filetype|-y) filetype=$2; shift ;;
    --wrfncxnj|-w) wrfncxnj=$2; shift ;;
    --dont-check-nofim) checknofim=0 ;;
    --help|-h) usage; exit ;;
    *) echo "Unknown option: $1"; exit ;;
  esac
  shift
done
#####################################################################
#
#  Functions
#
#####################################################################

function print_config(){
  # Prints the configuration options taken from the environment
  # or command line
  cat << End_of_conf
  SWPP_OUTPUT_PATTERN (--outpattern)
    ${outpattern}
  SWPP_TEMPORARY_STORAGE (--temporary-storage)
    ${tmpdirbase}
  geofile = ${geofile}
  fullfile = ${fullfile}
  wrfncxnj = ${wrfncxnj}
  xnjtable = ${xnjtable}
  xnjattr = ${xnjattr}
End_of_conf
}

function ipcc_varname() {
  awk '$1=="'$1'"{print $2}' $2
}

function filelist_cache(){
  idir=$1
  filetype=$2
  find ${idir}/ -name wrf${filetype}_d01_\*.nc | sort > cache.filelist
}

function get_yearmons(){
  idir=$1
  filetype=$2
  test -e cache.filelist || filelist_cache ${idir} ${filetype}
  cat cache.filelist \
    | sed -e 's/^.*wrf'${filetype}'_d01_\(.*\)..T......Z_\(.*\)..T......Z.nc$/\1/' \
    | sort -n | uniq
}

function get_filelist_yearmon(){
  exppath=$1
  yearmon=$2
  filetype=$3
  test -e cache.filelist || filelist_cache ${idir} ${filetype}
  # TODO: Esta forma de coger los ficheros es peligrosa. Ejemplo:
  # Buscar ficheros con yearmon=200101 cogeria un fichero de 20200101T...
  # que no le corresponde.
  grep -E "wrf${filetype}_d01_.*${yearmon}.*\.nc" cache.filelist | sort
}

function get_previousfile_yearmon(){
  exppath=$1
  yearmon=$2
  filetype=$3
  prevyearmon=$(date --utc '+%Y%m' -d "${yearmon:0:4}-${yearmon:4:2}-15 -1 month")
  get_filelist_yearmon ${exppath} ${prevyearmon} ${filetype} | tail -1
}

function get_nextfile_yearmon(){
  exppath=$1
  yearmon=$2
  filetype=$3
  nextyearmon=$(date --utc '+%Y%m' -d "${yearmon:0:4}-${yearmon:4:2}-15 1 month")
  get_filelist_yearmon ${exppath} ${nextyearmon} ${filetype} | head -1
}

function filelist_is_short(){
  flfile=$1
  ymon=$2
  s1=$(date --utc '+%s' -d "${ymon:0:4}-${ymon:4:2}-01 1 month")
  s2=$(date --utc '+%s' -d "${ymon:0:4}-${ymon:4:2}-01")
  daysinmonth=$(python -c "print ($s1-$s2)/60/60/24")
  test "${daysinmonth}" -ne $(cat ${flfile} | wc -l)
}

function valid_dates(){
  ymon=$1
  s1=$(date --utc '+%Y%m%d%H' -d "${ymon:0:4}-${ymon:4:2}-01")
  s2=$(date --utc '+%Y%m%d%H' -d "${ymon:0:4}-${ymon:4:2}-01 1 month")
  echo "${s1},${s2}"
}

function fix_mon_xtrm(){
  xtrmfile=$1
  mon=$2
  echo "Removing unwanted timesteps from ${xtrmfile}"
  cdo selmon,${mon} ${xtrmfile} ${tmpdir}/$(basename ${xtrmfile})
  mv ${tmpdir}/$(basename ${xtrmfile}) ${xtrmfile}
}

#####################################################################
#
#  Main program
#
#####################################################################

scriptdir=$( (cd `dirname $0` && echo $PWD) )
thisdir=$(pwd)
#
#  Get a private space to run
#
if test -n "${tmpdirbase}"; then
  tmpdir=${tmpdirbase}/tmpdir.`date +%Y%m%d%H%M%S%n`
else
  tmpdir=${thisdir}/tmpdir.`date +%Y%m%d%H%M%S%n`
fi
mkdir -p ${tmpdir} && cd ${tmpdir} || exit

wxajcmd="python ${wrfncxnj} -a ${xnjattr} -g ${geofile} -r 1950-01-01_00:00:00 --fullfile ${fullfile} --temp-dir ${tmpdir}/xnj" 

print_config

for yearmon in $(get_yearmons ${inputpath} ${filetype}); do
  year=${yearmon:0:4}
  month=${yearmon:4:2}
  outfiles=$(eval echo $outpattern) # expand year and month into the pattern
  rm -f filelist.txt
  get_filelist_yearmon ${inputpath} ${yearmon} ${filetype} \
    | while read f; do
#     if ! ncdump -h ${f} | grep -q num_metgrid_levels ; then
#       echo "$f was not filtered at running time. I'm aborting the processing of ${yearmon}"
#       echo "Run p_interp before extracting variables!"
#       continue
#     else
      echo ${f} >> filelist.txt
#     fi
  done
  if test ${checknofim} -eq 1; then
    if filelist_is_short filelist.txt ${yearmon}; then
      echo "[${yearmon}] Too few files for month ${month}: $(wc -l filelist.txt)"
      continue
    fi
  fi
  xnjflags=""
  prevfile=$(get_previousfile_yearmon ${inputpath} ${yearmon} ${filetype})
  nextfile=$(get_nextfile_yearmon     ${inputpath} ${yearmon} ${filetype})
  test -n "${prevfile}" && xnjflags="${xnjflags} --previous-file ${prevfile}"
  test -n "${nextfile}" && xnjflags="${xnjflags} --next-file ${nextfile}"
  test -n "${levlist}" && xnjflags="${xnjflags} --split-levels --plevs-filter ${levlist}"
  test -n "${slevlist}" && xnjflags="${xnjflags} --slevs-filter ${slevlist}"
  test "${filetype}" == "xtrm" && xnjflags="${xnjflags} --filter-times=$(valid_dates ${yearmon})"
  test -d $(dirname ${outfiles}) || mkdir -p  $(dirname ${outfiles})
  ${wxajcmd} \
    --from-file filelist.txt \
    --output-pattern ${outfiles} \
    -v ${varlist} \
    --split-variables \
    --output-format="NETCDF4_CLASSIC" \
    ${xnjflags} \
    >& xnj${yearmon}.log
  mv filelist.txt filelist_${yearmon}.txt
#  if [[ "${filetype}" == "xtrm" ]]
#  then
#    fix_mon_xtrm "$(echo ${outfiles} | sed s/'\[varcf\]'/tasmax/ | sed s/'\[level\]'//)" ${month}
#    fix_mon_xtrm "$(echo ${outfiles} | sed s/'\[varcf\]'/tasmin/ | sed s/'\[level\]'//)" ${month}
#  fi
done # yearmon
