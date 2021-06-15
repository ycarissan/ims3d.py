NICS_PREP_DIR = $HOME/prog/prepare_nics
for i in $*
do
        dir=${i/.xyz/}
        mkdir ${dir}
        cd ${dir}
            cp ../$i .
            python3 ${NICS_PREP_DIR}/nics_prep_angular.py $i -r 1 --depth 4
            subg09 *.com
        cd ../
done
