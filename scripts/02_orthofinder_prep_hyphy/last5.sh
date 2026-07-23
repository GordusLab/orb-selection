WD=/scratch4/agordus1/crunnel2/hyphy_wd
HOG_LIST=/home/crunnel2/orb-selection/data/N5.udiv.o75_list.txt
FG_LIST=/home/crunnel2/orb-selection/data/orbweavers-list.txt
FG_NAME=orb_fg

for NUM in 1050 1051 1048 1045 1047; do
    CURRENT_HOG=$(sed "${NUM}q;d" $HOG_LIST)
    echo "Working on ${CURRENT_HOG}..."
    FLTRD_TREE=${WD}/${CURRENT_HOG}.fltrd.tree

    #make a copy of the filtered tree to iteratively label
    LBLD_TREE=${WD}/${CURRENT_HOG}.${FG_NAME}.tree
    if [ -f ${LBLD_TREE} ]; then
        echo "The tree ${LBLD_TREE} already exists."
    else
        cp ${FLTRD_TREE} ${LBLD_TREE}

        while read p; do
            REGEX="^${p}"
            hyphy /home/crunnel2/bin/hyphy-analyses/LabelTrees/label-tree.bf \
            --tree ${LBLD_TREE} \
            --regexp $REGEX \
            --output ${LBLD_TREE}
        done < ${FG_LIST}

        #Add a semicolon to the end of the last line of the tree file to avoid errors in HyPhy
        sed -i '$s/$/;/' ${LBLD_TREE}

        mv ${LBLD_TREE} ${LBLD_TREE}.tmp

        # LabelTrees is adding 1E-10 branch lengths to every branch for some reason
        gotree brlen clear -i ${LBLD_TREE}.tmp -o ${LBLD_TREE} && rm ${LBLD_TREE}.tmp

    fi
    echo "Finished labeling ${CURRENT_HOG}."
done
