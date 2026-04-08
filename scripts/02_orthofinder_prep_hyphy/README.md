
### **OrthoFinder run log:**

### 2023-09-28
Initial run with 102 species (101 spiders + _Drosophila_).

```bash
orthofinder \
	-f PATH_TO_ORIG_FASTAS_DIR/ \
	-a 36 -t 36
```

### 2023-11-27
Added in _Antrodiaetus_roretzi_ (mygalomorph); started analysis from BLAST results in RESULTS_SEP28 and added one FASTA file in ADDL_FASTAS_DIR

```bash
orthofinder \
	-b PATH_TO_RESULTS_SEP28/ \
	-f PATH_TO_ADDL_FASTAS_DIR/ \
	-a 35 -t 35 -M msa
```

### 2024-05-12
The tree from the above OrthoFinder runs was edited to concur with Kulkarni et al. 2023. OrthoFinder was re-run from the beginning, using all species and an edited species tree.

```bash
orthofinder \
	-f PATH_TO_FINAL_FASTAS_DIR/ \
	-s PATH_TO_EDITED_ORTHO_TREE \
	-a 36 -t 36 -y 
```

### 2024-06-03
Restarted from BLAST results step (`-b`) after timeout.

```bash
orthofinder \
	-b PATH_TO_RESULTS_MAY12/ \
	-s PATH_TO_EDITED_ORTHO_TREE \
	-a 34 -t 34 -y
```

### 2024-06-04
Restart from orthogroups directory (`-fg`) after "too many files" error.

```bash
su USERNAME

ulimit -n 100000

orthofinder \
	-fg PATH_TO_RESULTS_JUN03 \
	-s PATH_TO_EDITED_ORTHO_TREE \
	--fewer-files -t 70 -y 
```

### 2025-01-06
Final re-run from tree inference step (`-ft`) with finalized species tree.

```bash
orthofinder \
	-ft PATH_TO_RESULTS_JUN04/ \
	-s PATH_TO_FINAL_SPECIES_TREE \
	--fewer-files -t 70 -y 
```
