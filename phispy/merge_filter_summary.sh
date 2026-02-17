(
  echo -e "idx\tgenome\tprophage\tstatus\tnote\tcontig\tlength\tIS_bp\tIS_pct"

  find /scratch/cw106/phispy_output -name '*_islen.tsv' | \
  while read file; do
    genome=$(basename "$(dirname "$file")")
    idx=$(grep -n -m1 -Fx "$genome" /projects/alvarez/genomes.list | cut -d: -f1)

    if [[ -z $idx ]]; then
      base=${genome%_genomic}
      idx=$(grep -n -m1 -Fx "$base" /projects/alvarez/genomes.list | cut -d: -f1)
    fi

    awk -F'\t' -v g="$genome" -v i="$idx" 'NR>1{print i"\t"g"\t"$0}' "$file"
  done
) > /scratch/cw106/all_filter_summary_islen.tsv
