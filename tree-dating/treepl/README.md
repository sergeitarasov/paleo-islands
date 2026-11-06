# Individual analyses


/treepl_run  treepl prime.ini > log/prime.log
/treepl_run  treepl cv.ini > log/cv.log
/treepl_run  treepl final.ini > log/final.log
/treepl_run  treepl cv.ini > log/cv.log
/treepl_run  treepl cv.ini > log/cv.log

---

# Batch analyses

chmod +x run-ini.sh
./run-ini.sh


## unparallel
./run-ini.sh  113.33s user 1.51s system 99% cpu 1:55.01 total
## parallel
./run-ini.sh  125.12s user 0.87s system 429% cpu 29.304 total

# Merge

cat timetrees/*.tre > summary-tree/treepl_merged.tre

---

# Next use tree annotator on merged
Create MCC tree with median heights

treeannotator \
  -burnin 0 \
  -heights median \
  summary-tree/treepl_merged.tre \
  summary-tree/treepl-316-mcc_median.nex

## Create MCC tree with common ancestor heights
treeannotator \
  -burnin 0 \
  -heights ca \
  summary-tree/treepl_merged.tre \
  summary-tree/treepl-316-ca.tre

treeannotator \
  -burnin 0 \
  -heights mean \
  summary-tree/treepl_merged.tre \
  summary-tree/treepl-316-mcc_mean.tre

