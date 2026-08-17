ARGIND==1 {
  line=$0
  gsub(/[ \t]*\|[ \t]*/,"\t", line)
  split(line,f,"\t")
  par[f[1]+0]   = f[2]+0
  trank[f[1]+0] = f[3]
  next
}

ARGIND==2 {
  line=$0
  gsub(/[ \t]*\|[ \t]*/,"\t", line)
  split(line,f,"\t")
  if (f[4]=="scientific name") taxname[f[1]+0] = f[2]
  next
}

{
  base = 0.2

  status=$1
  qid=$2
  assigned=$3+0

  # Rebuild full kmer list as ONE string, then split on spaces
  klist = ""
  for (i=5; i<=NF; i++) {
    if ($i != "") klist = klist (klist=="" ? "" : " ") $i
  }

  if (klist=="") {
    printf "%s\t%d\t%s\t%s\t%.2f\n", qid, assigned, "unclassified", "unclassified", 0
    next
  }

  tot_all=0
  zero=0
  delete clade

  n = split(klist, a, /[ ]+/)
  for (j=1; j<=n; j++) {
    if (a[j]=="" || a[j] !~ /:/) continue

    split(a[j], kv, ":")
    t = kv[1]+0
    c = kv[2]+0

    tot_all += c
    if (t==0) { zero += c; continue }

    cur = t
    while (cur > 0) {
      clade[cur] += c
      if (cur==1) break
      cur = par[cur]
      if (cur==0 || cur=="") break
    }
  }

  classified = tot_all - zero

  # Cap host resolution at species: if Kraken assigned a rank finer than species
  # (strain / subspecies / no-rank under a species), climb to the species ancestor
  # and report that. Assignments already at species or coarser (genus and up) are
  # kept unchanged.
  hostid = assigned
  up = assigned
  while (up > 1 && up != "" && up != 0) {
    if (trank[up] == "species") { hostid = up; break }
    up = par[up]
  }

  support = (hostid in clade ? clade[hostid] : 0)

  frac = (classified>0 ? support/classified : 0)
  if (frac > 1) frac = 1

  score = (classified>0 ? base + (1-base)*frac : 0)
  if (score > 1) score = 1

  name = (hostid in taxname ? taxname[hostid] : "NA")
  rank = (hostid in trank ? trank[hostid] : "NA")

  # HRGMv2 carries GTDB nomenclature: a rank prefix ("t__", "s__", "g__" ...) and, on leaf
  # taxa, a "/HRGMv2_NNNN" genome id — so a raw host reads "t__Megamonas funiformis/HRGMv2_3482"
  # instead of "Megamonas funiformis". This column is user-facing, so strip both.
  # A NO-OP on hgdb_unphaged, whose names are plain NCBI ("Megamonas") and match neither
  # pattern, which is why this is unconditional rather than switched on the database in use.
  sub(/^[a-z]__/, "", name)
  sub(/\/HRGMv2_[0-9]+$/, "", name)

  if (status=="U" || assigned==0) { score=0; name="unclassified"; rank="unclassified"; hostid=assigned }

  printf "%s\t%d\t%s\t%s\t%.2f\n", qid, hostid, rank, name, score
}

