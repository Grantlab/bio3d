library(biomaRt)


listMarts()

msd <- useMart("msd")
listDatasets(msd)


msd <- useMart("msd", dataset="msd")


filters = listFilters(msd)
attributes = listAttributes(msd)


getat <- c("id", "pdb_id", "chain_code",
           "sun_id",
           "scientific_name",
           "hetgroup_id","hetgroup_name",
           "uniprot_id", "pfam_id",
           "year")

getBM(attributes = getat,
      filters = "pdb_id_list",
      values = c("1bg2","2kin","4q21"),
      mart = msd)
