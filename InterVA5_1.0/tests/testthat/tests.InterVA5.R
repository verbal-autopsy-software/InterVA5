context("InterVA5 Algorithm")

# read in the dummydata created by InterVA5_dummy_v2.R (also in the test directory)
# dummy data is not created by the script here in case the probbase file changes unintentionally
load("InterVA5_dummy_v2.RData")

# A group
for (cod in rownames(dummydata)[grep("a_", rownames(dummydata))]){
  test_that(paste0(cod, " (", dummydata[cod, "ID"], ")"), {

    out <- InterVA5(dummydata[cod,, drop=FALSE], HIV = "l", Malaria = "l", write=F)

    expect_true(dummydata[cod, "ID"] %in% c(out$VA5[[1]]$PREGSTAT))
    rm(out)
  })
}

# B group
for (cod in rownames(dummydata)[grep("b_", rownames(dummydata))]){
  test_that(paste0(cod, " (", dummydata[cod, "ID"], ")"), {

    out <- InterVA5(dummydata[cod,, drop=FALSE], HIV = "l", Malaria = "l", write=F)

    expect_true(dummydata[cod, "ID"] %in% c(out$VA5[[1]]$CAUSE1))
    rm(out)
  })
}

# C group
for (cod in rownames(dummydata)[grep("c_", rownames(dummydata))]){
  test_that(paste0(cod, " (", dummydata[cod, "ID"], ")"), {

    out <- InterVA5(dummydata[cod,, drop=FALSE], HIV = "l", Malaria = "l", write=F)

    expect_true(dummydata[cod, "ID"] %in% c(out$VA5[[1]]$COMCAT))
    rm(out)
  })
}
