#!/app/easybuild/software/R/3.5.1-foss-2016b-fh1/bin/Rscript


#Jenny Smith 
#8/28/18

#Purpose: Run GSEA on NUP98-KDM5A

library(methods)

setwd('/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/2017.07.13_NUP98KDM5A_DEGs/')
source("~/scripts/RNAseq_Analysis/DifferentialExpn_PathwayAnalysis/GAGE_GSEA_Function.r")

DEGs <- get(load("TARGET_AML_NUP98.KDM5A_vs_OtherAML_IncludingNUPs_DEGs_List.RData"))
# C2.KEGG <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c2.cp.kegg.v6.0.symbols.RDS")
filename <- "TARGET_AML_allCohorts_NUP98.KDM5A_vs_OtherAML"


# print("starting1")
# 
# GSA <- lapply(DEGs, gage_from_pipeline,
#               method="trend",
#               type="expn",
#               geneset=C2.KEGG)
# 
# save(GSA,file=paste0(filename, "_expn_C2.KEGG.RData"))
# rm(GSA)
# gc()
# 
# print("done1")
# 
# 
# print("starting2")
# 
# GSA.FC <- lapply(DEGs, gage_from_pipeline,
#                  method="trend",
#                  type="FC",
#                  geneset=C2.KEGG)
# save(GSA.FC,file=paste0(filename, "_FC_C2.KEGG.RData"))
# rm(GSA.FC)
# gc()
# 
# 
# rm(C2.KEGG)
# print("done2")
# 
# print("starting3")
# 
# 
# C2.All <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c2.all.v6.0.symbols.RDS")
# 
# GSA.C2.All <- lapply(DEGs, gage_from_pipeline, 
#                      method="trend",
#                      type="expn",
#                      geneset=C2.All)
# save(GSA.C2.All,file=paste0(filename, "_expn_C2.All.RData"))
# rm(GSA.C2.All)
# gc()
# 
# 
# print("done3")
# 
# 
# print("starting4")
# 
# GSA.KEGG <- lapply(DEGs, gage_from_pipeline,
#                    method="trend",
#                    type="expn",
#                    geneset=NULL)
# save(GSA.KEGG,file=paste0(filename, "_expn_HSA.KEGG.RData"))
# rm(GSA.KEGG)
# gc()
# 
# print("done4")
# 
# C5 <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c5_list_SetSize_GE.50_LE.300_v6.1.symbols.RDS")
# 
# print("starting5")
# GSA.GO.BioProcess <- lapply(DEGs, gage_from_pipeline,
#                             method="trend",
#                             type="expn", 
#                             geneset=C5[["c5.bp"]])
# save(GSA.GO.BioProcess, file=paste0(filename, "_expn_C5.BioProcess_SetSize50to300.RData"))
# rm(GSA.GO.BioProcess)
# gc()
# print("done5")
# 
# 
# 
# print("starting6")
# GSA.GO.CellComp <- lapply(DEGs, gage_from_pipeline, 
#                           method="trend",
#                           type="expn", 
#                           geneset=C5[["c5.cc"]])
# save(GSA.GO.CellComp, file=paste0(filename, "_expn_C5.CellComp_SetSize50to300.RData"))
# rm(GSA.GO.CellComp)
# gc()
# print("done6")    
# 
# 
# 
# print("starting7")
# GSA.GO.MolFunc <- lapply(DEGs, gage_from_pipeline, 
#                          method="trend",
#                          type="expn",
#                          geneset=C5[["c5.mf"]])
# save(GSA.GO.MolFunc, file=paste0(filename, "_expn_C5.MolFunc_SetSize50to300.RData"))
# rm(GSA.GO.MolFunc)
# gc()
# print("done7")
# 
# 
# print("starting8")
# GSA.GO.BioProcess.FC <- lapply(DEGs, gage_from_pipeline, 
#                                method="trend",
#                                type="FC", 
#                                geneset=C5[["c5.bp"]])
# save(GSA.GO.BioProcess.FC, file=paste0(filename, "_FC_C5.BioProcess_SetSize50to300.RData"))
# rm(GSA.GO.BioProcess.FC)
# gc()
# print("done8")
# 
# 
# 
# print("starting9")
# GSA.GO.CellComp.FC <- lapply(DEGs, gage_from_pipeline, 
#                              method="trend",
#                              type="FC", 
#                              geneset=C5[["c5.cc"]])
# save(GSA.GO.CellComp.FC, file=paste0(filename, "_FC_C5.CellComp_SetSize50to300.RData"))
# rm(GSA.GO.CellComp.FC)
# gc()
# print("done9")    
# 
# 
# 
# print("starting10")
# GSA.GO.MolFunc.FC <- lapply(DEGs, gage_from_pipeline,
#                             method="trend",
#                             type="FC", 
#                             geneset=C5[["c5.mf"]])
# save(GSA.GO.MolFunc.FC, file=paste0(filename, "_FC_C5.MolFunc_SetSize50to300.RData"))
# rm(GSA.GO.MolFunc.FC)
# gc()
# print("done10")

C3 <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c3.tft.v6.2.symbols.RDS")

print("starting11")
GSA.TFB <- lapply(DEGs, gage_from_pipeline,
                  method="trend",
                  type="expn", 
                  geneset=C3)
save(GSA.TFB, file=paste0(filename, "_expn_C3_TFbindingSites.RData"))
rm(GSA.TFB)
gc()
print("done11")


print("starting12")
GSA.TFB.FC <- lapply(DEGs, gage_from_pipeline,
                     method="trend",
                     type="FC", 
                     geneset=C3)
save(GSA.TFB.FC, file=paste0(filename, "_FC_C3_TFbindingSites.RData"))
rm(GSA.TFB.FC)
gc()
print("done12")


