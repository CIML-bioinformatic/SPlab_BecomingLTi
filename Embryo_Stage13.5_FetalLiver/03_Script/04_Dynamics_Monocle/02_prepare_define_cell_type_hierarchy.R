
## @knitr define_cell_type_hierarchy

# Define the markers genes usd to classify the cells

Adgre1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Adgre1"))
Ccr6_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Ccr6"))
Cd14_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Cd14"))
Cd163_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Cd163"))
Cd4_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Cd4"))
Csf1r_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Csf1r"))
Cx3cr1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Cx3cr1"))
Cxcr4_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Cxcr4"))
Cxcr5_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Cxcr5"))
Cxcr6_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Cxcr6"))
Eomes_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Eomes"))
Fcgr2b_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Fcgr2b"))
Flt3_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Flt3"))
Gata3_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Gata3"))
Icos_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Icos"))
Id2_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Id2"))
Il13_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Il13"))
Il17rb_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Il17rb"))
Il18r1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Il18r1"))
Il1rl1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Il1rl1"))
Il22_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Il22"))
Il2ra_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Il2ra"))
Il2rb_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Il2rb"))
Il23r_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Il23r"))
Il33_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Il33"))
Il7r_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Il7r"))
Itga4_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Itga4"))
Itgb1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Itgb1"))
Itgb7_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Itgb7"))  
Klrg1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Klrg1"))
Lyve1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Lyve1"))
Lta_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Lta"))
Ltb_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Ltb"))
Ly6a_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Ly6a"))
Ncf1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Ncf1"))
Nfil3_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Nfil3"))
Psd2_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Psd2"))
Plekhs1_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Plekhs1"))
Rora_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Rora"))
Rorc_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Rorc"))
S100a8_id <- row.names(subset(fData( filtered_cds), gene_short_name == "S100a8"))
S100a9_id <- row.names(subset(fData( filtered_cds), gene_short_name == "S100a9"))
Siglecf_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Siglecf"))
Tbx21_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Tbx21"))
Tcf7_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Tcf7"))
Tnfsf11_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Tnfsf11"))
Tnfrsf11a_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Tnfrsf11a"))
Tox_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Tox"))
Zbtb16_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Zbtb16"))
Zfp541_id <- row.names(subset(fData( filtered_cds), gene_short_name == "Zfp541"))


# #####################################################################################
#  DIFFERENTIAL cell type hierarchy
# #####################################################################################

cat("<HR><BR<BR<BR<BR>Definition of cth_DIFFERENTIAL cell type hierarchy:")

cth_DIFFERENTIAL =  newCellTypeHierarchy()
cth_type_list = list()

cat("<BR>* CLP = Flt3+, Id2-")
cth_DIFFERENTIAL <- addCellType( cth_DIFFERENTIAL, "CLP", 
                             classify_func = function(x) { 
                               x[Flt3_id,] >= 1 & x[Id2_id,] < 1})
cth_type_list[[ "CLP"]] =  newCellTypeHierarchy()
cth_type_list[[ "CLP"]] <- addCellType( cth_type_list[[ "CLP"]], "CLP", 
                             classify_func = function(x) { 
                               x[Flt3_id,] >= 1 & x[Id2_id,] < 1})


cat("<BR>* AlphaLP = Id2+, Itga4+, Itgb7+, Rorc-")
cth_DIFFERENTIAL <- addCellType( cth_DIFFERENTIAL, "AlphaLP", 
                             classify_func = function(x) { 
                               x[ Id2_id,] >= 1 & x[ Itga4_id,] >= 1 & x[ Itgb7_id,] >= 1 &
                               x[ Rorc_id,] < 1})
cth_type_list[[ "AlphaLP"]] =  newCellTypeHierarchy()
cth_type_list[[ "AlphaLP"]] <- addCellType( cth_type_list[[ "AlphaLP"]], "AlphaLP", 
                             classify_func = function(x) { 
                               x[ Id2_id,] >= 1 & x[ Itga4_id,] >= 1 & x[ Itgb7_id,] >= 1 &
                               x[ Rorc_id,] < 1})


cat("<BR>* LTiP = Cxcr6+, Rorc+, Zbtb16-, Flt3-, Id2+, Tox+")
cth_DIFFERENTIAL <- addCellType( cth_DIFFERENTIAL, "LTiP", 
                                 classify_func = function(x) { 
                                   x[ Cxcr6_id,] >= 1 & x[ Rorc_id,] >= 1 & x[ Zbtb16_id,] < 1 & x[ Flt3_id,] < 1 &
                                   x[ Id2_id,] >= 1 & x[ Tox_id,] >= 1})
cth_type_list[[ "LTiP"]] =  newCellTypeHierarchy()
cth_type_list[[ "LTiP"]] <- addCellType( cth_type_list[[ "LTiP"]], "LTiP", 
                                 classify_func = function(x) { 
                                   x[ Cxcr6_id,] >= 1 & x[ Rorc_id,] >= 1 & x[ Zbtb16_id,] < 1 & x[ Flt3_id,] < 1 &
                                   x[ Id2_id,] >= 1 & x[ Tox_id,] >= 1})


cat("<BR>* LTi = Lta+, Ccr6+, Cd4+, Il7r+, Ltb+, Zbtb16-, Flt3-, Id2+")
cth_DIFFERENTIAL <- addCellType( cth_DIFFERENTIAL, "LTi", 
                                 classify_func = function(x) { 
                                   x[ Lta_id,] >= 1 & x[ Ccr6_id,] >= 1 & x[ Cd4_id,] >= 1 &
				   x[ Il7r_id,] >= 1 & x[ Ltb_id,] >= 1 & x[ Zbtb16_id,] < 1 &
 				   x[ Flt3_id,] < 1 & x[ Id2_id,] >= 1})
cth_type_list[[ "LTi"]] =  newCellTypeHierarchy()
cth_type_list[[ "LTi"]] <- addCellType( cth_type_list[[ "LTi"]], "LTi", 
                                 classify_func = function(x) { 
                                   x[ Lta_id,] >= 1 & x[ Ccr6_id,] >= 1 & x[ Cd4_id,] >= 1 &
				   x[ Il7r_id,] >= 1 & x[ Ltb_id,] >= 1 & x[ Zbtb16_id,] < 1 &
 				   x[ Flt3_id,] < 1 & x[ Id2_id,] >= 1})


cat("<BR>* ILCP = Zbtb16+, Rorc-, Ccr6-, Tox+")
cth_DIFFERENTIAL <- addCellType( cth_DIFFERENTIAL, "ILCP", 
                                 classify_func = function(x) { 
                                   x[ Zbtb16_id,] >= 1 & x[ Rorc_id,] < 1 & x[ Ccr6_id,] < 1 &
				   x[ Tox_id,] >= 1})
cth_type_list[[ "ILCP"]] =  newCellTypeHierarchy()
cth_type_list[[ "ILCP"]] <- addCellType( cth_type_list[[ "ILCP"]], "ILCP", 
                        classify_func = function(x) { 
                                   x[ Zbtb16_id,] >= 1 & x[ Rorc_id,] < 1 & x[ Ccr6_id,] < 1 &
				   x[ Tox_id,] >= 1})


cat("<BR>* ILC1 = Il2rb+, Rorc-, Tbx21+")
cth_DIFFERENTIAL <- addCellType( cth_DIFFERENTIAL, "ILC1", 
                                 classify_func = function(x) { 
                                   x[ Il2rb_id,] >= 1 & x[ Rorc_id,] < 1 & x[ Tbx21_id,] >= 1})
cth_type_list[[ "ILC1"]] =  newCellTypeHierarchy()
cth_type_list[[ "ILC1"]] <- addCellType( cth_type_list[[ "ILC1"]], "ILC1", 
                        classify_func = function(x) { 
                                   x[ Il2rb_id,] >= 1 & x[ Rorc_id,] < 1 & x[ Tbx21_id,] >= 1})


cat("<BR>* ILC2 = Il13+, Rora+, Rorc-, Il17rb+, Icos+")
cth_DIFFERENTIAL <- addCellType( cth_DIFFERENTIAL, "ILC2", 
                                 classify_func = function(x) { 
                                   x[ Il13_id,] >= 1 & x[ Rora_id,] >= 1 & x[ Rorc_id,] < 1 &
				   x[ Il17rb_id,] >= 1 & x[ Icos_id,] >= 1})
cth_type_list[[ "ILC2"]] =  newCellTypeHierarchy()
cth_type_list[[ "ILC2"]] <- addCellType( cth_type_list[[ "ILC2"]], "ILC2", 
                         classify_func = function(x) { 
                                   x[ Il13_id,] >= 1 & x[ Rora_id,] >= 1 & x[ Rorc_id,] < 1 &
				   x[ Il17rb_id,] >= 1 & x[ Icos_id,] >= 1})


cat("<BR>* ILC3 = Rorc+, Lta-, Cd4-, Il22+, Zbtb16+")
cth_DIFFERENTIAL <- addCellType( cth_DIFFERENTIAL, "ILC3", 
                                 classify_func = function(x) { 
                                   x[ Rorc_id,] >= 1 & x[ Lta_id,] < 1 & x[ Cd4_id,] < 1 & x[ Il22_id,] >= 1 & x[ Zbtb16_id,] >= 1})
cth_type_list[[ "ILC3"]] =  newCellTypeHierarchy()
cth_type_list[[ "ILC3"]] <- addCellType( cth_type_list[[ "ILC3"]], "ILC3", 
                                 classify_func = function(x) { 
                                   x[ Rorc_id,] >= 1 & x[ Lta_id,] < 1 & x[ Cd4_id,] < 1 & x[ Il22_id,] >= 1 & x[ Zbtb16_id,] >= 1})



# #####################################################################################
#  ORIGINAL cell type hierarchy
# #####################################################################################
cat("<BR<BR<BR<BR>Definition of cth_ORIGINAL cell type hierarchy:")

cth_ORIGINAL =  newCellTypeHierarchy()
cth_type_list = list()

cat("<BR>* CLP = Flt3+ Id2-")
cth_ORIGINAL <- addCellType( cth_ORIGINAL, "CLP", 
                                 classify_func = function(x) { 
                                   x[Flt3_id,] >= 1 & x[Id2_id,] < 1})
cth_type_list[[ "CLP"]] =  newCellTypeHierarchy()
cth_type_list[[ "CLP"]] <- addCellType( cth_type_list[[ "CLP"]], "CLP", 
                                 classify_func = function(x) { 
                                   x[Flt3_id,] >= 1 & x[Id2_id,] < 1})


cat("<BR>* AlphaLP = Id2+, Itga4+, Itgb1+, Rorc-")
cth_ORIGINAL <- addCellType( cth_ORIGINAL, "AlphaLP", 
                                 classify_func = function(x) { 
                                   x[ Id2_id,] >= 1 & x[ Itga4_id,] >= 1 & x[ Itgb1_id,] >= 1 &
                                   x[ Rorc_id,] < 1})
cth_type_list[[ "AlphaLP"]] =  newCellTypeHierarchy()
cth_type_list[[ "AlphaLP"]] <- addCellType( cth_type_list[[ "AlphaLP"]], "AlphaLP", 
                                 classify_func = function(x) { 
                                   x[ Id2_id,] >= 1 & x[ Itga4_id,] >= 1 & x[ Itgb1_id,] >= 1 &
                                   x[ Rorc_id,] < 1})


cat("<BR>* LTiP = Cxcr6+, Rorc+, Zbtb16-, Flt3- Id2+")
cth_ORIGINAL <- addCellType( cth_ORIGINAL, "LTiP", 
                                 classify_func = function(x) { 
                                   x[ Cxcr6_id,] >= 1 & x[ Rorc_id,] >= 1 & x[ Zbtb16_id,] < 1 & x[ Flt3_id,] < 1 &
                                   x[ Id2_id,] >= 1})
cth_type_list[[ "LTiP"]] =  newCellTypeHierarchy()
cth_type_list[[ "LTiP"]] <- addCellType( cth_type_list[[ "LTiP"]], "LTiP", 
                                  classify_func = function(x) { 
                                   x[ Cxcr6_id,] >= 1 & x[ Rorc_id,] >= 1 & x[ Zbtb16_id,] < 1 & x[ Flt3_id,] < 1 &
                                   x[ Id2_id,] >= 1})


cat("<BR>* LTi = Lta+, Ccr6+, Cd4+")
cth_ORIGINAL <- addCellType( cth_ORIGINAL, "LTi", 
                                 classify_func = function(x) { 
                                   x[ Lta_id,] >= 1 & x[ Ccr6_id,] >= 1 & x[ Cd4_id,] >= 1})
cth_type_list[[ "LTi"]] =  newCellTypeHierarchy()
cth_type_list[[ "LTi"]] <- addCellType( cth_type_list[[ "LTi"]], "LTi", 
                                 classify_func = function(x) { 
                                   x[ Lta_id,] >= 1 & x[ Ccr6_id,] >= 1 & x[ Cd4_id,] >= 1})

cat("<BR>* ILCP = Zbtb16+, Rorc-, Ccr6-")
cth_ORIGINAL <- addCellType( cth_ORIGINAL, "ILCP", 
                                 classify_func = function(x) { 
                                   x[ Zbtb16_id,] >= 1 & x[ Rorc_id,] < 1 & x[ Ccr6_id,] < 1})
cth_type_list[[ "ILCP"]] =  newCellTypeHierarchy()
cth_type_list[[ "ILCP"]] <- addCellType( cth_type_list[[ "ILCP"]], "ILCP", 
                                 classify_func = function(x) { 
                                   x[ Zbtb16_id,] >= 1 & x[ Rorc_id,] < 1 & x[ Ccr6_id,] < 1})


cat("<BR>* ILC1 = Il2rb+, Rorc-, Eomes+")
cth_ORIGINAL <- addCellType( cth_ORIGINAL, "ILC1", 
                             classify_func = function(x) { 
                               x[ Il2rb_id,] >= 1 & x[ Rorc_id,] < 1 & x[ Eomes_id,] >= 1})
cth_type_list[[ "ILC1"]] =  newCellTypeHierarchy()
cth_type_list[[ "ILC1"]] <- addCellType( cth_type_list[[ "ILC1"]], "ILC1", 
                                         classify_func = function(x) { 
                                           x[ Il2rb_id,] >= 1 & x[ Rorc_id,] < 1 & x[ Eomes_id,] >+ 1})


cat("<BR>* ILC2 = Il13+, Rora+, Rorc-")
cth_ORIGINAL <- addCellType( cth_ORIGINAL, "ILC2", 
                                 classify_func = function(x) { 
                                   x[ Il13_id,] >= 1 & x[ Rora_id,] >= 1 & x[ Rorc_id,] < 1})
cth_type_list[[ "ILC2"]] =  newCellTypeHierarchy()
cth_type_list[[ "ILC2"]] <- addCellType( cth_type_list[[ "ILC2"]], "ILC2", 
                                 classify_func = function(x) { 
                                   x[ Il13_id,] >= 1 & x[ Rora_id,] >= 1 & x[ Rorc_id,] < 1})


cat("<BR>* ILC3 = Rorc+, Lta-, Cd4- Il22+")
cth_ORIGINAL <- addCellType( cth_ORIGINAL, "ILC3", 
                                 classify_func = function(x) { 
                                   x[ Rorc_id,] >= 1 & x[ Lta_id,] < 1 & x[ Cd4_id,] < 1 & x[ Il22_id,] >= 1})
cth_type_list[[ "ILC3"]] =  newCellTypeHierarchy()
cth_type_list[[ "ILC3"]] <- addCellType( cth_type_list[[ "ILC3"]], "ILC3", 
                                 classify_func = function(x) { 
                                   x[ Rorc_id,] >= 1 & x[ Lta_id,] < 1 & x[ Cd4_id,] < 1 & x[ Il22_id,] >= 1})


