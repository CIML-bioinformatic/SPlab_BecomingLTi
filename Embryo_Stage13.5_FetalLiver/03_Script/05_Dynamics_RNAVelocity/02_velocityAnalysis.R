# ####################################################
# This script aim to use the Velocyto R pipeline
# to analyse the RNA velocity on LysoDC data
# ####################################################

## @knitr estimate_velocity

# Take cluster labels from Pagoda2 pre-processing
#cluster.label <- r$clusters$PCA[[1]]
#cell.colors <- pagoda2:::fac2col(cluster.label)
cluster.label = as.factor( SEURAT_CLUSTERS[ colnames( emat), "Idents.sc10x.rna.seurat."])
names( cluster.label) = colnames( emat)
cell.colors.cluster = fac2col( cluster.label)

# Filter genes based on the minimum average expresion magnitude 
# (in at least one of the clusters), output total number of resulting valid genes
emat <- filter.genes.by.cluster.expression( emat, cluster.label, min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression( nmat, cluster.label, min.max.cluster.average = 0.05)
cat("<BR>Number of validated genes for analysis =",length(intersect(rownames(emat),rownames(nmat))), "<BR>")

# Take t-SNE embedding from Pagoda2 pre-processing
emb_tSNE <- as.matrix( SEURAT_TSNE_EMBEDDING_WITHCELLTYPE[ colnames( emat), c("tSNE_1", "tSNE_2")])
emb_UMAP <- as.matrix( SEURAT_UMAP_EMBEDDING_WITHCELLTYPE[ colnames( emat), c("UMAP_1", "UMAP_2")])

# In addition to clustering and the t-SNE embedding, from the pagoda2 processing 
# we will also take a cell-cell distance, which will be better than the default 
# whole-transcriptome correlation distance that velocyto.R would normally use.
#cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))
cell.dist = as.dist(1-armaCor( t( SEURAT_PCA_EMBEDDING[ colnames( emat), ])))

# Estimate RNA velocity (using gene-relative model with k=20 cell kNN pooling and 
# using top/bottom 2% quantiles for gamma fit)
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)

cat("<H4>Velocity flow on PCA map colored by cluster</H4>")
pca.velocity.plot( rvel.cd,
                   cell.colors=ac( cell.colors.cluster,alpha=0.5),
                   scale='sqrt',
                   cex=0.8,
                   arrow.scale=1,
                   show.grid.flow=TRUE,
                   min.grid.cell.mass=0.5,
                   grid.n=40,
                   arrow.lwd=1,
                   do.par=F,
                   cell.border.alpha = 0.1)


# Visualize velocity on the t-SNE embedding, using velocity vector fields
cat("<HR>")
par( mfrow = c(1,1))
# 
# cat("<H4>Velocity flow on t-SNE map</H4>")
# 
# show.velocity.on.embedding.cor(emb_tSNE,
#                                rvel.cd,
#                                n=100,
#                                scale='sqrt',
#                                cell.colors=ac(cell.colors.cluster,alpha=0.5),
#                                cex=0.8,
#                                arrow.scale=3,
#                                show.grid.flow=FALSE,
#                                min.grid.cell.mass=0.5,
#                                grid.n=40,
#                                arrow.lwd=1,
#                                do.par=F,
#                                cell.border.alpha = 0.1)
# 
# show.velocity.on.embedding.cor(emb_tSNE,
#                                rvel.cd,
#                                n=100,
#                                scale='sqrt',
#                                cell.colors=ac(cell.colors.cluster,alpha=0.5),
#                                cex=0.8,
#                                arrow.scale=3,
#                                show.grid.flow=TRUE,
#                                min.grid.cell.mass=0.5,
#                                grid.n=40,
#                                arrow.lwd=1,
#                                do.par=F,
#                                cell.border.alpha = 0.1)

cat("<HR>")
cat("<H4>Velocity flow on UMAP map with one arrow per cell and colored by cluster</H4>")

show.velocity.on.embedding.cor(emb_UMAP,
                               rvel.cd,
                               n=100,
                               scale='sqrt',
                               cell.colors=ac(cell.colors.cluster,alpha=0.5),
                               cex=0.8,
                               arrow.scale=3,
                               show.grid.flow=FALSE,
                               min.grid.cell.mass=0.5,
                               grid.n=40,
                               arrow.lwd=1,
                               do.par=F,
                               cell.border.alpha = 0.1)

cat("<H4>Velocity flow on UMAP map with grid flow and colored by cluster</H4>")

show.velocity.on.embedding.cor(emb_UMAP,
                               rvel.cd,
                               n=100,
                               scale='sqrt',
                               cell.colors=ac(cell.colors.cluster,alpha=0.5),
                               cex=0.8,
                               arrow.scale=3,
                               show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,
                               grid.n=40,
                               arrow.lwd=1,
                               do.par=F,
                               cell.border.alpha = 0.1)

# Visualize a fit for maker genes (we reuse rvel.cd to save on calcualtions here):

cat("<HR>")
cat("<H4>Velocity flow details for selected genes on UMAP map</H4>")
par( mfrow = c(2,2))
excluded_genes = vector()
for( gene_name in sort( c( "Ccr6", "Cd4", "Cxcr6", "Eomes", "Flt3", "Id2", "Il13", "Il22", "Il2rb", "Il33", "Inta4", "Intb1",
                     "Itgam", "Lta", "Plekhs1", "Rora", "Rorc", "Zbtb16", "Plppr3", "Gfra1", "Flt1", "Hic1"))){
  if( gene_name %in% rownames( emat) && gene_name %in% rownames( nmat)){
    gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells = 20, kGenes=1,
                                     fit.quantile=fit.quantile, cell.emb=emb_UMAP, cell.colors=cell.colors.cluster,
                                     cell.dist=cell.dist, old.fit=rvel.cd, do.par=FALSE,
                                     show.gene= gene_name)
  }else{
    cat("<HR><BR>The gene", gene_name, "is not present in the matrices<HR>")
  }
}
par( mfrow = c(1,1))
