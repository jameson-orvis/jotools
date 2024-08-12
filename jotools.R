####DUMPING ALL MY LITTLE TOOLS INTO A PACKAGE


library(skitools)
##library(GxG)
library(data.table)
library(Hmisc)
library(gUtils)
library(devtools)
library(skitools)
library(MASS)
library(EnhancedVolcano)
library(purrr)
##library(simplextree)
devtools::load_all('/gpfs/commons/home/jorvis/git/chromunity_experimental_branch') ##without zero weight edges
devtools::load_all('/gpfs/commons/home/jorvis/git/GxG_fork/GxG') ####MAKING THE CHANGE WHERE WE REMOVE dedupe in GXG
                                        #devtools::load_all('/gpfs/commons/home/jorvis/git/chromunity_new')
##devtools::load_all('/gpfs/commons/home/jorvis/git/chromunity')
##devtools::load_all('/gpfs/commons/home/jorvis/git/chromunity_fastversion')

setDTthreads(1)


#########throwing this function up here 
load_EP_features <- function(){
    ##########
    chain19to38 = rtracklayer::import.chain("~/DB/UCSC/hg19ToHg38.over.chain")
    eigen.porec = fread("~/projects/PoreC/db/f4dbd_GM12878_NlaIII_native_50kb_cooltools_dev_defaultnan.cis.vecs.tsv")
    setnames(eigen.porec, "chrom", "seqnames")
    eigen.porec[, bin := .I]
    gr.eigen = dt2gr(eigen.porec)
    k27 = import("/gpfs/commons/groups/imielinski_lab/data/ChIPseq/GM12878.ENCODE/GRCh38/H3K27ac/ENCFF218QBO.bigWig")
    track.k27 = gTrack(k27, y.field = 'score', bars = TRUE)
    peak = fread("/gpfs/commons/groups/imielinski_lab/data/ChIPseq/GM12878.ENCODE/GRCh38/H3K27ac/ENCFF211DND.bed")
    colnames(peak) <- c("chr","start","end","peak","val","A","B","C","D","E")
    gr.peak = dt2gr(peak)
    peak.signal.k27 = gr.peak %*% k27
    peak.eigen = peak.signal.k27 %*% gr.eigen
    peak.eigen = gr2dt(peak.eigen)
    pos.e1 = c("chr2", "chr6", "chr9", "chr13", "chr14", "chr16", "chr18", "chr19", "chr20", "chr22", "chrX")
    eigen.porec = na.omit(eigen.porec)
    eigen.porec[, E1 := ifelse(seqnames %in% pos.e1, E1*-1, E1)]
    eigen.porec[ , active := (E1 >=0) ]
    eigen.porec[ , in.active := (E1 <0)]
    tiles.en = gr.tile(hg_seqlengths(), 5e5)
    compartments = eigen.porec %>% dt2gr  %>% sort
########
    compartments = gr2dt(compartments)
    comp.a = compartments[E1 > 0] %>% dt2gr()
    comp.a = sortSeqlevels(comp.a)
    comp.a = sort(comp.a)
    comp.a = gr.reduce(comp.a)
    comp.a = gr2dt(comp.a)
    comp.a[, idx := seq_len(.N), by = seqnames]
    comp.a[, name := paste0("A", idx)]
    comp.b = compartments[E1 < 0] %>% dt2gr()
    comp.b = sortSeqlevels(comp.b)
    comp.b = sort(comp.b)
    comp.b = gr.reduce(comp.b)
    comp.b = gr2dt(comp.b)
    comp.b[, idx := seq_len(.N), by = seqnames]
    comp.b[, name := paste0("B", idx)]
    message("Comp done")
####
#####
    gm12878_exp = readRDS("/gpfs/commons/groups/imielinski_lab/DB/ENCODE/GM12878_gene_exp/gm12878_exp.rds")
    genes.gencode.all = import("/gpfs/commons/groups/imielinski_lab/DB/GENCODE/hg38/v32/gencode.v32.annotation.gtf")
    genes.gencode = genes.gencode.all %Q% (gene_type %in% "protein_coding") %Q% (type %in% "gene")
    gm12878_exp = merge(gm12878_exp, gr2dt(genes.gencode)[,.(seqnames, start, end, gene_name)], by.x = "SYMBOL", by.y = "gene_name", all.x = T)
    gm12878_exp = na.omit(gm12878_exp)
    gm12878_exp = gm12878_exp[, .(SYMBOL, TPM, FPKM, seqnames, start, end)]
    ## act and inact genes
    gm12878_inact = dt2gr(gm12878_exp[TPM == 0])

    gm12878_act = gm12878_exp[TPM > 0]
    gm12878_act[ , quantile := cut(TPM,
                                   breaks = quantile(TPM, probs = 0:4/4),
                                   labels = 1:4, right = FALSE)]
    gm12878_act = dt2gr(gm12878_act)
    message("Exp done")

    ## ctcf
    ctcf = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/files/GM12878_chip/GM12878_CTCF_hg38_ENCFF796WRU.bed.gz")
    setnames(ctcf, c("V1", "V2", "V3", "V7"), c("seqnames", "start", "end", "score"))
    ctcf = dt2gr(ctcf[order(score, decreasing = TRUE)])#[1:5000])

    ## h27ac
    h27ac = fread("/gpfs/commons/groups/imielinski_lab/data/ChIPseq/GM12878.ENCODE/GRCh38/H3K27ac/ENCFF211DND.bed")
    setnames(h27ac, c("V1", "V2", "V3", "V4"), c("seqnames", "start", "end", "peak"))
    h27ac = dt2gr(h27ac[order(V9, decreasing = TRUE)])#[1:5000])

    ## h4me3
    h4me3 = fread("/gpfs/commons/groups/imielinski_lab/data/ChIPseq/GM12878.ENCODE/GRCh38/H3K4me3/ENCFF190OOC.bed")
    setnames(h4me3, c("V1", "V2", "V3", "V4"), c("seqnames", "start", "end", "peak"))
    h4me3 = dt2gr(h4me3[order(V9, decreasing = TRUE)][1:5000])

    ## h4me1
    h4me1 = fread("/gpfs/commons/groups/imielinski_lab/data/ChIPseq/GM12878.ENCODE/GRCh38/H3K4me1/ENCFF006FIO.bed")
    setnames(h4me1, c("V1", "V2", "V3", "V4"), c("seqnames", "start", "end", "peak"))
    h4me1 = dt2gr(h4me1[order(V9, decreasing = TRUE)][1:5000])

    ## h36me3
    h36me3 = fread("/gpfs/commons/groups/imielinski_lab/data/ChIPseq/GM12878.ENCODE/GRCh38/H3K36me3/ENCFF986RLS.bed")
    setnames(h36me3, c("V1", "V2", "V3", "V4"), c("seqnames", "start", "end", "peak"))
    h36me3 = dt2gr(h36me3[order(V9, decreasing = TRUE)][1:5000])

    ## h27me3
    h27me3 = fread("/gpfs/commons/groups/imielinski_lab/data/ChIPseq/GM12878.ENCODE/GRCh38/H3K27me3/ENCFF490UWI.bed")
    setnames(h27me3, c("V1", "V2", "V3", "V4"), c("seqnames", "start", "end", "peak"))
    h27me3 = dt2gr(h27me3[order(V9, decreasing = TRUE)][1:5000])

    ## h9me3
    h9me3 = fread("/gpfs/commons/groups/imielinski_lab/data/ChIPseq/GM12878.ENCODE/GRCh38/H3K9me3/ENCFF833BJA.bed")
    setnames(h9me3, c("V1", "V2", "V3", "V4"), c("seqnames", "start", "end", "peak"))
    h9me3 = dt2gr(h9me3)##[order(V9, decreasing = TRUE)][1:5000])

    h9me3 = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/files/GM12878_chip/GSM733664_hg19_wgEncodeBroadHistoneGm12878H3k9me3StdPk.broadPeak.gz")
    setnames(h9me3, c("V1", "V2", "V3", "V4"), c("seqnames", "start", "end", "peak"))
    h9me3 = na.omit(h9me3) %>% dt2gr() %>% rtracklayer::liftOver(chain19to38)%>% grl.unlist %Q% (!duplicated(grl.ix))
    h9me3 = dt2gr(gr2dt(h9me3)[order(V7, decreasing = TRUE)][1:5000])

    ## Dnase
    dnase = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/files/GM12878_chip/GM12878_DNAseq.bed.gz")
    setnames(dnase, c("V1", "V2", "V3", "V7"), c("seqnames", "start", "end", "score"))
    dnase = dt2gr(dnase[order(score, decreasing = TRUE)][1:5000])

    ## chmm
    chmm.19 = import("/gpfs/commons/groups/imielinski_lab/DB/ENCODE/HMM/wgEncodeBroadHmmGm12878HMM.bed")
    chmm = chmm.19 %>% rtracklayer::liftOver(chain19to38)%>% grl.unlist %Q% (!duplicated(grl.ix))
    chmm.ep = dt2gr(gr2dt(chmm)[grepl("Promoter", name) | grepl("Enhancer", name)])
    chmm.e = dt2gr(gr2dt(chmm)[grepl("Strong_Enhancer", name)])

    ## SE
    sum.en = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/GM12878.super.enhancer.bed.txt") %>% setnames(., c("V1", "V2", "V3"), c("seqnames", "start", "end")) %>% dt2gr()
    sup.en.38 = grl.unlist(rtracklayer::liftOver(sum.en, chain19to38))

     ## Prom
    gm_prom_annot = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/GenoSTAN_promoters.bed.gz")
    gm_prom_annot = gm_prom_annot[grepl("E116", V4)] %>% setnames(., c("V1", "V2", "V3"), c("seqnames", "start", "end")) %>% dt2gr()

    ## Enh
    gm_enh_annot = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/GenoSTAN_enhancers.bed.gz")
    gm_enh_annot = gm_enh_annot[grepl("E116", V4)] %>% setnames(., c("V1", "V2", "V3"), c("seqnames", "start", "end")) %>% dt2gr()

    ## TFs
    p300 = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/files/GM12878_chip/TFBS/GM12878_ep300_narrowpeaks.bed.gz") %>% setnames(c("V1", "V2", "V3", "V4", "V5", "V7"), c("seqnames", "start", "end", "strand", "visual", "qscore")) ## %>% dt2gr()
    p300 = dt2gr(p300[order(qscore, decreasing = T)][1:5000])

    pol2 = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/files/GM12878_chip/TFBS/GM12878_pol2RA_narrowpeaks.bed.gz") %>% setnames(c("V1", "V2", "V3", "V4", "V5", "V7"), c("seqnames", "start", "end", "strand", "visual", "qscore")) ## %>% dt2gr()
    pol2 = dt2gr(pol2[order(qscore, decreasing = T)][1:5000])

    pol2s2 = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/files/GM12878_chip/TFBS/GM12878_pol2S2_narrowpeaks.bed.gz") %>% setnames(c("V1", "V2", "V3", "V4", "V5", "V7"), c("seqnames", "start", "end", "strand", "visual", "qscore")) ## %>% dt2gr()
    pol2s2 = dt2gr(pol2s2[order(qscore, decreasing = T)][1:5000])

    pol2s5 = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/files/GM12878_chip/TFBS/GM12878_pol2S5_narrowpeaks.bed.gz") %>% setnames(c("V1", "V2", "V3", "V4", "V5", "V7"), c("seqnames", "start", "end", "strand", "visual", "qscore")) ## %>% dt2gr()
    pol2s5 = dt2gr(pol2s5[order(qscore, decreasing = T)][1:5000])

    ## hcfc1 = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/files/GM12878_chip/GM12878_hcfc1.bed.gz")%>% setnames(c("V1", "V2", "V3", "V4", "V5", "V7"), c("seqnames", "start", "end", "strand", "visual", "qscore")) ## %>% dt2gr()
    ## hcfc1 = dt2gr(hcfc1)##[order(qscore, decreasing = T)][1:5000])

    ## chd2 = fread("/gpfs/commons/groups/imielinski_lab/projects/PoreC/files/GM12878_chip/GM12878_chd2.bed.gz")%>% setnames(c("V1", "V2", "V3", "V4", "V5", "V7"), c("seqnames", "start", "end", "strand", "visual", "qscore")) ##%>% dt2gr()
    ## chd2 = dt2gr(chd2[order(qscore, decreasing = T)][1:5000])

###############

#####
    mega.list = list(sup.en.38 = sup.en.38, gm_enh_annot = gm_enh_annot, gm_prom_annot = gm_prom_annot,
                     gm12878_act = gm12878_act, gm12878_inact = gm12878_inact,
                     h27ac = h27ac, h4me1 = h4me1, h4me3 = h4me3, h36me3 = h36me3, h27me3 = h27me3, h9me3 = h9me3,
                     dnase = dnase, 
                     pol2 = pol2, pol2s2 = pol2s2, pol2s5 = pol2s5, p300 = p300, ctcf = ctcf,
                     comp.a = comp.a, comp.b = comp.b)
    return(mega.list)
}

.make.volcano <- function(syn, facet.on = "annotation", or.thresh = 1, fdr.thresh= 0.1,
                          xl = c(-10, 10), yl = c(0, 25), x.ax = "log.est", y.ax = "fdr"){
  num.facets = nrow(unique(syn[, facet.on, with = F]))
  unique.facets = unique(syn[, facet.on, with = F])
  this.envol.list = list()
  for (i in 1:num.facets){
    this.facet = unique.facets[[1]][i]
    this.envol = EnhancedVolcano(syn[annotation == this.facet],
                                 pCutoff = fdr.thresh,
                                 FCcutoff = log2(1),
                                 subtitle = NULL,
                                 xlab = "Log2 Relative Risk",
                                 ylab = "-Log10 FDR",
                                 lab = NA,
                                 x = x.ax,
                                 y = y.ax,
                                 colAlpha = 0.5,
                                 legendLabSize = 7,
                                 legendLabels = c("NS", paste0("OR > ", or.thresh),
                                                  paste0("fdr < ", fdr.thresh),
                                                  paste0("fdr < ", fdr.thresh, " & OR > ", or.thresh)),
                                 drawConnectors = TRUE,
                                 widthConnectors = 0.05,
                                 boxedLabels = TRUE,
                                 title = this.facet,
                                 xlim = xl,
                                 ylim = yl,
                                 axisLabSize = 30,
                                 pointSize = 4
    )
    print(this.envol)
    this.envol.list[[this.facet]] <- this.envol
  }
  return(this.envol.list)
}


rebin_community = function(concatemers, this.chrom.w, resolution = 5e4, rebin_thresh=0.85) {
    tiles = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), resolution)
    this.chrom = gr2dt(concatemers %Q% (chid %in% this.chrom.w))
    this.list.chrom = pbmclapply(1:length(this.chrom.w), function(j){
        this.pr = dt2gr(this.chrom[chid %in% this.chrom.w[j]])
        sum.this.com = gr.sum((this.pr)+1e4)
        sum.this.com = gr2dt(sum.this.com)
        sum.this.com[, q := quantile(score, rebin_thresh), by = seqnames]
        sum.this.com[, q := ifelse(q < 5, 5, q)]
        active.cont = tryCatch((tiles %&% dt2gr(sum.this.com[score > q])), error = function(e) NULL)
        this.clust = gr2dt(gr.reduce(active.cont))
        this.clust[, chid := this.chrom.w[j]]
        return(this.clust)
    }, mc.cores  = 10)
    if(length(this.chrom.w) == 1){
        this.list.chrom = this.list.chrom$value
    }
    this.chrom.dt = rbindlist(this.list.chrom, fill = TRUE)
    return(this.chrom.dt)
}


load_bad_regions = function(chromosome) {
    this.chr = chromosome
    bands.td = gTrack::karyogram(file = "/gpfs/commons/groups/imielinski_lab/DB/UCSC/hg38.cytoband.txt")
    bands = bands.td@data
    bands = grl.unlist(do.call(`GRangesList`, bands))
    cen = bands %Q% (stain=="acen")
    if (this.chr != 'chrX') {
        chr.ind = as.numeric(sub("chr*","",this.chr))
    } else {
        chr.ind = 'X'
    }
    this.max = GRanges(paste0(this.chr, ":", hg_seqlengths()[chr.ind]-1e6, "-",  hg_seqlengths()[chr.ind]))         
    this.min = GRanges(paste0(this.chr, ":", "1-1e6"))                                                                                                                       
    this.cen = (cen %Q% (seqnames == this.chr))+1e6
    this.bad = c(this.min, this.cen, this.max) 
    return(this.bad)
}
    


####specifically for running giga chromunity a little bit
evaluate_synergy = function(res, leave_out_concatemers, chid.to.test, chromosome = NULL, filter_binsets = TRUE, folder = NULL, rebin_thresh = 0.85, mc.cores = 20, numchunks = mc.cores*200 + 1, dont_rebin=FALSE, sub.binset.order=3, resolution=1e4) {


    if (!dir.exists(folder)) {
        stop("output folder does not exist")
    }
    if(is.null(chromosome)){
        chromosome = c(paste0("chr", c(as.character(1:22), "X")))
        all.bad = pbmclapply(chromosome, mc.preschedule=FALSE, function(chr) {
            bad.gr = muffle(load_bad_regions(chr))
            bad.dt = gr2dt(bad.gr)
            return(bad.dt)
        })
        this.bad = rbindlist(all.bad) %>% dt2gr
    } else {
        this.bad = load_bad_regions(chromosome)
    }
    tiles = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), resolution)
    gc5b = readRDS("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/gc.38.rds")
    ## create a list of covariates
    cov_list = list(gc5b)
    ## Specify what kind of covariate it is. Score will be aggregated over bins while the number of intervals will be calculated otherwise.
    names(cov_list) <- c("score:gc.cov")
    ## Make the covariate object
    gc_cov = covariate(name = c("gc"), type = c("numeric"), field = c("score"), data = cov_list)
    gc.cov = gc_cov

    ##handle all that rebinning shit UPSTREAM please thank you 
    
    ## this.dat.chrom = data.table()
    ## if(dont_rebin) this.chrom = gr2dt(concatemers)
    ## if (filter_binsets) {
    ##     this.chrom = this.chrom[support > summary(unique(this.chrom[, .(support, chid)])$support)[3]]####filtering out all bins below median support
    ## }
    ## this.chrom.w = unique(chid.to.test)
    ## if(dont_rebin) {
    ##     print('yay')
    ##     this.chrom.dt = this.chrom
    ## } else {
    ##     this.chrom.dt = rebin_community(res$concatemers, this.chrom.w, resolution=resolution)
    ## }

    ####

    ##browser()
    this.chrom.dt = gr2dt(res$binsets)
    
    this.chrom.dt[, cardinality := .N, by = chid]
    this.chrom.dt = na.omit(this.chrom.dt)
    ##
    this.all.dat = copy(this.chrom.dt)
    this.all.dat = this.all.dat[cardinality > 2]
    this.all.dat[, bid := chid]
    this.all.dat = this.all.dat[seqnames %in% chromosome]
    this.all.dat[, overall.cardinality := cardinality, by = bid]
    this.chrom.card = unique(this.all.dat[, .(overall.cardinality, bid)])
    ######HELLL NAH
    this.all.dat = this.all.dat[cardinality < 100]
    ####
    this.sub.parq = leave_out_concatemers
    this.sub.parq$cid = this.sub.parq$read_idx
    this.all.dat[, binid := .I]
    ##
    ###filtering out the bad regions
    ###this makes sense why it wouldn't be here for RE chromunity cause you're only looking at annotated regions anyway

###remove only bins which intersect bad region, not all

    this.all.dat$bid = this.all.dat$chid
    this.bad.chrom = unique(as.character((dt2gr(this.all.dat) %&% (this.bad))$bid))                                                 
    print(this.bad.chrom)
    this.all.dat = this.all.dat[!bid %in% this.bad.chrom]  
                                        #browser()
    #debug(annotate)
    chrom.annot.output = annotate(binsets = dt2gr(this.all.dat[, bid := chid]),
                              k = sub.binset.order,
                              concatemers = this.sub.parq,
                              covariates = gc.cov, resolution = resolution,
                              mc.cores = mc.cores, numchunks = numchunks)

    this.chrom.dat = chrom.annot.output[[1]]
    this.chrom.dat.2 = this.chrom.dat
    this.chrom.dat = merge(this.chrom.dat, this.chrom.card[, bid := as.factor(bid)], by = "bid")
    this.chrom.dat[, annotation := "chromunity"]
    set.seed(198)
    message("generating background binsets for the model")
    ##
                                        #browser()

    n.chrom = length(unique(this.chrom.dat$bid))
    if(length(chromosome) > 1){
        back.dt = re_background(binsets = dt2gr(this.all.dat), n = (n.chrom*3), resolution = resolution)#, mc.cores=mc.cores)
    } else {

        back.dt = sliding_window_background(chromosome = chromosome, binsets = dt2gr(this.all.dat), n = (n.chrom*3), resolution = resolution)#, mc.cores=mc.cores)
    }
    
    upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T)
    setkeyv(back.dt, c("seqnames", "start"))
    back.dt[, V1 := NULL]
    back.dt = na.omit(back.dt)
    back.dt = back.dt[!bid %in% back.dt[width < (resolution-1)]$bid]
    back.dt = gr2dt(gr.reduce(dt2gr(back.dt), by = "bid"))
    back.dt$bid <- as.factor(back.dt$bid)
    back.dt = merge(back.dt, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
    back.dt = back.dt[end < V2][start < V2]
    back.dt[, overall.cardinality := .N, by = bid]
    back.dt = back.dt[overall.cardinality > 1]
    ##
    this.card = unique(back.dt[, .(bid, overall.cardinality)])
    back_gr = dt2gr(back.dt)
    message("extracting background binsets distances")



    back.train.output = annotate(binsets = dt2gr(back.dt),
                                   interchromosomal.table = NULL, #all.hr.dt.mean,
                                   gg = NULL,
                                   concatemers = this.sub.parq, k = sub.binset.order,
                                   covariates = gc.cov, resolution = resolution,
                                   mc.cores = mc.cores, numchunks = numchunks)
    this.back.train.dat = back.train.output[[1]]
    this.back.train.dat = merge(this.back.train.dat, this.card[, bid := as.factor(bid)], by = "bid")
    this.back.train.dat[, annotation := "random"]
    this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][width <= resolution]$bid)]
    this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][min.dist < resolution]$bid)]
    back.dt[, binid := .I]
    this.bad.train = unique(as.character((back_gr %&% (this.bad))$bid))
    this.back.train.dat = this.back.train.dat[!bid %in% this.bad.train]
    this.back.train.dat = this.back.train.dat[!bid %in% this.back.train.dat[, .(sum(count)), by = bid][V1 == 0]$bid]
####    

    ##browser()
    
    back.model = fit(na.omit(this.back.train.dat)[sum.counts > 0][, setdiff(names(this.back.train.dat), c('overall.cardinality', 'chr', 'annotation')), with = F])


    
    this.chrom.dat = sscore(this.chrom.dat, model = back.model)
    message("generating random binsets for testing")
    n.chrom = length(unique(this.chrom.dat$bid))
    this.all.dat = this.all.dat[seqnames %in% chromosome]

    if(length(chromosome) > 1){
        back.test = re_background(binsets = dt2gr(this.all.dat),
                                          n = n.chrom,
                                          resolution = resolution)
     } else {
        back.test = sliding_window_background(binsets = dt2gr(this.all.dat), chromosome = chromosome,
                                            n = n.chrom,
                                            resolution = resolution)
    }                      

    back.test[, V1 := NULL]
    back.test = na.omit(back.test)
    setkeyv(back.test, c("seqnames", "start"))
    back.test[, start := ifelse(start < 0, 0, start)]
    back.test = back.test[!bid %in% back.test[width < (resolution-1)]$bid]
    back.test = gr2dt(gr.reduce(dt2gr(back.test), by = "bid"))
    back.test = merge(back.test, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
    back.test = back.test[end < V2][start < V2]
    back.test[, overall.cardinality := .N, by = bid]
    back.test = back.test[overall.cardinality > 1]
    back_test = dt2gr(back.test)
###
    back.test$bid <- as.factor(back.test$bid)
    this.back.card = unique(back.test[, .(bid, overall.cardinality)])
    bid.back = unique(back.test$bid)
    message("extracting random binsets distances") 
    back.test.output = annotate(binsets = dt2gr(back.test), k=sub.binset.order,
                                  concatemers = this.sub.parq,
                                  covariates = gc.cov, resolution = resolution,
                                  mc.cores = mc.cores, numchunks=numchunks)
    this.back.test.dat = back.test.output[[1]]
    this.back.test.dat = merge(this.back.test.dat, this.back.card, by = "bid")
    this.back.test.dat[, annotation := "random"]
    this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][width <= resolution]$bid)]
###
    this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][min.dist < resolution]$bid)]
####
####
    this.bad.dat = unique(as.character((back_test %&% (this.bad))$bid))
    this.back.test.dat = this.back.test.dat[!bid %in% this.bad.dat]
    back_test = gr2dt(back_test)[bid %in% unique(this.back.test.dat$bid)]
    #this.back.test.dat = sscore(this.back.test.dat, model = back.model)
####
    back.test[, binid := .I]
#####
    set.seed(178)
    all.bid = unique(this.all.dat$bid)
###
    ##browser()
    message("starting random walks")
    message("extracting shuffled binsets distances")

    ##browser()
    ##nr = 1
    
    ##resolution=1e4

    ##all.bid
    ##
    ##browser()
    s.chrom = synergy(binsets = dt2gr(this.all.dat), #theta = back.model$model$theta,
                      annotated.binsets = this.chrom.dat, model = back.model)
    s.chrom$fdr = signif(p.adjust(s.chrom$p, "BH"), 2)
    print('the hit rate!!!')
    print(s.chrom[, .N])
    print(s.chrom[fdr<0.1, .N])
    

    this.all.dat.shuff5 = pbmclapply(1:length(all.bid), function(nr){
         this.clust = dt2gr(this.all.dat[bid == all.bid[nr]])
         this.chrs = .chr2str(as.character(unique(seqnames(this.clust))))
         this.clust.wind = gr.reduce(this.clust+2e6)
         upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T)
         this.clust.wind = (gr2dt(this.clust.wind)[, start := ifelse(start < 0, 1, start)])
         this.clust.wind = merge(this.clust.wind, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
         this.clust.wind[, end := ifelse(end > V2, V2, end)]
         this.clust.wind = dt2gr(this.clust.wind)

         ##this.clust.wind
         this.sub.win = gr2dt(this.sub.parq %&% this.clust.wind)
         ##this.sub.win[, new.count := .N, by = read_idx]

         ###throw this in there
         

         
         this.tiles.orig = gr.tile(this.clust.wind, resolution)
         
         this.tiles.orig$binid = 1:length(this.tiles.orig)
         this.sub.win = bin_concatemers(dt2gr(this.sub.win), this.tiles.orig)
         this.sub.win = unique(this.sub.win, by=c('read_idx','binid'))
         this.sub.win[, new.count := .N, by = read_idx]
         ##
         card = unique(this.sub.win[, .(read_idx, new.count)])
         this.steps = sum(card$new.count)

         if (nrow(this.sub.win) > 0){
             this.tgm = cocount(dt2gr(this.sub.win), bins = (this.tiles.orig), by = 'read_idx', full = T)
             ##debug(shuffle_concatemers)
             ##shuff.concats = shuffle_concatemers(this.sub.win, this.tgm)
             ##shuff.concats$cid = shuff.concats$read_idx

             ##ppdf(plot(c(shuff.contacts$gtrack(name='shuff', clim=c(0,100)), this.tgm$gtrack(name='normal',clim=c(0,100))), gr.reduce(this.tiles.orig) + 3e5))
             ##ppdf(plot(this.tgm$gtrack(clim=c(0,100)), gr.reduce(this.tiles.orig) + 3e5))
             ##shuff.contacts = cocount(dt2gr(shuff.concats), bins = this.tiles.orig, by='read_idx', full=T)
             A = this.tgm$mat %>% as.matrix
             rownames(A) <- NULL
             colnames(A) <- NULL
             A[cbind(1:nrow(A), 1:nrow(A))] = 0
             A = A + t(A)
             An = A 
             An = round(1+10*An/min(An[An>0]))
             edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
             ## ##
             G = graph.edgelist(edges[, cbind(row, col)])
             RW = random_walk(G, start = 1, steps = sum(card$new.count)) %>% as.numeric
             ## ##
             rm(G)
             rm(edges)
             gc()
             out = this.tgm$gr[RW]%>% gr2dt()
             out$read_idx = card[, rep(read_idx, new.count)] 
             out[, bid := all.bid[nr]]
             out[, cid := read_idx]

             ##debug(annotate)
             sh.output = annotate(binsets = this.clust, verbose = F,
                                          k = sub.binset.order,
                                          concatemers = dt2gr(out),
                                          covariates = gc.cov, resolution = resolution, 
                                          mc.cores = 5, numchunks = numchunks)
             ## shuff.concats$cid = shuff.concats$read_idx
             ## sh.output.2 = annotate(binsets = this.clust, verbose = F,
             ##                              k = 5,
             ##                              concatemers = dt2gr(shuff.concats),
             ##                              covariates = gc.cov, resolution = resolution, 
             ##                              mc.cores = 10, numchunks = numchunks)
             
             this.chrom.sh.dat = sh.output[[1]]
         } else {
             this.chrom.sh.dat = data.table(NA)
         }
         return(this.chrom.sh.dat)
    }, mc.cores = 5, mc.preschedule = T)
    #this.shuff.chrom = tryCatch(rbindlist(the.shuff.list, fill = T), error = function(e) NULL)


    ##browser()

    this.all.dat.shuff5.ne = this.all.dat.shuff5[sapply(this.all.dat.shuff5, function(x) !inherits(x, "try-error"))]
    this.all.sh = rbindlist(this.all.dat.shuff5.ne, fill =T)
###
    this.all.sh = sscore(this.all.sh, model = back.model)
    sh.chrom = tryCatch(synergy(binsets = dt2gr(this.all.dat), annotated.binsets = this.all.sh, model = back.model), error = function(e) NULL)
    sh.chrom$fdr = signif(p.adjust(sh.chrom$p, "BH"), 2)
    ##
#####
    this.back.test.dat = sscore(this.back.test.dat, model = back.model)
    theta = back.model$model$theta
    s.chrom = synergy(binsets = dt2gr(this.all.dat), #theta = back.model$model$theta,
                      annotated.binsets = this.chrom.dat, model = back.model)
    s.chrom$fdr = signif(p.adjust(s.chrom$p, "BH"), 2)
    b.chrom = synergy(binsets = dt2gr(back.test), #theta = back.model$model$theta,
                      annotated.binsets = na.omit(this.back.test.dat), model = back.model)
    b.chrom$fdr = signif(p.adjust(b.chrom$p, "BH"), 2)
    ##
    synergy.inter.EP = rbind(s.chrom[, annotation := "chromunity"],
                             b.chrom[, annotation := "random"], 
                             sh.chrom[, annotation := "shuffled"], fill = T)
    if (!is.null(folder)) {
        saveRDS(this.chrom.dat, paste0(folder,'chrom_annotate.rds'))
        saveRDS(this.back.test.dat, paste0(folder,'back_annotate.rds'))
        saveRDS(this.all.sh, paste0(folder,'shuffled_annotate.rds'))
        saveRDS(this.all.dat, paste0(folder,'binsets.rds'))
        saveRDS(back.model, paste0(folder,'back_model.rds'))
        saveRDS(synergy.inter.EP, paste0(folder,'synergy_results.rds'))
    }
    return(synergy.inter.EP)
}

make_bin_graph = function(concatemers, bins, hyperedge.thresh = 3, genes, methylation = NULL, h27ac = NULL, features=NULL) {
                                        #bins = bins %*% GRanges(view_range)
    max.slice = 1e6
    mc.cores=5
    verbose=TRUE
    concatemers$binid = gr.match(concatemers, bins, max.slice = max.slice, mc.cores =  mc.cores, verbose = verbose)

    ## maybe NA need to be removed
    concatemers = concatemers %Q% (!is.na(binid))


                                        #REALLY INTERESTED IN MYC
    reads = as.data.table(concatemers)[, `:=`(count, .N), by = cid]    
    reads2 = reads[count >= hyperedge.thresh, ]
    reads2$cid = factor(reads2$cid)
    reads2[, cidi := as.integer(cid)]
    
    unique_reads2 = unique(reads2, by = c('binid', 'cidi'))
    numbins = unique_reads2[, length(binid), by='cidi']
    unique_reads2_trim = merge(unique_reads2, numbins, by='cidi')[V1 >= hyperedge.thresh,]   ###only keep hyperedges that hit at least 3 bins
    rm(unique_reads2)

    ##trimming unnecessary fields
    unique_reads2_trim[, setdiff(colnames(unique_reads2_trim), c('cidi','binid')) := NULL]

    ##unique_binid = unique(unique_reads2_trim$binid)
    ##num.bins = length(unique_binid)
    ##concats.per.bin = sum(unique_reads2_trim[, .N, by='binid']$N) / num.bins
    ##concats.per.bin = sum(bincounts$N) / dim(bincounts)[1]
    ##total.pairs = num.concats * bins.per.concat * concats.per.bin
    ##numchunks = ceiling(total.pairs / 1e8)   ####~100 million pairs per chunk
    ##if (numchunks < 10) numchunks = 10

    numchunks = 10
    unique_binids = unique(unique_reads2_trim$binid)
    ubidl = split(unique_binids, ceiling(runif(length(unique_binids))*numchunks)) ## randomly chop up cids
    setkey(unique_reads2_trim, 'binid')
    
    #dt = gr2dt(concatemers[,c('read_idx','binid')])
    #dt[, numbins := .N, by=read_idx]
    #dt = dt[numbins > 2]

    #concatm = sparseMatrix(unique_reads2_trim$cidi %>% as.integer, unique_reads2_trim$binid %>% as.integer, x=1, repr='R')
    edgelist = pbmclapply(ubidl, mc.cores=mc.cores, mc.preschedule=TRUE, function(ubids) {    
        dt = unique_reads2_trim[.(ubids), .(binid1 = binid, cidi = cidi), by=binid]
        dt[, binid := NULL]
        dt.merged = merge(dt, unique_reads2_trim, by='cidi', how='left', allow.cartesian = TRUE)[binid1 < binid]
        edgelist = dt.merged[, .N, by=c('binid1','binid')]
        edgelist = edgelist[N >= hyperedge.thresh]  ###only bin pairs if there are 3 concatemers supporting it
        #colnames(edgelist)[colnames(edgelist) == "cidi1"] = "bx1"
        #colnames(edgelist)[colnames(edgelist) == "cidi"] = "bx2"
        #colnames(edgelist)[colnames(edgelist) == "N"] = "mat"
        #p1 = concatm[edgelist[, bx1],]
        #p2 = concatm[edgelist[, bx2],]
        #total = Matrix::rowSums(p1 | p2)
        #edgelist$tot = total
        #edgelist[, frac := mat/tot]
        #edgelist = edgelist[mat > 2 & tot > 2]
        return(edgelist)
    })
    edgelist = rbindlist(edgelist)

    #bincombs = dt[, t(combn(binid, 2)) %>% as.data.table, by = read_idx]
    #edgelist = bincombs
    #edgelist = bincombs[, .N, by=c('V1','V2')]

    edgelist_2 = copy(edgelist)
    edgelist_2$binid = edgelist$binid1
    edgelist_2$binid1 = edgelist$binid
    edgelist_3 = rbind(edgelist, edgelist_2)
    final.edgelist = edgelist_3[, N := sum(N), by=c('binid1','binid')]

    final.edgelist = final.edgelist[N >= hyperedge.thresh]
    A = sparseMatrix(final.edgelist$binid1, final.edgelist$binid, x=final.edgelist$N)

    G_bins = graph.adjacency(A, weighted = TRUE, mode = "undirected")

    ##genes = genes.gencode.all[, c('gene_name')]
    ##G_bins = set_vertex_attr(G_bins, 'genomic_pos', V(G_bins), (V(G_bins) %>% as.integer) / max(V(G_bins) %>% as.integer))

    bins$binid = 1:length(bins)
    annotated.targets = bins %$% genes

    if(!is.null(features)){
        compartment.a = annotated.targets %O% dt2gr(features$comp.a)[, 'active']
        compartment.b = annotated.targets %O% dt2gr(features$comp.b)[, 'active']
        mask.a = compartment.a > 0
        mask.b = compartment.b > 0
        G_bins = set_vertex_attr(G_bins, 'active', value=0)
        G_bins = set_vertex_attr(G_bins, 'active', annotated.targets$binid[mask.a], 1)
        G_bins = set_vertex_attr(G_bins, 'active', annotated.targets$binid[mask.b], -1)

        ac = annotated.targets %$% features$h27ac[, 'V9']
        G_bins = set_vertex_attr(G_bins, 'h27ac', value=0)
        mask = !is.na(ac$V9)
        G_bins = set_vertex_attr(G_bins, 'h27ac', ac$binid[mask], ac$V9[mask])
    }

    annotated.targets$genomic_pos = annotated.targets$binid / max(annotated.targets$binid)
    annotated.targets.dt = gr2dt(annotated.targets)
    coords = annotated.targets.dt[, .(binid, paste0(seqnames, ':', start, '-', end))]
    G_bins = set_vertex_attr(G_bins, 'gene_coords', coords$binid, coords$V2)
    
    G_bins = set_vertex_attr(G_bins, 'gene', annotated.targets$binid, annotated.targets$gene_name)
    G_bins = set_vertex_attr(G_bins, 'seqnames', annotated.targets$binid, gr2dt(annotated.targets)$seqnames)
    G_bins = set_vertex_attr(G_bins, 'genomic_pos', annotated.targets$binid, annotated.targets$genomic_pos)

    if (!is.null(methylation)) {
        methylation = methylation[,'meth']
        methyl = bins %$% methylation
        methyl = gr2dt(methyl)
        methyl[is.na(meth), meth := 0]
        methyl$meth = as.double(methyl$meth) + 0.00000001
        G_bins = set_vertex_attr(G_bins, 'methylation', methyl$binid, methyl$meth)
    }
    #browser()
    if (!is.null(h27ac)) {
        acetyl = bins %$% h27ac
        acetyl = gr2dt(acetyl)
        acetyl[is.na(V9), V9 := 0]
        acetyl$V9 = as.double(acetyl$V9) + 0.0000001
        G_bins = set_vertex_attr(G_bins, 'acetylation', methyl$binid, acetyl$V9)
    }
    return(G_bins)
}

bin_concatemers = function(concatemers, bins, max.slice = 1e6, mc.cores=5, verbose=TRUE, hyperedge.thresh=NULL) {
    max.slice = 1e6
    mc.cores=5
    verbose=TRUE
    concatemers$binid = gr.match(concatemers, bins, max.slice = max.slice, mc.cores =  mc.cores, verbose = verbose)

    concatemers$cid = concatemers$read_idx
    ## maybe NA need to be removed

    concatemers = concatemers %Q% (!is.na(binid))

                                        #REALLY INTERESTED IN MYC
    reads = as.data.table(concatemers)[, `:=`(count, .N), by = cid]    
    reads[, cidi := as.integer(cid)]
    return(reads)
}

make_line_graph = function(concatemers, bins, hyperedge.thresh = 3, prebinned = FALSE, mc.cores=5, order=NULL) {
    
    if(!prebinned){
        reads2 = bin_concatemers(concatemers, bins, max.slice = 1e6, mc.cores=5, verbose=TRUE, hyperedge.thresh = hyperedge.thresh)
    } else {
        reads2 = concatemers
    }
    ##browser()
    unique_reads2 = unique(reads2, by = c('binid', 'cidi')) ##DEDUPINGTON
    numbins = unique_reads2[, length(binid), by='cidi']
    unique_reads2_trim = merge(unique_reads2, numbins, by='cidi')[V1 >= hyperedge.thresh,]   ###only keep hyperedges that hit at least 3 bins
    rm(unique_reads2)

    ##trimming unnecessary fields
    unique_reads2_trim[, setdiff(colnames(unique_reads2_trim), c('cidi','binid')) := NULL]

    ##unique_binid = unique(unique_reads2_trim$binid)
    ##num.bins = length(unique_binid)
    ##concats.per.bin = sum(unique_reads2_trim[, .N, by='binid']$N) / num.bins
    ##concats.per.bin = sum(bincounts$N) / dim(bincounts)[1]
    ##total.pairs = num.concats * bins.per.concat * concats.per.bin
    ##numchunks = ceiling(total.pairs / 1e8)   ####~100 million pairs per chunk
    ##if (numchunks < 10) numchunks = 10
    bincounts = unique_reads2_trim[, .N, by=binid] #number of concatemers that hits every bin
    bincounts = bincounts[N > 1]

    

    unique_cidi = unique(unique_reads2_trim$cidi)
        
    concat.counts = unique_reads2_trim[, .N, by='cidi']
    unique_reads2_trim_merge = copy(unique_reads2_trim)
    if(!is.null(order)){
        concats.to.keep = concat.counts[N==order]$cidi
        unique_cidi = unique_cidi[unique_cidi %in% concats.to.keep]
        sub.concats = concat.counts[N==(order-1)]$cidi
        unique_reads2_trim_merge = unique_reads2_trim[cidi %in% sub.concats]
    }    
    unique_binid = unique(unique_reads2_trim$binid)
    concats.per.bin = mean(unique_reads2_trim[, .N, by='binid']$N)
    bins.per.concat = mean(unique_reads2_trim[, .N, by='cidi']$N)
    total.pairs = length(unique_cidi) * bins.per.concat * concats.per.bin
    numchunks = ceiling(total.pairs / 1e8)   ####~100 million pairs per chunk
    if (numchunks < 10) numchunks = 10
    print(numchunks)

    ucidl = split(unique_cidi, ceiling(runif(length(unique_cidi))*numchunks)) ## randomly chop up cids
    setkey(unique_reads2_trim, 'cidi')
    #numchunks = 10
    
    
    
    #dt = gr2dt(concatemers[,c('read_idx','binid')])
    #dt[, numbins := .N, by=read_idx]
                                        #dt = dt[numbins > 2]
    ##browser()
    concatm = sparseMatrix(unique_reads2_trim$cidi %>% as.integer, unique_reads2_trim$binid %>% as.integer, x=1, repr='R')

    ##browser()
    ##cidis = ucidl[[1]]
    edgelist = pbmclapply(ucidl, mc.cores=mc.cores, mc.preschedule=TRUE, function(cidis) {
        dt = unique_reads2_trim[.(cidis), .(cidi1 = cidi, binid = binid), by=cidi]
        dt[, cidi := NULL]
        if(!is.null(order)){
            dt.merged = merge(dt, unique_reads2_trim_merge, by='binid', how='left', allow.cartesian = TRUE)
        } else {
            dt.merged = merge(dt, unique_reads2_trim_merge, by='binid', how='left', allow.cartesian = TRUE)[cidi1 < cidi]
        }
        edgelist = dt.merged[, .N, by=c('cidi1','cidi')]
        edgelist = edgelist[N >= hyperedge.thresh]
        colnames(edgelist)[colnames(edgelist) == "cidi1"] = "bx1"
        colnames(edgelist)[colnames(edgelist) == "cidi"] = "bx2"
        colnames(edgelist)[colnames(edgelist) == "N"] = "mat"
        p1 = concatm[edgelist[, bx1],]
        p2 = concatm[edgelist[, bx2],]
        total = Matrix::rowSums(p1 | p2)
        edgelist$tot = total
        edgelist[, frac := mat/tot]
        edgelist = edgelist[mat >= hyperedge.thresh & tot >= hyperedge.thresh]
        return(edgelist)
    })

    edgelist = rbindlist(edgelist)
    print(edgelist)

    p1 = concatm[edgelist[, bx1],]
    p2 = concatm[edgelist[, bx2],]
    mat = p1 & p2
    ##browser()
    df.t = summary(t(mat)) %>% as.data.table
    list.agg = df.t[, .(list(i)), by=j]

                                        #saveRDS(edgelist, 'honker_run/edgelist.rds')
    #browser()
    edgelist_2 = copy(edgelist)
    edgelist_2$bx2 = edgelist$bx1
    edgelist_2$bx1 = edgelist$bx2
    edgelist_3 = rbind(edgelist, edgelist_2)
    
    knn.dt = edgelist_3

    ###This will go ahead and factor the strings that represent each unique 3 way overlap 
    str.identifiers = list.agg$V1 %>% paste0
    str.identifiers.factor = factor(str.identifiers)
    list.agg$hashes = str.identifiers.factor
    
    return(list(knn.dt, reads2, list.agg))
}


annotate_graph = function(G, concatemers, gene_range = NULL) {
    seqcount = gr2dt(concatemers)[, .N, by=c('cid', 'seqnames')]
    seqcount[, is_multichromosomal := length(seqnames) > 1, by = cid]
    seqcount = unique(seqcount[order(N, decreasing=T)], by='cid')
    seqcount = seqcount[cid %in% as.integer(igraph::V(G)$cid),]
    setkey(seqcount, 'cid')
    G = igraph::set_vertex_attr(G, 'seqnames', value=as.character(seqcount[cid %in% concatemers$cid]$seqnames))
    if(!is.null(gene_range)) {
        concats.gr = concatemers
        gene.concats = concats.gr %*% gene_range
        gene.mask = V(G)$cid %in% gene.concats$cid
        G = igraph::set_vertex_attr(G, 'hits_gene', value=as.integer(gene.mask))
    }
    return(G)
}

plot_gene = function(res, chid.to.use, gene_range, gene_name, view_range=NULL, plot_filename, rebin_thresh=NULL, features) {
    if(!is.null(rebin_thresh)){
        binsets = rebin_community(res, chid.to.use, rebin_thresh=0.85)
    } else {
        binsets = gr2dt(res$binsets %Q% (chid == chid.to.use))
    }

    print(binsets)
    if(is.null(view_range)) {
                                        #view_range = gene_range + 1e5
        view_range = reduce(dt2gr(binsets)) + 1e6
    }
    #browser()
    dt = gr2dt(view_range)
    this.chrom = gr2dt(res$concatemers %Q% (chid==chid.to.use))
    this.pr = dt2gr(this.chrom)
    pr = GenomicRanges::split((this.pr), this.pr$cid)
    sum.this.com = gr.sum((this.pr)+1e3)

    tiles.tgm = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), 2.5e4)
    this.tgm = cocount(this.pr, bins = unique(tiles.tgm %&% (view_range)), by = 'cid')

    ppdf(plot(c(gTrack(dt2gr(binsets), name='Binset', height=5), gTrack(gene_range, name=gene_name, height=5), gTrack(features$sup.en.38, height=1, name='Superenhancer'), gTrack(features$gm_prom_annot, height=1, name='Promoter'), gTrack(features$gm_enh_annot, height=1, name='Enhancer'), gTrack(features$ctcf, y1 = 100, y.field = "score", bars = T, height = 2, name='CTCF peaks'), gTrack(features$h27ac, y1 = 100, y.field = "V9", bars = T, height = 2, name='H3K27ac peaks'), gTrack(pr %>% unname, name='concatemers'), this.tgm$gtrack(cmap.max = 120), gTrack(sum.this.com, height = 2,  y.field = "score", bars = T)), view_range), height=40, cex=0.7, plot_filename)
}

make_chromunities_from_line_graph = function(edgelist, reads2, bins) {
    reads2$binid = factor(reads2$binid)
    setkey(reads2, binid)
    ubid = reads2[, .N, by = binid][N>1 , binid]
 
    
    reads2 = reads2[.(ubid), ]
    ## added for subsamp
    reads2$binid = factor(reads2$binid) 

    reads2$cid = factor(reads2$cid)
    ucid = levels(reads2$cid)
    edgelist_3 = edgelist
    #edgelist_2 = copy(edgelist)
    #edgelist_2$bx2 = edgelist$bx1
    #edgelist_2$bx1 = edgelist$bx2
    #edgelist_3 = rbind(edgelist, edgelist_2)
    edgelist_3$nmat = edgelist_3$mat
    edgelist_3$nfrac = edgelist_3$frac
    #setkeyv(edgelist_3, c("nfrac", "nmat"))
    edgelist_3 = unique(edgelist_3)
    #edgelist_3 = edgelist_3[order(nfrac, nmat, -bx2, decreasing = T)]

    #fuck knn all my homies hate knn
    #k=k.knn
    #knn.dt = edgelist_3[mat > 2 & tot > 2, .(knn = bx2[1:k]), by = bx1][!is.na(knn), ]

    A = sparseMatrix(edgelist_3$bx1, edgelist_3$bx2, x = 1, dims = c(max(edgelist_3), max(edgelist_3)))

    #knn.shared = knn %*% knn
    
    #KMIN = k.min
    #A = knn.shared * sign(knn.shared > KMIN)
    A[cbind(1:nrow(A), 1:nrow(A))] = 0
                                        #  A <- as(A, "sparseMatrix")
    A = A + t(A)
    A[1,1] = 0
    G = graph.adjacency(A, weighted = TRUE, mode = "undirected")

    #print(sum(E(G)$weight == 0))
    #cidis = as.integer(colnames(A))
    cidis = as.integer(igraph::V(G))
    #seqcount = reads2[cidi %in% as.integer(igraph::V(G)),]
    #seqcount = seqcount[, .N, by=c('cidi', 'seqnames')]
    #seqcount[, is_multichromosomal := length(seqnames) > 1, by = cidi]
    #seqcount = unique(seqcount[order(N, decreasing=T)], by='cidi')
    G = igraph::set_vertex_attr(G, 'cidi', value=cidis)
    #setkey(seqcount, 'cidi')
    
    #G = igraph::set_vertex_attr(G, 'seqnames', value=as.character(seqcount[cidis]$seqnames))
    #G = igraph::set_vertex_attr(G, 'is_multichromosomal', value=as.integer(seqcount[cidis]$is_multichromosomal))

                                        #cl.l = cluster_fast_greedy(G)
    isolated = which(degree(G)==0)
    G = delete.vertices(G, isolated)

    cl.l = cluster_leiden(G, objective_function='modularity')


    
    cl = cl.l$membership
    ## rename so chid has most support
    cls = 1:max(cl)
    names(cls) = cl %>% table %>% sort %>% rev %>% names
    cl = cls[as.character(cl)]

    cidi.cid = setkey(unique(reads2[, c('cid','cidi')], by=c('cid','cidi')), 'cidi')
    dt = data.table(nodes=V(G)$cidi) 
    cid.dt = cidi.cid[dt$nodes]
    memb.dt = data.table(cid=as.numeric(as.character(cid.dt$cid)), chid = cl, cidi=cid.dt$cidi, V=1:dim(cid.dt)[1])
    
    G = igraph::set_vertex_attr(G, 'chid', index=memb.dt$V, value=memb.dt$chid)
    G = igraph::set_vertex_attr(G, 'cid', index=memb.dt$V, value=memb.dt$cid)
    #G_asdf = igraph::set_vertex_attr(G, 'cid', index=231098, value=9040456)
    #if (!is.null(filename))
    #{
    #    cmessage('writing network file')
    #    write_graph(G, filename, format='gml')
    #}
    browser()
    reads2[, cid := as.numeric(as.character(cid))]
    reads2 = merge(reads2, memb.dt[, c('cid','chid')], by = "cid")
    reads2[, `:=`(support, length(unique(cid))), by = chid]
    reads2 = dt2gr(reads2)

    cc = reads2
    min.support = 5
    if (length(cc))
      cc = cc %Q% (support>=min.support)
      
  
    uchid = unique((cc %Q% (support >= min.support))$chid)


    mc.cores=10
    
#####INCLUDE SOME SAVING FUNCTION THING HERE
    #for(chid.iter in uchid) {
    #    G_sub = subgraph(G, V(G)[V(G)$chid==chid.iter])
    #    eigen = eigen_centrality(G_sub)
    #    eigen.dt = data.table(cidis=V(G_sub)$cid, centrality=eigen$vector)
    #    cid.central = eigen.dt[centrality>0]$cid
    #    this.cc = cc %Q% (cid %in% cid.central)
    #    peaks = gr.sum(this.cc + 1e3) %>% gr.peaks('score')
    #    binset = bins[, c()] %&% (peaks[peaks$score > quantile(peaks$score, 0.85)])    
    #    binset = gr.stripstrand(binset)
    #    ppdf(plot(c(gTrack(split(this.cc, this.cc$cid) %>% unname), gTrack(binset), gTrack(myc_range)), 'chr8:127035434-128342951'))
    #}
    binsets = pbmclapply(uchid, mc.cores = mc.cores, function(this.chid)
    {
        suppressWarnings({
            this.cc = cc %Q% (chid == this.chid)
            peaks = gr.sum(this.cc + 1e3) %>% gr.peaks('score')
            binset = bins[, c()] %&% (peaks[peaks$score > quantile(peaks$score, 0.85)])
            if (length(binset))
            {
                binset$chid = this.chid
                binset$bid = this.chid
                binset$winid = this.cc$winid[1]
                binset = gr.reduce(binset)
            }
        })
        binset
    })  %>% do.call(grbind, .)
    chrom.object = Chromunity(concatemers = cc[cc$chid %in% binsets$chid], binsets = binsets, meta = data.table())
    return(chrom.object)
}

make_chromunities_from_knn_graph = function(edgelist, reads2, k.knn=10, k.min=1, jaccard.cutoff=NULL, bins) {
    #reads2$binid = factor(reads2$binid)
    setkey(reads2, binid)
    ubid = reads2[, .N, by = binid][N>1 , binid]
 
    
    reads2 = reads2[.(ubid), ]
    ## added for subsamp
    reads2$binid = factor(reads2$binid) 

    reads2$cid = factor(reads2$cid)
    ucid = levels(reads2$cid)

    edgelist_2 = copy(edgelist)
    edgelist_2$bx2 = edgelist$bx1
    edgelist_2$bx1 = edgelist$bx2
    edgelist_3 = rbind(edgelist, edgelist_2)
    edgelist_3$nmat = edgelist_3$mat
    edgelist_3$nfrac = edgelist_3$frac
    setkeyv(edgelist_3, c("nfrac", "nmat"))
    edgelist_3 = unique(edgelist_3)
    edgelist_3 = edgelist_3[order(nfrac, nmat, -bx2, decreasing = T)]

    #fuck knn all my homies hate knn
    k=k.knn
    #browser()
    if(!is.null(jaccard.cutoff)){
        knn.dt = edgelist_3[frac > jaccard.cutoff, .(knn = bx2, weight=(frac * (1/jaccard.cutoff))), by = bx1][!is.na(knn), ]    
    }
    else {
        knn.dt = edgelist_3[mat > 2 & tot > 2, .(knn = bx2[1:k]), by = bx1][!is.na(knn), ]
    }
    knn = sparseMatrix(knn.dt$bx1, knn.dt$knn, x = 1, dims = c(max(knn.dt), max(knn.dt)))

    knn.shared = knn %*% knn
    
    KMIN = k.min
    A = knn.shared * sign(knn.shared > KMIN)
    A[cbind(1:nrow(A), 1:nrow(A))] = 0
                                        #  A <- as(A, "sparseMatrix")
    A = A + t(A)
    A[1,1] = 0
    rownames(A) = 1:dim(A)[1]
    G = graph.adjacency(A, weighted = TRUE, mode = "undirected", add.rownames=NULL)

    print(sum(E(G)$weight == 0))
    #cidis = as.integer(colnames(A))
    cidis = as.integer(igraph::V(G))
    #seqcount = reads2[cidi %in% as.integer(igraph::V(G)),]
    #seqcount = seqcount[, .N, by=c('cidi', 'seqnames')]
    #seqcount[, is_multichromosomal := length(seqnames) > 1, by = cidi]
    #seqcount = unique(seqcount[order(N, decreasing=T)], by='cidi')
    
    #setkey(seqcount, 'cidi')
    
    #G = igraph::set_vertex_attr(G, 'seqnames', value=as.character(seqcount[cidis]$seqnames))
    #G = igraph::set_vertex_attr(G, 'is_multichromosomal', value=as.integer(seqcount[cidis]$is_multichromosomal))

                                        #cl.l = cluster_fast_greedy(G)
    isolated = which(degree(G)==0)
    G = delete.vertices(G, isolated)
    #print(length(unique(reads2_god[(reads2_god$cidi %in% V(G)$name)]$cidi)))
    cl.l = cluster_leiden(G, objective_function='modularity')


    
    cl = cl.l$membership
    ## rename so chid has most support
    cls = 1:max(cl)
    names(cls) = cl %>% table %>% sort %>% rev %>% names
    cl = cls[as.character(cl)]

    cidi.cid = setkey(unique(reads2[, c('cid','cidi')], by=c('cid','cidi')), 'cidi')
    dt = data.table(nodes=V(G)$name) 
    cid.dt = cidi.cid[as.integer(dt$nodes)]
    memb.dt = data.table(cid=cid.dt$cid, chid = cl)

    G = igraph::set_vertex_attr(G, 'chid', value=memb.dt$chid)
    G = igraph::set_vertex_attr(G, 'cid', value=memb.dt$cid)
    browser()
    #if (!is.null(filename))
    #{
    #    cmessage('writing network file')
    #    write_graph(G, filename, format='gml')
    #}

    reads2[, cid := as.character(cid)]
    reads2 = merge(reads2, memb.dt, by = "cid")
    reads2[, `:=`(support, length(unique(cid))), by = chid]
    reads2 = dt2gr(reads2)

    cc = reads2
    min.support = 5
    if (length(cc))
      cc = cc %Q% (support>=min.support)
      
    #browser()
    uchid = unique((cc %Q% (support >= min.support))$chid)
    #print(length(unique(reads2_god[(reads2_god$cidi %in% cc[cc$chid %in% uchid]$cidi)]$cidi)))



#####INCLUDE SOME SAVING FUNCTION THING HERE
    mc.cores=10
    binsets = pbmclapply(uchid, mc.cores = mc.cores, function(this.chid)
    {
        suppressWarnings({
            this.cc = cc %Q% (chid == this.chid)
            peaks = gr.sum(this.cc + 1e3) %>% gr.peaks('score')
            binset = bins[, c()] %&% (peaks[peaks$score > quantile(peaks$score, 0.85)])
            if (length(binset))
            {
                binset$chid = this.chid
                binset$bid = this.chid
                binset$winid = this.cc$winid[1]
                binset = gr.reduce(binset)
            }
        })
        binset
    })  %>% do.call(grbind, .)
    chrom.object = Chromunity(concatemers = cc[cc$chid %in% binsets$chid], binsets = binsets, meta = data.table())
    return(list(chrom.object, G))
}

## aggregate_synergy_results = function(folder) {
##     all.chr = c(paste0("chr", c(as.character(1:22), "X")))
##     syn.list = list()
##     for(chr in all.chr) {
##         folder_name = sub('chromosome', chr, folder)
##         print(paste0(folder_name, 'synergy_results.rds'))
##         if(!file.exists(paste0(folder_name, 'synergy_results.rds'))){
##             next
##         }
##         synergy_result = readRDS(paste0(folder_name, 'synergy_results.rds'))
##         synergy_result[, sig := ifelse(fdr < 0.1, "yes", "no")]
##         synergy_result[, log.est := log2(estimate)]
##         synergy_result$chr = chr
##         syn.list = append(syn.list, list(synergy_result))
##     }
##     return(rbindlist(syn.list))
## }

##let's try to run giga chromunity interchromosomally
evaluate_synergy_interchr = function(concatemers, leave_out_concatemers, chid.to.test, chromosome = NULL, filter_binsets = TRUE, folder = NULL, rebin_thresh = 0.85, mc.cores = 20, numchunks = mc.cores*200 + 1) {
    this.chr = chromosome
    bands.td = gTrack::karyogram(file = "/gpfs/commons/groups/imielinski_lab/DB/UCSC/hg38.cytoband.txt")
    bands = bands.td@data
    bands = grl.unlist(do.call(`GRangesList`, bands))
    cen = bands %Q% (stain=="acen")
    this.max = GRanges(paste0(this.chr, ":", hg_seqlengths()[sub("chr*","",this.chr)]-1e6, "-",  hg_seqlengths()[sub("chr*","",this.chr)]))         
    this.min = GRanges(paste0(this.chr, ":", "1-1e6"))                                                                                                                       
    this.cen = (cen %Q% (seqnames %in% this.chr))+1e6
    this.bad = c(this.min, this.cen, this.max) 
    if (!dir.exists(folder)) {
        stop("output folder does not exist")
    }
    if(is.null(chromosome)){
        chromosome = c(paste0("chr", c(as.character(1:22), "X")))
    }
    resolution = 1e4
    tiles = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), resolution)
    gc5b = readRDS("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/gc.38.rds")
    ## create a list of covariates
    cov_list = list(gc5b)
    ## Specify what kind of covariate it is. Score will be aggregated over bins while the number of intervals will be calculated otherwise.
    names(cov_list) <- c("score:gc.cov")
    ## Make the covariate object
    gc_cov = covariate(name = c("gc"), type = c("numeric"), field = c("score"), data = cov_list)
    gc.cov = gc_cov
    this.dat.chrom = data.table()
    this.chrom = gr2dt(concatemers)
    if (filter_binsets) {
        this.chrom = this.chrom[support > summary(unique(this.chrom[, .(support, chid)])$support)[3]]####filtering out all bins below median support
    }
    this.chrom.w = unique(chid.to.test)
    this.chrom.dt = rebin_community(concatemers, this.chrom.w, resolution=resolution)
    this.chrom.dt[, cardinality := .N, by = chid]
    this.chrom.dt = na.omit(this.chrom.dt)
    ##
    this.all.dat = copy(this.chrom.dt)
    this.all.dat = this.all.dat[cardinality > 2]
    this.all.dat[, bid := chid]
    this.all.dat = this.all.dat[seqnames %in% chromosome]
    this.all.dat[, overall.cardinality := cardinality, by = bid]
    this.chrom.card = unique(this.all.dat[, .(overall.cardinality, bid)])
    ######HELLL NAH
    #this.all.dat = this.all.dat[cardinality < 100]
    ####
    this.sub.parq = leave_out_concatemers
    this.sub.parq$cid = this.sub.parq$read_idx
    this.all.dat[, binid := .I]
    ##
    ###filtering out the bad regions
    ###this makes sense why it wouldn't be here for RE chromunity cause you're only looking at annotated regions anyway
    this.all.dat$bid = this.all.dat$chid

    ##only drop the binids in bad regions, keep the binsets as a whole perhaps
    this.bad.binids = unique(as.character((dt2gr(this.all.dat) %&% (this.bad))$binid))                                                 
    this.all.dat = this.all.dat[!binid %in% this.bad.binids]  
    #browser()
    debug(annotate)
    this.chrom.dat = annotate(binsets = dt2gr(this.all.dat[, bid := chid]),
                              k = 3,
                              concatemers = this.sub.parq,
                              covariates = gc.cov, resolution = resolution,
                              mc.cores = mc.cores, numchunks = numchunks)
    this.chrom.dat.2 = this.chrom.dat
    this.chrom.dat = merge(this.chrom.dat, this.chrom.card[, bid := as.factor(bid)], by = "bid")
    this.chrom.dat[, annotation := "chromunity"]
    set.seed(198)
    message("generating background binsets for the model")
##
    #back.dt = sliding_window_background(chromosome = chromosome, binsets = dt2gr(this.all.dat), n = 1000, resolution = resolution)#, mc.cores=mc.cores)
    browser()
    back.dt = re_background(binsets = dt2gr(this.all.dat), resolution = resolution, n=dim(this.all.dat)[1])#, mc.cores=mc.cores)
    back.dt = back.dt[start != 1]
    upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T)
    setkeyv(back.dt, c("seqnames", "start"))
    back.dt[, V1 := NULL]
    back.dt = na.omit(back.dt)
    back.dt = back.dt[!bid %in% back.dt[width < (resolution-1)]$bid]
    back.dt = gr2dt(gr.reduce(dt2gr(back.dt), by = "bid"))
    back.dt$bid <- as.factor(back.dt$bid)
    back.dt = merge(back.dt, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
    back.dt = back.dt[end < V2][start < V2]
    back.dt[, overall.cardinality := .N, by = bid]
    back.dt = back.dt[overall.cardinality > 1]
    ##
    this.card = unique(back.dt[, .(bid, overall.cardinality)])
    back_gr = dt2gr(back.dt)
    message("extracting background binsets distances")
    this.back.train.dat = annotate(binsets = dt2gr(back.dt),
                                   interchromosomal.table = NULL, #all.hr.dt.mean,
                                   gg = NULL,
                                   concatemers = this.sub.parq, k = 3,
                                   covariates = gc.cov, resolution = resolution,
                                   mc.cores = mc.cores, numchunks = numchunks)
    this.back.train.dat = merge(this.back.train.dat, this.card[, bid := as.factor(bid)], by = "bid")
    this.back.train.dat[, annotation := "random"]
    this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][width <= resolution]$bid)]
    this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][min.dist < resolution]$bid)]
    back.dt[, binid := .I]
    this.bad.train = unique(as.character((back_gr %&% (this.bad))$bid))
    this.back.train.dat = this.back.train.dat[!bid %in% this.bad.train]
    this.back.train.dat = this.back.train.dat[!bid %in% this.back.train.dat[, .(sum(count)), by = bid][V1 == 0]$bid]
####    
    back.model = fit(na.omit(this.back.train.dat)[sum.counts > 0][, setdiff(names(this.back.train.dat), c('overall.cardinality', 'chr', 'annotation')), with = F])

    this.chrom.dat = sscore(this.chrom.dat, model = back.model.ep)
    message("generating random binsets for testing")
    n.chrom = length(unique(this.chrom.dat$bid))
    this.all.dat = this.all.dat[seqnames %in% chromosome]
    back.test = re_background(binsets = dt2gr(this.all.dat), resolution = resolution, n=dim(this.all.dat)[1]*3)#, mc.cores=mc.cores)
    #back.test = sliding_window_background(binsets = dt2gr(this.all.dat), chromosome = chromosome,
    #                                      n = n.chrom*3,
    #                                      resolution = resolution)
    back.test[, V1 := NULL]
    back.test = na.omit(back.test)
    setkeyv(back.test, c("seqnames", "start"))
    back.test[, start := ifelse(start < 0, 0, start)]
    back.test = back.test[!bid %in% back.test[width < (resolution-1)]$bid]
    back.test = gr2dt(gr.reduce(dt2gr(back.test), by = "bid"))
    back.test = merge(back.test, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
    back.test = back.test[end < V2][start < V2]
    back.test[, overall.cardinality := .N, by = bid]
    back.test = back.test[overall.cardinality > 1]
    back_test = dt2gr(back.test)
###
    back.test$bid <- as.factor(back.test$bid)
    this.back.card = unique(back.test[, .(bid, overall.cardinality)])
    bid.back = unique(back.test$bid)
    message("extracting random binsets distances") 
    this.back.test.dat = annotate(binsets = dt2gr(back.test), k=sub.binset.order,
                                  concatemers = this.sub.parq,
                                  covariates = gc.cov, resolution = resolution,
                                  mc.cores = mc.cores, numchunks=numchunks)
    this.back.test.dat = merge(this.back.test.dat, this.back.card, by = "bid")
    this.back.test.dat[, annotation := "random"]
    this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][width <= resolution]$bid)]
###
    this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][min.dist < resolution]$bid)]
####
####
    this.bad.dat = unique(as.character((back_test %&% (this.bad))$bid))
    this.back.test.dat = this.back.test.dat[!bid %in% this.bad.dat]
    back_test = gr2dt(back_test)[bid %in% unique(this.back.test.dat$bid)]
    #this.back.test.dat = sscore(this.back.test.dat, model = back.model)
####
    back.test[, binid := .I]
#####
    set.seed(178)
    all.bid = unique(this.all.dat$bid)
###
    #browser()
    message("starting random walks")
    message("extracting shuffled binsets distances")
    this.all.dat.shuff5 = pbmclapply(1:length(all.bid), function(nr){
         this.clust = dt2gr(this.all.dat[bid == all.bid[nr]])
         this.chrs = .chr2str(as.character(unique(seqnames(this.clust))))
         this.clust.wind = gr.reduce(this.clust+2e6)
         upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T)
         this.clust.wind = (gr2dt(this.clust.wind)[, start := ifelse(start < 0, 1, start)])
         this.clust.wind = merge(this.clust.wind, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
         this.clust.wind[, end := ifelse(end > V2, V2, end)]
         this.clust.wind = dt2gr(this.clust.wind)
         this.sub.win = gr2dt(this.sub.parq %&% this.clust.wind)
         this.sub.win[, new.count := .N, by = read_idx]
         ##
         card = unique(this.sub.win[, .(read_idx, new.count)])
         this.steps = sum(card$new.count)
         this.tiles.orig = gr.tile(this.clust.wind, resolution)
         if (nrow(this.sub.win) > 0){
             this.tgm = cocount(dt2gr(this.sub.win), bins = (this.tiles.orig), by = 'read_idx', full = T)
             A = this.tgm$mat %>% as.matrix
             rownames(A) <- NULL
             colnames(A) <- NULL
             A[cbind(1:nrow(A), 1:nrow(A))] = 0
             A = A + t(A)
             An = A 
             An = round(10*An/min(An[An>0]))  ##can you remove this background 1 value? does this change anything. peculiar
             An = round(1+10*An/min(An[An>0]))
             edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
             ##
             G = graph.edgelist(edges[, cbind(row, col)])
             RW = random_walk(G, start = 200, steps = sum(card$new.count)) %>% as.numeric
             ##
             rm(G)
             rm(edges)
             gc()
             out = this.tgm$gr[RW]%>% gr2dt()
             out$read_idx = card[, rep(read_idx, new.count)] 
             out[, bid := all.bid[nr]]
             out[, cid := read_idx]
             
             this.chrom.sh.dat = annotate(binsets = this.clust, verbose = F,
                                          k = 3,
                                          concatemers = dt2gr(out.old),
                                          covariates = gc.cov, resolution = resolution, 
                                          mc.cores = 10, numchunks = numchunks)
         } else {
             this.chrom.sh.dat = data.table(NA)
         }
         return(this.chrom.sh.dat)
    }, mc.cores = 10, mc.preschedule = T)
####
##     bid.chrom = unique(this.all.dat$bid)
##     chrom.this.chr.shuff5 = pbmclapply(1:length(bid.chrom), function(nr){
##         this.clust = dt2gr(this.all.dat[bid == bid.chrom[nr]])
##         #this.clust.wind = .collapse_gr(this.clust) + ((window.size))
##         #this.sub.win = gr2dt(this_gr_testing %&% this.clust.wind)
##         this.sub.win[, new.count := .N, by = read_idx]
##         ##
##         card = unique(this.sub.win[, .(read_idx, new.count)])
##         this.steps = sum(card$new.count)
##         this.tiles.orig = gr.tile(this.clust.wind, resolution)
##         if (nrow(this.sub.win) > 0){
##             this.tgm = cocount(dt2gr(this.sub.win), bins = (this.tiles.orig), by = 'read_idx', full = T)
##             A = this.tgm$mat %>% as.matrix
##             rownames(A) <- NULL
##             colnames(A) <- NULL
##             A[cbind(1:nrow(A), 1:nrow(A))] = 0
##             A = A + t(A)
##             An = A
##             An = round(1+10*An/min(An[An>0]))
##             edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
##             ##
##             G = graph.edgelist(edges[, cbind(row, col)])
##             RW = random_walk(G, start = 1, steps = sum(card$new.count)) %>% as.numeric
##             ##
##             out = this.tgm$gr[RW]%>% gr2dt()
##             out$read_idx = card[, rep(read_idx, new.count)]
##             out[, bid := bid.chrom[nr]]
##         } else {
##             out = data.table(NA)
##         }
##         return(out)
##     }, mc.cores = 20)
##     kill.zombies()
##     message("running synergy")
## ###
##     chrom.this.chr.shuff5.ne = chrom.this.chr.shuff5[sapply(chrom.this.chr.shuff5, function(x) !inherits(x, "try-error"))]
##     chrom.this.chr.shuff.dt = rbindlist(chrom.this.chr.shuff5.ne, fill = T)
##     chrom.this.chr.shuff.dt[, `:=`(V1 = NULL)]
##     chrom.this.chr.shuff.dt = na.omit(chrom.this.chr.shuff.dt)
##     chrom.this.chr.shuff.dt[, unique.chr.comm := paste0(seqnames, "_", bid)]
##     chrom.this.chr.shuff.dt[, cid := read_idx]
##     set.seed(198)
## ############
##     message("extracting shuffled binsets distances")
##     the.shuff.list = pbmclapply(1:length(bid.chrom), function(x){
##         this.clust = dt2gr(this.all.dat[bid == bid.chrom[x]])
##         chromosome = unique(as.character(seqnames(this.clust)))
##         if (length(this.clust) %in% c(3:10)){
##             this.clust.wind = (.collapse_gr(this.clust)+((window.size)))
##             if (length(this.clust.wind %&% this.bad)  == 0){
##                 this.chrom = tryCatch(annotate(concatemers = dt2gr(chrom.this.chr.shuff.dt[bid == bid.chrom[x]]), binsets = this.clust, covariates = gc_cov, verbose = F, resolution = resolution, k = 5, mc.cores = 1), error = function(e) NULL)
##                 this.chrom = sscore(this.chrom, model = back_model)
##                 }
##             }
##         },  mc.cores = 20)
##         ##kill.zombies()
## ###
    #this.shuff.chrom = tryCatch(rbindlist(the.shuff.list, fill = T), error = function(e) NULL)
    #browser()
    this.all.dat.shuff5.ne = this.all.dat.shuff5[sapply(this.all.dat.shuff5, function(x) !inherits(x, "try-error"))]
    this.all.sh = rbindlist(this.all.dat.shuff5.ne, fill =T)
###
    this.all.sh = sscore(this.all.sh, model = back.model)
    sh.chrom = tryCatch(synergy(binsets = dt2gr(this.all.dat), annotated.binsets = this.all.sh, model = back.model), error = function(e) NULL)
    sh.chrom$fdr = signif(p.adjust(sh.chrom$p, "BH"), 2)
    ##
#####
    this.back.test.dat = sscore(this.back.test.dat, model = back.model)
    theta = back.model$model$theta
    s.chrom = synergy(binsets = dt2gr(this.all.dat), #theta = back.model$model$theta,
                      annotated.binsets = this.chrom.dat, model = back.model)
    s.chrom$fdr = signif(p.adjust(s.chrom$p, "BH"), 2)
    b.chrom = synergy(binsets = dt2gr(back.test), #theta = back.model$model$theta,
                      annotated.binsets = na.omit(this.back.test.dat), model = back.model)
    b.chrom$fdr = signif(p.adjust(b.chrom$p, "BH"), 2)
    ##
    synergy.inter.EP = rbind(s.chrom[, annotation := "chromunity"],
                             b.chrom[, annotation := "random"], 
                             sh.chrom[, annotation := "shuffled"], fill = T)
    if (!is.null(folder)) {
        saveRDS(this.chrom.dat, paste0(folder,'chrom_annotate.rds'))
        saveRDS(this.back.test.dat, paste0(folder,'back_annotate.rds'))
        saveRDS(this.all.sh, paste0(folder,'shuffled_annotate.rds'))
        saveRDS(this.all.dat, paste0(folder,'binsets.rds'))
        saveRDS(back.model, paste0(folder,'back_model.rds'))
        saveRDS(synergy.inter.EP, paste0(folder,'synergy_results.rds'))
    }
    return(synergy.inter.EP)
}




#####This function takes the output of the make_line_graph function and 
#####constructs a pairwise contact map out of the virtual concatemers corresponding
#####to edges in the line graph.
hyper_concatemers = function(edgelist, list.agg, bins, set.counts=NULL, binid = NULL) {
    ##edgelist = line.graph.output[[1]]
    ##reads2 = line.graph.output[[2]]
    ##list.agg = line.graph.output[[3]]

    edgelist = edgelist[1:dim(list.agg)[1]]
    edgelist$j = list.agg$j

    ##list.agg[, hashes := hash(V1)]
    list.agg$hashes = unlist(list.agg$hashes)
    list.agg.count = list.agg[, .(count = .N), by='hashes']
    if(!is.null(set.counts)){
        list.agg.count = set.counts[, c('hashes','num.concats')]
        colnames(list.agg.count)[2] = 'count'
        
    }
    
    unique.list.agg = unique(list.agg, by='hashes') ##with the actual bin labels
    unique.list.agg = merge(unique.list.agg, list.agg.count, by='hashes')
    
    unlist.bins = unique.list.agg[, unlist(V1), by='hashes'] ###nominate "binsets"

    binsets = bins[unlist.bins$V1]
    binsets$hashes = unlist.bins$hashes
    binsets$bid = binsets$hashes    

    unique.bins = unique(unlist.bins$V1)
    dt = data.table(bins=unique.bins)
    dt = dt[order(bins)]

    unique.list.agg[, num := 1:dim(unique.list.agg)[1]]

    asdf = unique.list.agg[, .(unlist(V1), count), by='num']
    if(!is.null(binid)){
        asdf = asdf[V1 != binid]
    }
    
    ##concatm = sparseMatrix(asdf$num, asdf$V1, x=1)


    ## numchunks=100
    ## unique.num = unique(asdf$num)
    ## unum = split(unique.num, ceiling(runif(length(unique.num))*numchunks)) ## randomly chop up cids
    ## mc.cores = 10
    ## nums = unum[[1]]
    ## setkey(asdf,'num')
    ## hyper.edgelist = pbmclapply(unum, mc.cores=mc.cores, mc.preschedule=TRUE, function(nums) {
    ##     dt = asdf[.(nums)]
    ##     dt.merged = merge(dt, asdf, by='num', how='left', allow.cartesian = TRUE)
    ##     hyper.edgelist = dt.merged[, .N, by=c('V1.x','V1.y')]
    ##     hyper.edgelist = hyper.edgelist[V1.x != V1.y]
    ##     #colnames(hyper.edgelist)[colnames(hyper.edgelist) == "V1.x"] = "bx1"
    ##     #colnames(hyper.edgelist)[colnames(hyper.edgelist) == "num.y"] = "bx2"
    ##     #colnames(hyper.edgelist)[colnames(hyper.edgelist) == "N"] = "mat"
    ##     #p1 = concatm[edgelist[, bx1],]
    ##     #p2 = concatm[edgelist[, bx2],]
    ##     #total = Matrix::rowSums(p1 | p2)
    ##     #edgelist$tot = total
    ##     #edgelist[, frac := mat/tot]
    ##     #edgelist = edgelist[mat >= hyperedge.thresh & tot >= hyperedge.thresh]
    ##     return(hyper.edgelist)
    ## })

    ##hyper.edgelist = rbindlist(hyper.edgelist)
    #hyper.edgelist
    bin.adj = merge(asdf, asdf[, .(num,V1)], by='num', allow.cartesian=TRUE)
    bin.adj[, num := NULL]
    bin.adj = bin.adj[V1.x != V1.y]
    bin.adj = bin.adj[, .(full.count = sum(count)), by=c('V1.x','V1.y')]
    bin.sparse = sparseMatrix(bin.adj$V1.x, bin.adj$V1.y, x=bin.adj$full.count)
    matrix = gM(bins, dat=bin.sparse)
    return(matrix)
}


contact_matrix = function(line.graph.output, bins) {
    edgelist = line.graph.output[[1]]
    reads2 = line.graph.output[[2]]
    list.agg = line.graph.output[[3]]

    edgelist = edgelist[1:dim(list.agg)[1]]
    edgelist$j = list.agg$j

    list.agg[, hashes := hash(V1)]
    list.agg$hashes = unlist(list.agg$hashes)
    list.agg.count = list.agg[, .N, by='hashes']
    unique.list.agg = unique(list.agg, by='hashes')

    unlist.bins = unique.list.agg[, unlist(V1), by='hashes']
    binsets = bins[unlist.bins$V1]
    binsets$hashes = unlist.bins$hashes
    binsets$bid = binsets$hashes    

    unique.bins = unique(unlist.bins$V1)
    dt = data.table(bins=unique.bins)
    dt = dt[order(bins)]

    unique.list.agg[, num := 1:dim(unique.list.agg)[1]]

    asdf = unique.list.agg[, .(unlist(V1)), by='num']
    ##concatm = sparseMatrix(asdf$num, asdf$V1, x=1)


    ## numchunks=100
    ## unique.num = unique(asdf$num)
    ## unum = split(unique.num, ceiling(runif(length(unique.num))*numchunks)) ## randomly chop up cids
    ## mc.cores = 10
    ## nums = unum[[1]]
    ## setkey(asdf,'num')
    ## hyper.edgelist = pbmclapply(unum, mc.cores=mc.cores, mc.preschedule=TRUE, function(nums) {
    ##     dt = asdf[.(nums)]
    ##     dt.merged = merge(dt, asdf, by='num', how='left', allow.cartesian = TRUE)
    ##     hyper.edgelist = dt.merged[, .N, by=c('V1.x','V1.y')]
    ##     hyper.edgelist = hyper.edgelist[V1.x != V1.y]
    ##     #colnames(hyper.edgelist)[colnames(hyper.edgelist) == "V1.x"] = "bx1"
    ##     #colnames(hyper.edgelist)[colnames(hyper.edgelist) == "num.y"] = "bx2"
    ##     #colnames(hyper.edgelist)[colnames(hyper.edgelist) == "N"] = "mat"
    ##     #p1 = concatm[edgelist[, bx1],]
    ##     #p2 = concatm[edgelist[, bx2],]
    ##     #total = Matrix::rowSums(p1 | p2)
    ##     #edgelist$tot = total
    ##     #edgelist[, frac := mat/tot]
    ##     #edgelist = edgelist[mat >= hyperedge.thresh & tot >= hyperedge.thresh]
    ##     return(hyper.edgelist)
    ## })

    ##hyper.edgelist = rbindlist(hyper.edgelist)
    #hyper.edgelist
    bin.adj = merge(asdf, asdf, by='num', allow.cartesian=TRUE)
    bin.adj[, num := NULL]
    bin.adj = bin.adj[V1.x != V1.y]
    bin.adj = bin.adj[, .(count = .N), by=c('V1.x','V1.y')]
    bin.sparse = sparseMatrix(bin.adj$V1.x, bin.adj$V1.y, x=(bin.adj$count))
    matrix = gM(bins, dat=bin.sparse)
    return(matrix)
}



shuffle_concatemers = function(concatemers, contact_matrix) {
    A = contact_matrix$mat %>% as.matrix
    rownames(A) <- NULL
    colnames(A) <- NULL
    A[cbind(1:nrow(A), 1:nrow(A))] = 0
    A = A + t(A)
    An = A 
    ##browser()
    An = round(10*An/min(An[An>0]))  ##can you remove this background 1 value? does this change anything. peculiar     An = round(1+10*An/min(An[An>0]))
    edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
     ##
    G = graph.edgelist(edges[, cbind(row, col)])
    ##G.2 = graph.adjacency(An)#, directed=FALSE)
    
##    browser()
    concats.binned = bin_concatemers(concatemers, contact_matrix$gr)
    concats.dt = concats.binned
    concats.dt = unique(concats.dt, by=c('binid','cidi')) ###dedupe!!

    concat.counts = concats.dt[, new.count := .N, by='read_idx']


    card = unique(concat.counts[, .(read_idx, new.count)])

    this.steps = sum(card$new.count)

    row.labels = 1:dim(An)[1]
    start.index = row.labels[rowSums(An)>1]
    
    RW = random_walk(G, start = start.index, steps = sum(card$new.count)) %>% as.numeric
##    browser()
    ##
    ##rm(G)
    ##rm(edges)
    ##gc()
    out = contact_matrix$gr[RW] %>% gr2dt()
    out$read_idx = card[, rep(read_idx, new.count)] 


    ## ####handling multi bin hits
    ## ##browser()
    ## bins = contact_matrix$gr
    ## bins$binid = 1:length(bins)
    ## bins.dt = gr2dt(bins)
    ## out = merge(out, bins.dt[, c('start','end','binid')], by=c('start','end'))
    ## rewalk = out[, .N, by=c('read_idx','binid')][N>1]
    ## RW.rewalk = RW
    ## out.old = out

    
    ## while(dim(rewalk)[1] > 1) { ###literally just keep walking if you have a multi bin hit
    ##     print(rewalk)
    ##     out = out[!(read_idx %in% rewalk$read_idx)]
    ##     card.rewalk = merge(rewalk, card, by='read_idx')
    ##     end.index = RW.rewalk[length(RW.rewalk)]
    ##     RW.rewalk = random_walk(G, start = end.index, steps = sum(card.rewalk$new.count)) %>% as.numeric
    ##     ##out.new =
        
    ##     out.new = contact_matrix$gr[RW.rewalk] %>% gr2dt()
    ##     out.new$read_idx = card.rewalk[, rep(read_idx, new.count)] 
    ##     out.new = merge(out.new, bins.dt[, c('start','end','binid')], by=c('start','end'))
    ##     print(out.new)
    ##     rewalk = out.new[, .N, by=c('read_idx','binid')][N>1]
    ##     out.new = out.new[!(read_idx %in% rewalk$read_idx)]
    ##     out = rbind(out, out.new)    
    ## }
    ##out[, bid := all.bid[nr]]
    ##out[, cid := read_idx]
    return(out)
}

###takes 
make_overlap_concatemers = function(edgelist, list.agg, pair.thresh = 0) {    
    edgelist.half = edgelist[1:dim(list.agg)[1]]
    edgelist.half$j = list.agg$j
    edgelist.half$hashes = list.agg$hashes
    edgelist.half[, .N, by='hashes']
    edgelist.half[, length(unique(bx1)), by='hashes']


    edgelist_2 = copy(edgelist.half)
    edgelist_2$bx2 = edgelist.half$bx1
    edgelist_2$bx1 = edgelist.half$bx2
    edgelist_3 = rbind(edgelist.half, edgelist_2)

    concat.count = edgelist_3[, .(num.concats = length(unique(bx1))), by='hashes'] ###you can ALSO GET the number of raw concatemers hitting this from the line graph
    pair.count = list.agg[, .N, by='hashes']
    all.counts = merge(pair.count, concat.count, by='hashes')
    all.counts[, diff := N - num.concats]

    unique.list.agg = unique(list.agg, by='hashes')
    unique.unlist.agg = unique.list.agg[, unlist(V1), by='hashes']
    ##FABULOUS FILTER 
    overlap.concats = merge(all.counts[N >= pair.thresh], unique.unlist.agg, by='hashes')
##VERY IMPORTANT:: KEEP ONLY THE OVERLAP CONCATEMERS WITH A MINIMUM SUPPORT OF 3 CONCATEMER PAIRS
##IF THIS MAKES SENSE

    overlap.concats = overlap.concats[, c('hashes', 'V1')]
    overlap.concats$hashes = factor(overlap.concats$hashes)
    hash.cidi.table = data.table(hashes = overlap.concats$hashes, cidi = as.integer(overlap.concats$hashes))
    overlap.concats$hashes = as.integer(overlap.concats$hashes)
    colnames(overlap.concats) = c('cidi','binid')
    return(list(overlap.concats, hash.cidi.table))
}

subset_edges = function(overlap.concats, overlap.line.graph.output) {
    edgelist = overlap.line.graph.output[[1]]
    edgelist.half = edgelist[1:(dim(edgelist)[1] / 2)]
    
    edgelist.copy = copy(edgelist.half)
    edgelist.copy$bx1 = edgelist.half$bx2
    edgelist.copy$bx2 = edgelist.half$bx1
    edgelist.half = rbind(edgelist.half, edgelist.copy)
    edgelist.half$j = 1:dim(edgelist.half)[1]
    list.overlap.concats = overlap.concats[, .(list(binid)), by='cidi']
    
    overlap.edgelist.merge.1 = merge(edgelist.half, list.overlap.concats, by.x='bx1',by.y='cidi', allow.cartesian=TRUE) ###When the binids equal bx1

    overlap.edgelist.merge.2 = merge(edgelist.half, list.overlap.concats, by.x='bx2',by.y='cidi', allow.cartesian=TRUE) ###When the binids equal bx2
    overlap.edgelist.merge.p0 = merge(overlap.edgelist.merge.1, overlap.edgelist.merge.2[, c('j','V1')], by='j', allow.cartesian=TRUE)
    ##overlap.edgelist.merge.p0.unique = unique(overlap.edgelist.merge.p0, by=c('j','binid.x','binid.y'))
    
    subset.comp = overlap.edgelist.merge.p0[, .(sub = unlist(V1.x)), by='j']
    is.subset.dt = overlap.edgelist.merge.p0[, .(is.subset = all(unlist(V1.x) %in% unlist(V1.y))), by=c('j')]
    unique.subsets = is.subset.dt[is.subset==TRUE,] ##%>% unique(by=c('bx1','bx2'))
    subset.edges = merge(edgelist.half, unique.subsets, by='j')
    return(subset.edges)
}

##basal.subsets
create_missing_subsets = function(subset.edges, unique.combo.registry, line.graph.edgelist, list.agg, superset.order=4, subset.order=3) {
    ###unique.unlist.count$hashes = unique.unlist.count$hashes %>% factor ### WE KNOW ALREADY HOW MANY TIMES EACH SUPERSET HAPPENS
    hash.cidi.table = unique.combo.registry[, c('hashes','cidi')]
    ##hash.cidi.table$hashes = factor(hash.cidi.table$hashes)
    ##hash.cidi.table = hash.cidi.table %>% unique(by=c('hashes','cidi'))
    ##cidi.count = merge(unique.combo.registry, hash.cidi.table, by='hashes')
    cidi.count = unique.combo.registry
    unique.cidi.count = unique(cidi.count, by='hashes')
    print("made unique cidi count")
    ##browser()
    trimmed.subsets = subset.edges[, c('bx1','bx2')]
    colnames(trimmed.subsets) = c('subset','superset')
    frequency = unique.cidi.count[, c('cidi','frequency', 'N')]
    super.frequency = frequency[N == superset.order]
    sub.frequency = frequency[N == subset.order]
    
    trimmed.subsets.count = merge(trimmed.subsets, super.frequency, by.x='superset', by.y='cidi', all.y=TRUE)
    colnames(trimmed.subsets.count) = c('superset','subset','superset.frequency','superset.N')
    trimmed.subsets.count = merge(trimmed.subsets.count, sub.frequency, by.x='subset', by.y='cidi', all.x=TRUE)
    
    trimmed.subsets.count[, aggregated.frequency := frequency + superset.frequency]
    
    ##browser()
    trimmed.subsets.count[, num.expected.subsets := choose(superset.order,subset.order)] ###ONLY BETWEEN LAYERS
    trimmed.subsets.count.N3 = trimmed.subsets.count[N==subset.order & superset.N==superset.order]
    trimmed.subsets.count.N3[, num.subsets := .N, by='superset']
    trimmed.subsets.count.N3[, remaining.subsets := num.expected.subsets - num.subsets]
    ##browser()

###brother you need to do a BIG COMBN
    
    
###add condition of trimmed subsets equaling the newly added ones

    ##2147199
    ##browser()
    if(all(trimmed.subsets.count$N %>% is.na)) {
         agg.combos = list()
         big.combn = trimmed.subsets.count$superset %>% as.data.table
         colnames(big.combn) = 'superset'
    } else {
        unique.unlist.count = unique.combo.registry[, unlist(V1), by='hashes']
        unique.unlist.count.cidi = merge(hash.cidi.table, unique.unlist.count[,c('hashes','V1')], by='hashes')
        subset.trimmed.N3 = merge(trimmed.subsets.count.N3[remaining.subsets>0], unique.unlist.count.cidi, by.y='cidi', by.x='subset')
        bincounts = subset.trimmed.N3[, .(bincount = .N), by=c('V1', 'superset')] ##yes yes YES
###For a superset, count the number of times each bin appears amongst subsets 
        unique.supersets = unique(bincounts, by='superset')
        unique.supersets[, bincount := NULL]
        colnames(bincounts)[1] = 'super.bins'
        super.bins = merge(unique.supersets[, c('superset')], unique.unlist.count.cidi, by.x='superset', by.y='cidi')
        colnames(super.bins)[3] = 'super.bins'
        super.comp = merge(super.bins, bincounts, by=c('super.bins','superset'), all.x=TRUE, allow.cartesian=TRUE)
        super.comp[is.na(bincount), bincount:=0]
####YOU CAN CONSTRUCT THE MISSING TRIPLE FROM THIS DATA HERE

        super.comp[, max.bincount := max(bincount), by='superset']
        combo.bins = super.comp[bincount == max.bincount]
        combo.bins[, num.left := .N - 1, by='superset'] ##this is incorrect
###combos = combo.bins[, t(combn(super.bins, num.left[1])) %>% as.data.table, by='superset']
        combos = combo.bins[, .(combn(super.bins, num.left[1], simplify=FALSE)), by='superset'] ###YESYESYESYESYES
        non.combo.bins = super.comp[bincount < max.bincount]
        non.combos = non.combo.bins[, .(list(super.bins)), by='superset']
        agg.combos = merge(combos, non.combos, by='superset')
        print('combos made')
        agg.combos$combs = do.call(Map, c(f = c, agg.combos[, c('V1.x','V1.y')]))
        ##browser()
                                        #this segfaults??? lmao
        ##unlist.agg = agg.combos[, .(unlist(combs)), by='superset']
        ##setorder(unlist.agg, V1)

        ##unlist.agg[, .(list(V1)), by='superset']
            ##browser()
        mask = trimmed.subsets.count$superset %in% agg.combos$superset
    
        big.combn = trimmed.subsets.count[!mask]$superset %>% unique %>% as.data.table

        colnames(big.combn) = 'superset'
        big.combn.merge = merge(big.combn, trimmed.subsets.count.N3, by='superset', all.x=TRUE)
        big.combn.merge[is.na(remaining.subsets), remaining.subsets := 1] ##ensuring we don't do a big combn on sets which have ALL subsets filled already
        big.combn = big.combn.merge[remaining.subsets != 0]$superset %>% as.data.table
        colnames(big.combn) = 'superset'

    }
    print('Doing a BIG COMBN')
    ##browser()

    
    big.combn.bins = merge(big.combn, unique.combo.registry, by.x='superset', by.y='cidi')
    big.combn.bins = big.combn.bins[, unlist(V1), by='superset']
    big.combn.bins = unique(big.combn.bins, by=c('superset','V1')) ##brother you must dedupe
    big.combn = big.combn.bins[, .(combn(V1, subset.order, simplify=FALSE)), by='superset']
    
    ##fuckr()
    if(all(trimmed.subsets.count$N %>% is.na)) {
        colnames(big.combn)[2] = 'combs.sort'
        str.identifiers = big.combn$combs.sort %>% paste0
        str.identifiers.factor = factor(str.identifiers)
        big.combn$hashes = str.identifiers.factor
        big.combn$N = subset.order
        big.combn = unique(big.combn, by='hashes')
        agg.combos = big.combn
    } else {
        
        agg.combos$combs = unname(agg.combos$combs)
        agg.combos$combs.sort = map(agg.combos$combs, sort)
        agg.combos = agg.combos[, c('superset','combs.sort')]
        agg.combos$combs.sort = map(agg.combos$combs.sort, unname)
        ##browser()


        str.identifiers = agg.combos$combs.sort %>% paste0
        str.identifiers.factor = factor(str.identifiers)
        agg.combos$hashes = str.identifiers.factor

        all.combos = agg.combos %>% copy

###browser()
        
        mask = agg.combos$hashes %in% unique.combo.registry$hashes
        agg.combos = agg.combos[!mask]
        ##browser()
        ##agg.combos[, hashes := hash(combs.sort)]

        ##agg.combos$hashes = unlist(agg.combos$hashes)
        ##agg.combos$hashes = agg.combos$hashes %>% factor

        ##xbrowser()

        agg.combos$N = subset.order
        ##agg.combos$j = 1:dim(agg.combos)[1]
        ##agg.combos[, frequency := .N, by='hashes']
        ##agg.combos = unique(agg.combos, by='hashes')
        ##agg.combos.copy = agg.combos %>% copy
        ##agg.combos = agg.combos.copy %>% copy
        
        ##agg.combos[, superset := NULL]
        colnames(agg.combos) = c('superset', 'combs.sort', 'hashes', 'N')

        colnames(big.combn)[2] = 'combs.sort'
        str.identifiers = big.combn$combs.sort %>% paste0
        str.identifiers.factor = factor(str.identifiers)
        big.combn$hashes = str.identifiers.factor
        big.combn$N = subset.order
        
        agg.combos = rbind(agg.combos, big.combn)
    }    

    agg.combos$j = 1:dim(agg.combos)[1]
    agg.combos[, frequency := .N, by='hashes']
    
    
    ##browser()
    edgelist.half = line.graph.edgelist[1:(dim(line.graph.edgelist)/2)]
    edgelist.half$hashes = list.agg$hashes
    edgelist.half = merge(edgelist.half, hash.cidi.table, by='hashes')
    
    edgelist_2 = copy(edgelist.half)
    edgelist_2$bx2 = edgelist.half$bx1
    edgelist_2$bx1 = edgelist.half$bx2
    edgelist_3 = rbind(edgelist.half, edgelist_2)


    ###MERGE UNTIL YOU REACH THE BASAL LAYER

    ###OH YES I REMEMBER

    basal.number = edgelist_3$cidi %>% max
    ##debug(merge_to_basal)

    ##browser()
    basal.subsets = merge_to_basal(basal.number, trimmed.subsets)
    
    ##basal.number = set.counts$cidi %>% max
    ##subsets.basal = trimmed.subsets[superset <= basal.number]
    ##extraneous.subsets = trimmed.subsets[superset > basal.number]
    ##subsets.again = merge(extraneous.subsets, subsets.basal, by.x='superset', by.y='subset')
    agg.combos.backup = agg.combos %>% copy
    agg.combos = agg.combos.backup %>% copy
    
    agg.combos = merge(basal.subsets, agg.combos.backup, by.x='subset', by.y='superset', all.y=TRUE)
    edge.case = agg.combos[subset<basal.number][subset != superset]
    edge.case.copy = edge.case %>% copy
    edge.case$subset = edge.case.copy$superset
    edge.case$superset = edge.case.copy$subset
    edge.case = rbind(edge.case, edge.case.copy)
    agg.combos = rbind(agg.combos, edge.case)

    agg.combos[is.na(superset), superset := subset]
    agg.combos = unique(agg.combos, by=c('superset','hashes'))
    
    
    
    supersets.count = merge(agg.combos[, c('superset','hashes', 'combs.sort')], edgelist_3, by.x='superset', by.y='cidi')
    supersets.count[, num.concats := length(unique(bx1)), by='hashes.x']
    new.combos.count = supersets.count[, c('hashes.x','num.concats', 'combs.sort', 'superset')]
    colnames(new.combos.count)[1] = 'hashes'
    new.combos.count = unique(new.combos.count, by=c('hashes','superset'))
    new.combos.unlist = new.combos.count[, unlist(combs.sort), by='hashes']
    new.combos.unlist = new.combos.unlist[, .N, by='hashes']
    new.combos.count = merge(new.combos.count, new.combos.unlist, by='hashes')
    new.combos.count$j = 1:dim(new.combos.count)[1]
    colnames(new.combos.count)[3] = 'V1'
    new.combos.count[, frequency := choose(num.concats, 2)]
    new.combos.count$cidi = 0
    
    ## overlap.overlap = overlap.line.graph.output[[3]]
    ## unique.overlap.overlap = unique(overlap.overlap, by='hashes')
    ## overlap.overlap.unlist = unique.overlap.overlap[, unlist(V1), by='hashes']
    ## asdf = overlap.overlap.unlist[, .N, by='hashes']
   
    ## overlap.overlap.count = overlap.overlap[, .(overlap.overlap.count = .N), by='hashes']
    ## asdf = merge(overlap.overlap.count, merged, by='hashes', all.y=TRUE)
    ## ##THIS segfaults
    ## ##agg.combos$hashes = unlist(agg.combos$hashes)
    ## agg.combos = unique(agg.combos, by='hashes')

    ## unique.overlap.overlap = overlap.overlap.count
    ## masker = unique.overlap.overlap$hashes %in% list.agg.chr8$hashes
    ## overlap.overlap.trimmed = unique.overlap.overlap[!masker]
    ## mask = overlap.overlap.trimmed$hashes %in% agg.combos$hashes
    
    ## overlap.overlap.trimmed[mask][hashes=='c(1076, 1622, 1623)']
    ## merged[hashes=='c(52, 70, 78)']


    ## overlap.edgelist = overlap.line.graph.output[[1]]
    ## overlap.edgelist.half = overlap.edgelist[1:(dim(overlap.edgelist)[1]/2)]
    ## overlap.edgelist.half$hashes = overlap.overlap$hashes

    ##browser()
    
    

    new.combos.count$cidi = (max(unique.combo.registry$cidi) + 1):(max(unique.combo.registry$cidi) + dim(new.combos.count)[1])
    new.combos.count[, cidi := min(cidi), by='hashes']
    super.combos = new.combos.count %>% copy
    new.combos.count[, superset := NULL]
    new.combos.count$N = subset.order
    new.combos.count = new.combos.count %>% unique(by='hashes') ##cases where you have more than one superset
    unique.combo.registry.agg = rbind(unique.combo.registry, new.combos.count)
    ##unique.combo.registry.agg$cidi = 1:(dim(unique.combo.registry.agg)[1])

    ##super.combos$cidi = max(unique.combo.registry$cidi):(max(unique.combo.registry$cidi) + dim(super.combos)[1] - 1)
    
    ###append to subsets here
    ##new.combos.agg = unique.combo.registry.agg[hashes %in% new.combos.count$hashes]
    ##new.combos.agg = merge(new.combos.agg, super.combos[, c('hashes','superset')], by='hashes')

    all.supers.count = merge(super.combos,agg.combos.backup[, c('hashes','superset')], by='hashes')
    all.supers.count.2 = all.supers.count %>% copy
    all.supers.count.2$superset.x = all.supers.count$superset.y
    all.supers.count.2$superset.y = all.supers.count$superset.x
    all.supers = rbind(all.supers.count, all.supers.count.2)
    all.supers = unique(all.supers, by=c('hashes','superset.x'))
    
    new.subsets = all.supers[, c('j', 'cidi', 'superset.x', 'N')]
    colnames(new.subsets)[2:4] = c('bx1','bx2','mat')
    new.subsets$mat = subset.order
    new.subsets$tot = superset.order
    new.subsets[, frac := mat / tot]
    new.subsets$is.subset = TRUE
    subset.append = rbind(subset.edges, new.subsets)
    
    return(list(unique.combo.registry.agg, subset.append))
}


merge_to_basal = function(basal.number, trimmed.subsets) {
    subsets.basal = trimmed.subsets[superset <= basal.number]
    extraneous.subsets = trimmed.subsets[superset > basal.number]
    subsets.again = merge(extraneous.subsets, subsets.basal, by.x='superset', by.y='subset')
    subsets.again[, superset := NULL]
    colnames(subsets.again)[2] = 'superset'
    if(max(subsets.again$superset) > basal.number){
        extraneous.subsets = subsets.again[superset > basal.number]
        subsets.chunk = merge_to_basal(basal.combos, extraneous.subsets) ##RECURSION!!!
        subsets.again = rbind(subsets.again, subsets.chunk)
    }
    return(rbind(subsets.basal, subsets.again))
}


identify_existing_subsets = function(iterative.registry, old.iterative.registry, subsets, overlap.concats, subset.order=3, superset.order=4) {
    
    new.cidis = iterative.registry$cidi[!(iterative.registry$cidi %in% old.iterative.registry$cidi)]
    new.concats = iterative.registry[cidi %in% new.cidis, c('V1','cidi')]
    new.concats = new.concats[, .(binid = unlist(V1)), by='cidi']
    ##overlap.concats[, yeet := .N, by='cidi']
    sub.concats = overlap.concats[numbins==subset.order]
    sub.concats[, numbins := NULL]
    concats = rbind(new.concats, sub.concats)

    lineg = make_line_graph(concats, bins.chr8, prebinned=TRUE, order=superset.order)
    
    subsets.dunka = subset_edges(concats, lineg)
    additional.subsets = rbind(subsets, subsets.dunka)

    return(additional.subsets)
}


                                    
make_weighted_simplicial_complex = function(concatemers, bins, lineg.output=NULL, overlap.line.graph=NULL) {

    ##start with the line graph
    if(is.null(lineg.output)){
        lineg.output = make_line_graph(concatemers, bins)
##        saveRDS(lineg.output, 'line_graphs/wholegenome_line_graph_50k.rds')
        print("Concatemer line graph made")
    }
    edgelist = lineg.output[[1]]
    reads2 = lineg.output[[2]]
    list.agg = lineg.output[[3]]

    
    overlap.concats.output = make_overlap_concatemers(edgelist, list.agg, pair.thresh=0)
    overlap.concats = overlap.concats.output[[1]]

    if(is.null(overlap.line.graph)){
        overlap.line.graph = make_line_graph(overlap.concats, bins, prebinned = TRUE)
  ##      saveRDS(overlap.line.graph, 'line_graphs/wholegenome_overlap_line_graph_50k.rds')
        print("Overlap concatemer line graph made")
    }
    subset.edges = subset_edges(overlap.concats, overlap.line.graph)

    list.agg.count = list.agg[, .(frequency = .N), by='hashes']
    beginning.combo.registry = unique(list.agg, by='hashes')
    beginning.combo.registry.unlist = beginning.combo.registry[, unlist(V1), by='hashes']
    beginning.combo.registry.unlist[, N := .N, by='hashes']
    beginning.combo.registry = merge(beginning.combo.registry, beginning.combo.registry.unlist[, c('hashes','N')], by='hashes') %>% unique(by='hashes')
    beginning.combo.registry$hashes = beginning.combo.registry$hashes %>% factor
    beginning.combo.registry = merge(beginning.combo.registry, list.agg.count, by='hashes')

    hash.cidi.table = overlap.concats.output[[2]]
    unique.hash.cidi.table = unique(hash.cidi.table, by='hashes')

    beginning.combo.registry.cidi = merge(beginning.combo.registry, unique.hash.cidi.table, by='hashes')
    dt = merge(subset.edges, beginning.combo.registry.cidi[, c('frequency','cidi')], by.x='bx1',by.y='cidi')
    colnames(dt)[8] = 'sub.frequency'
    dt = merge(dt, beginning.combo.registry.cidi[, c('frequency','cidi')], by.x='bx2',by.y='cidi')
    colnames(dt)[9] = 'super.frequency'

    set.counts = beginning.combo.registry.cidi %>% copy

    dt = merge(subset.edges, set.counts[, c('frequency','cidi')], by.x='bx1',by.y='cidi')
    colnames(dt)[8] = 'sub.frequency'
    dt = merge(dt, set.counts[, c('frequency','cidi')], by.x='bx2',by.y='cidi')
    colnames(dt)[9] = 'super.frequency'

    sums = dt[, .(sub.frequency, sum(super.frequency)), by='bx1']
    sums = unique(sums, by='bx1')
    sums[, total.pairs := sub.frequency + V2, by='bx1']
    sums[, num.concats := inverse.choose(total.pairs,2), by='bx1']

    set.counts[sums$bx1, frequency := sums$total.pairs]
    set.counts[, num.concats := inverse.choose(frequency,2), by='cidi'] #very cute
    print("Overlap concatemer Omega values calculated")
    
    beginning.combo.registry.cidi = set.counts

    overlap.concats[, numbins := .N, by='cidi']

    ##we are SO BACK
    iterative.registry = set.counts %>% copy
    iterative.subsets = subset.edges %>% copy
    ##table = iterative.registry[, c('hashes','cidi')]

    beginning.order = iterative.subsets$mat %>% max
    
    print("Beginning iterative expansion of missing subsets")
    for(order in beginning.order:4){
        print(order)
        iteration = create_missing_subsets(iterative.subsets, iterative.registry, edgelist, list.agg, superset.order=order, subset.order=order-1)
        old.iterative.registry = iterative.registry
        old.iterative.subsets = iterative.subsets
        iterative.registry = iteration[[1]]
        iterative.subsets = iteration[[2]]
        print(iterative.registry)
        if((order-2)==2){
            break
        }
        iterative.subsets = identify_existing_subsets(iterative.registry, old.iterative.registry, iterative.subsets, overlap.concats, subset.order=order-2, superset.order=order-1)
    }


    ###time to handle the big concats
    
    big.concats = overlap.concats[numbins >= beginning.order]

    ##real coding

    print("Expanding particularly large overlap concatemers")
    big.combns = pbmclapply((beginning.order-1):3, function(i) {
        combn = big.concats[, .(combn(binid, i, simplify=FALSE)), by='cidi']
        combn$N = i
        return(combn)
    }, mc.cores = 3)

    ##browser()
    big.combns = rbindlist(big.combns)

    ##big.combns$combs = unname(big.combns$V1)
    ##big.combns$combs.sort = map(big.combns$combs, sort)
    
    str.identifiers = big.combns$V1 %>% paste0
    str.identifiers.factor = factor(str.identifiers)
        
    big.combns$hashes = str.identifiers.factor
    big.combns[, count := .N, by='hashes']
    big.combns = unique(big.combns, by='hashes')


    ##big.combns.backup = big.combns %>% copy
    ###OOOOOOOOOOHHHHHHHOHOHOHOHOHO
###HERE IS MY LOGIC MISTAKE
    ##NON-TRIVIAL

    
    mask = big.combns$hashes %in% iterative.registry$hashes
    big.combns = big.combns[!mask]
    
    merged = merge(big.combns[count==1], iterative.registry[, c('cidi', 'num.concats')], by='cidi')
    merged[, cidi := NULL]
    merged[, count := NULL]
    merged[, frequency := choose(num.concats, 2)]
    merged$cidi = 1:dim(merged)[[1]] + (max(iterative.registry$cidi))
    merged$j = 1:dim(merged)[[1]]

    ##iterative.registry.backup = iterative.registry %>% copy
    iterative.registry = rbind(iterative.registry, merged)

    return(list(iterative.registry, iterative.subsets, bins))
}

simplex_degree_filtration = function(iterative.registry, cardinality.input) {
    max.omega = iterative.registry[N==cardinality.input]$num.concats %>% max
    filtration = pbmclapply(2:max.omega, function(x){
        iterative.registry.unlist.filtered = iterative.registry[num.concats >= x, .(unlist(V1)), by='hashes']
        top.degree.dt.filt = iterative.registry.unlist.filtered[, cardinality := .N, by='hashes']
        simplex.3.filt = top.degree.dt.filt[cardinality == cardinality.input, .(degree = .N), by='V1']
        simplex.3.filt$value = x
        return(simplex.3.filt)
    }, mc.cores = 5)

    filtration.dt = rbindlist(filtration)
    filtration.dt$number.bin = filtration.dt$V1
    filtration.dt$V1 = factor(filtration.dt$V1)
    filtration.dt$log.degree = log(filtration.dt$degree)
    filtration.dt$log.omega = log(filtration.dt$value)
    return(filtration.dt)
}

create_degree_track = function(iterative.registry, simplex.cardinality, omega.thresh, bins) {
    iterative.registry.unlist.filtered = iterative.registry[num.concats >= omega.thresh & N == simplex.cardinality, .(unlist(V1)), by='hashes']
    top.degree.dt.filt = iterative.registry.unlist.filtered[, cardinality := .N, by='hashes']
    simplex.filt = top.degree.dt.filt[cardinality == simplex.cardinality, .(degree = .N), by='V1']
    simplex.filt.gr = bins[simplex.filt$V1]
    simplex.filt.gr$degree = simplex.filt$degree
    return(simplex.filt.gr)
}



annotate_simplicies = function(iterative.registry, bin.width, bins, simplex.order = 3) {
    
    iterative.sub = iterative.registry[N==simplex.order]
    iterative.registry.unlist = iterative.sub[, .(unlist(V1), num.concats, N), by='hashes']
    iterative.registry.unlist[, span := max(V1) - min(V1), by='hashes']

##iterative.registry.unlist[V1==2555 & span > 10 & num.concats>3]
##print(setorder(iterative.registry.unlist[V1==2555 & span > 10 & num.concats>3], -num.concats))

##    mins = iterative.registry.unlist[1:10][, .(diff(combn(V1,2,simplify=FALSE))), by='hashes']

    mins = iterative.registry.unlist[, .(combn(V1,2,simplify=FALSE)), by='hashes']

    diffs = iterative.registry.unlist[, .(diffs = diff(V1)), by='hashes']
    min.dist = diffs[, .(min.dist = min(diffs)), by='hashes']
    max.dist = iterative.registry.unlist[, c('hashes', 'span')]
    max.dist = unique(max.dist, by='hashes')

    iterative.sub = merge(iterative.sub, max.dist, by='hashes')
    iterative.sub = merge(iterative.sub, min.dist, by='hashes')


    iterative.sub.unlist = iterative.sub[, .(unlist(V1), num.concats, N, span, min.dist), by='hashes']
    ##iterative.sub.unlist[V1==2555]

    ##mycskis = iterative.registry.unlist[span>10 & min.dist>5 & V1==2555]

#' jorvis Tuesday, Dec 05, 2023 09:25:40 AM
#' aaa

##    setorder(mycskis, -num.concats)

    ##Let's annotate some things

    ##bin.width = 50000

    iterative.sub[, min.dist.genomic := (min.dist-1)*bin.width] 
    iterative.sub[, max.dist := (span-1)*bin.width]    


    ##create this format of table from iterative registry

    iterative.sub[, min.dist := min.dist.genomic]
    
    colnames(iterative.sub)[colnames(iterative.registry)=='N'] = 'cardinality'
    iterative.sub[, width := bin.width * cardinality]

    iterative.sub.unlist = iterative.sub[, .(unlist(V1)), by='hashes']
    ##iterative.registry.unlist
    unlisted = bins[iterative.sub.unlist$V1]
    unlisted$hashes = iterative.sub.unlist$hashes




    ##gc.chr8 = gc5b %Q% (seqnames=='chr8')
    gc.bins = bins %$% gc5b
    gc.bins.dt = gr2dt(gc.bins)
    

    ##iterative.registry.unlist = iterative.registry[, .(unlist(V1)), by='hashes']
    unlisted = gc.bins[iterative.sub.unlist$V1]
    unlisted$hashes = iterative.sub.unlist$hashes
    unlisted.dt = gr2dt(unlisted)
    gc.dt = unlisted.dt[, .(value = mean(score)), by='hashes']
    iterative.sub[, cardinality := NULL]
    iterative.sub[, cardinality := simplex.order]
    iterative.sub = merge(iterative.sub, gc.dt, by='hashes')

    iterative.sub[, min.dist.genomic := NULL]
    return(iterative.sub)
}


make_weighted_simplicial_complex_simpler = function(dt.concats, numchunks=50, cardinality=3, cores=5) {

    
    unique_cidi = dt.concats$cidi %>% unique
    ucidl = split(unique_cidi, ceiling(runif(length(unique_cidi))*numchunks))


    dt.concats.sort = dt.concats[order(binid, cidi)]
    dt.concats.sort[, count := .N, by='cidi']
    dt.concats.sort = dt.concats.sort[count >= cardinality]

   ## cidis = ucidl[[1]]

    combinatorics = pbmclapply(ucidl, mc.cores = cores, function(cidis) {
        combn.chunk = dt.concats.sort[cidi %in% cidis, .(combn(binid,cardinality,simplify=FALSE)), by='cidi']
        combn.chunk$hashes = combn.chunk$V1 %>% paste0
        ##combn.chunk$comb.fact = factor(asdf$comb.str)
        combn.chunk.count = combn.chunk[, .(count = .N, V1), by='hashes']


        ##combn.chunk.count = combn.chunk %>% copy
        ##combn.chunk.count[, count := .N, by='hashes']
        combns = unique(combn.chunk.count, by='hashes')
        ##combns = combn.chunk.count %>% copy
        return(combns)
    })

    dt = rbindlist(combinatorics)
    simplicial = dt[, .(num.concats = sum(count)), by='hashes']

    hashlist = dt[, c('hashes','V1')] %>% unique(by='hashes')
    simplicial.registry = merge(simplicial, hashlist, by='hashes')
    simplicial.registry$N = cardinality
    simplicial.registry$sidi = 1:dim(simplicial.registry)[[1]]
    return(list(simplicial.registry, dt))
}


annotate_distance_decay = function(concats.gr, bins.gr, min.value=50, window.size=50, min.simplex = 1, simplicial.complex = NULL, pairwise=NULL, iterative.registry.unlist.filterable=NULL, return.training=FALSE, numchunks=500) {

    simplex.order = 3
    if (is.null(simplicial.complex) & is.null(pairwise)){
        binned.concats = bin_concatemers(concats.gr, bins = bins.gr, max.slice = 1e6, mc.cores=5, verbose=TRUE, hyperedge.thresh = 3)
        
        binned.concats = unique(binned.concats, by=c('cidi','binid')) 
        ##dt.concats = unique(dt.concats, by=c('cidi','binid'))
        dt.concats = binned.concats[, c('cidi','binid')]
        
        concats.gr = dt2gr(binned.concats)
        simplicial.complex = make_weighted_simplicial_complex_simpler(dt.concats)
        contact_matrix = cocount(concats.gr, bins = bins.gr, by = 'cid')

        pairwise = contact_matrix$dat
        pairwise[, dist := j - i]
        pairs.to.test = pairwise[dist > 1 & value > min.value]

     
        concat.decomps = simplicial.complex[[2]]
        simplicial.complex = simplicial.complex[[1]]
        iterative.registry.unlist.filterable = simplicial.complex[N==simplex.order, .(unlist(V1), num.concats, N, sidi), by='hashes']
    }

    ## else {
##         iterative.registry.unlist.filterable = simplicial.complex[N==simplex.order, .(unlist(V1), num.concats, N, sidi), by='hashes']
##         pairwise = contact_matrix$dat
##         pairwise[, dist := j - i]
##         pairs.to.test = pairwise[dist > 1 & value > min.value]
##     }
    
        

##     colnames(pairwise)[3] = 'pair.value'
##     pairwise = pairwise[i != j]
##     ##pairwise.unlist = melt(pairwise[, c('i','j','id','pair.value')], id.vars=c('id','pair.value'))
    

## ####this is the line graph algorithm again

    
###actually this is way more straightforward than line graph shit 

    pairs.to.test = pairwise[pair.value > min.value & dist > 1]
    str.identifiers = do.call(Map, c(f = c, pairs.to.test[, c('i','j')]))
    pairs.to.test$agg = str.identifiers
    pairs.to.test$pair.hashes = map(str.identifiers, unname) %>% paste0


    
    print('doing a big combn')
    iterative.registry.unlist.filterable$label = rep(c('i','j','k'), dim(iterative.registry.unlist.filterable)[[1]]/3)
    
    ##asdfasdf
    casting = data.table::dcast(iterative.registry.unlist.filterable, hashes ~ label, value.var='V1')
    subset.1 = merge(casting, pairs.to.test, by=c('i','j'))[, c('hashes','pair.hashes')]
    subset.2 = merge(casting, pairs.to.test, by.x=c('i','k'), by.y=c('i','j'))[, c('hashes','pair.hashes')]
    subset.3 = merge(casting, pairs.to.test, by.x=c('j','k'), by.y=c('i','j'))[, c('hashes','pair.hashes')]               ###okay this is kind of elegant


    combn.city = rbind(subset.1, subset.2, subset.3)



    
    ## print('doing a big combn')
    ## if(is.null(combn.city)){
    ##     iterative.registry.filt.min = iterative.registry.unlist.filterable[num.concats >= min.simplex]
    ##     combn.city = iterative.registry.filt.min[, .(combn(V1,2,simplify=FALSE), num.concats), by='hashes']
    ##     combn.city$pair.hashes = factor(combn.city$V1 %>% paste0)
    ## }
    ## ## unique.sidi = iterative.registry.unlist.filterable[num.concats > 1]$sidi %>% unique
    ## ## sidl = split(unique.sidi, ceiling(runif(length(unique.sidi))*numchunks)) ## randomly chop up cids
    ## ## setkey(iterative.registry.unlist.filterable, by='sidi')

    ## ## ##sidl = sidl[1:2]
    ## ## ##sidis = sidl[[2]]
    ## ## big.combn = pbmclapply(sidl, mc.cores=10, function(sidis) {
    ## ##     new.comb = iterative.registry.unlist.filterable[.(sidis), .(combn(V1,2,simplify=FALSE), num.concats), by='hashes'] ##this shit takes forever
    ## ##     new.comb$pair.hashes = factor(new.comb$V1 %>% paste0)
    ## ##     new.comb[, V1 := NULL]
    ## ##     return(new.comb)
    ## ## })

    ## ## combn.city = rbindlist(big.combn)
    
    ## ##pairs.to.test = pairwise
    ## str.identifiers = do.call(Map, c(f = c, pairs.to.test[, c('i','j')]))
    ## pairs.to.test$agg = str.identifiers
    ## pairs.to.test$pair.hashes = factor(map(str.identifiers, unname) %>% paste0)

    ## dt.pair.merge = merge(pairs.to.test, combn.city, by='pair.hashes')
    
    ## dt.pair.merge[, id := NULL]
    ## dt.pair.merge[, agg := NULL]
    ## dt.pair.merge[, V1 := NULL]

    ##rerun

    ##iterative.registry.unlist.filterable.2 = gm12878_registry[, .(unlist(V1), num.concats, N), by='hashes']

    dt.pair.merge = merge(pairs.to.test, combn.city, by='pair.hashes')
    
    dt.pair.merge[, id := NULL]
    dt.pair.merge[, agg := NULL]
    ##dt.pair.merge[, V1 := NULL]


    ##get third bin C

    ##dt = merge(dt, simplicial.complex[, c('hashes','sidi')], by='hashes')
    dt = merge(dt.pair.merge, iterative.registry.unlist.filterable, allow.cartesian=TRUE, by='hashes')

    dt[, dist := j - i]
    
    dt[, num.concats.y := NULL]
    ##colnames(dt)[7] = 'num.concats'
    dt = dt[V1 != i & V1 != j]



    
    ###browser()


    dt[, dist.i := abs(V1 - i)]
    dt[, dist.j := abs(V1 - j)]



    ###actually you could literally just tack this on at the end I think...

##     bins$binid = 1:length(bins)
##     bins.dt = gr2dt(bins)
##     setkey(bins.dt, 'binid')

##     dt$chr.i = bins.dt[dt$i]$seqnames
##     dt$chr.j = bins.dt[dt$j]$seqnames
##     dt$chr.V1 = bins.dt[dt$V1]$seqnames

## ###interchromosomal distance

##     dt[chr.i != chr.V1, dist.i := 2000]
##     dt[chr.j != chr.V1, dist.j := 2000]


##     dt[, chr.i := NULL]
##     dt[, chr.j := NULL]
##     dt[, chr.V1 := NULL]

####return.training = false: obliterate window,
    ##check every 3 way pair within bin range if they have a concatemer
    ##return only nonzero

    if(return.training==TRUE) {
        numchunks = 500


        pairs.to.test.sub = pairs.to.test[pair.value > min.value & i != j]
        
        
        ##window.size = length(bins.gr)
        ##}


        pairs.to.test.sub[, min.window := NULL]
        pairs.to.test.sub[, max.window := NULL]

        pairs.to.test.sub[, min.window.i := i - window.size]
        pairs.to.test.sub[, max.window.i := i + window.size]

        pairs.to.test.sub[, min.window.j := j - window.size]
        pairs.to.test.sub[, max.window.j := j + window.size]

        pairs.to.test.sub[min.window.i < 1, min.window.i := 1]
        pairs.to.test.sub[min.window.j < 1, min.window.j := 1]

        pairs.to.test.sub[max.window.i > length(bins.gr), max.window.i := length(bins.gr)]
        pairs.to.test.sub[max.window.j > length(bins.gr), max.window.j := length(bins.gr)]
        
        
        print('creating window bins')

        ##if(return.training==TRUE) {

        ###can vary sampling procedure

        allbins.yeet = 1:length(bins.gr)

        ##pairs.to.test.sub[, windowbins := list(c(min.window.i:max.window.i, min.window.j:max.window.j) %>% unique), by='pair.hashes']
        sampling = pairs.to.test.sub[, .(windowbins = list(c(min.window.i:max.window.i, min.window.j:max.window.j) %>% unique), agg, min.window.i, max.window.i, min.window.j, max.window.j, dist, i, j), by='pair.hashes']

        wide.sample = sampling[dist > window.size, .(samplebins = list(sample(c(1:min.window.i, max.window.i+1:min.window.j-1, max.window.j:length(bins.gr)), 100)) %>% unique), by='pair.hashes']
        normal.sample = sampling[dist <= window.size, .(samplebins = list(sample(c(1:min.window.i, max.window.j:length(bins.gr)), 100)) %>% unique), by='pair.hashes']

        all.samples = rbind(wide.sample, normal.sample)
        windowbins = merge(sampling, all.samples, by='pair.hashes', all.x=TRUE)
        ##sampling[, allbins := list(c(windowbins, samplebins)), by='pair.hashes']
        
        
        ##windowbins = pairs.to.test.sub[200:210, .(windowbins = list(c(min.window.i:max.window.i, sample(c(1:(min.window-1), (max.window+1):length(bins.gr)), 100)) %>% unique), pair.hashes, agg, i, j), by='pair.hashes']

        ##} else {
       ##     windowbins = pairs.to.test.sub[, .(windowbins = list(min.window:max.window), pair.hashes, agg, i, j), by='pair.hashes']
       ## }
        windowbins$pairid = 1:dim(windowbins)[[1]]
        upair = windowbins$pairid %>% unique
        upairl = split(upair, ceiling(runif(length(upair))*numchunks)) ## randomly chop up cids
        setkey(windowbins, by='pairid')

        pairids = upairl[[1]]
        ##browser()



#####THIS IS ONLY NECESSARY IF YOU CREATING A TRAINING SET
        
        big.decay = pbmclapply(upairl, mc.cores = 10, function(pairids) {
            dt.window.1 = unique(windowbins[.(pairids), .(as.integer(unlist(windowbins)), agg, i, j), by='pair.hashes'], by=c('pair.hashes','V1'))
            dt.window.2 = unique(windowbins[.(pairids), .(as.integer(unlist(samplebins)), agg, i, j), by='pair.hashes'], by=c('pair.hashes','V1'))
            dt.window = rbind(dt.window.1, dt.window.2)
            mergerious = merge(dt.window, dt[, c("hashes","pair.hashes","dist","num.concats","V1","N","sidi","dist.i","dist.j")], by=c('pair.hashes','V1'), all.x=TRUE)
            mergerious[is.na(num.concats), num.concats := 0]
            mergerious = mergerious[V1 != i & V1 != j]
            mergerious[, dist.i := abs(V1 - i)]
            mergerious[, dist.j := abs(V1 - j)]
            mergerious[, c('hashes','dist','N','sidi') := NULL]
            return(mergerious)
        })
        
        mergerious = rbindlist(big.decay)

        
    } else {
        mergerious = dt[, c('pair.hashes','V1','i','j','num.concats','dist.i','dist.j')]
        mergerious = merge(mergerious, pairs.to.test[, c('pair.hashes','agg')], by='pair.hashes')
    }
    
#####I think you jump here 

    ###need to merge here with all pairs, not just subset you are interrogating

    print('filling in pair gaps')

    mergerious = unique(mergerious, by=c('pair.hashes','V1'))
    dt.sub = mergerious[, .(sub.bin = unlist(agg), num.concats, i, j), by=c('pair.hashes','V1')]

    ##dt.sub = dt[1:100]
    ##dt.sub[, .(sub = unlist(agg), num.concats, value, pair.hashes), by=c('hashes','V1')]
    ###swap sub bin with V1 if sub.bin is less than V1

    dt.sub[sub.bin < V1, c('sub.bin','V1') := .(V1, sub.bin)]

    ##dt.yeet = dt.sub[1:100]
    ##dt.yeet[sub.bin < V1, c('sub.bin','V1') := .(V1, sub.bin)]

    ##colnames(dt.sub)[5] = 'pair.value'



    pairwise$agg = do.call(Map, c(f = c, pairwise[, c('i','j')]))
    pairwise$pair.hashes = factor(map(pairwise$agg, unname) %>% paste0)


    ##colnames(pairwise)[3] = 'pair.value'
    dt.sub = merge(dt.sub, pairwise[, c('i','j','pair.value')], by.x=c('V1','sub.bin'), by.y=c('i','j'), all.x=TRUE)
    dt.sub[is.na(pair.value), pair.value := 0]

    dt.sub = dt.sub[V1 != sub.bin]
    ##dt.sub[, sum.pairwise.contact := sum(value), by=c('pair.hashes','V1')]

    ##dt.sub[10000:10100][, .(sum.pairwise.contact = sum(value)), by=c('pair.hashes','V1')]


    ##dt.sub.backup = dt.sub %>% copy
    ##dt.sub = dt.sub.backup %>% copy

    ##dt.sub = merge(dt.sub, pairwise[, c('i','j','pair.hashes','pair.value')], by='pair.hashes', all.x=TRUE)
    dt.sub$V1.isinter = !(dt.sub[, V1 == i] | dt.sub[, V1 == j])
    dt.sub$sub.isinter = !(dt.sub[, sub.bin == i] | dt.sub[, sub.bin == j])


    dt.sub$binterrogate = 0

    ##dt.sub[1:100][V1.isinter==TRUE, binter
    ##dt.subbington = dt.sub[1:100]

    ##dt.subbington[, binterrogate := NULL]
    ##dt.subbington$binterrogate = 0


    dt.sub[V1.isinter==TRUE, binterrogate := V1]
    dt.sub[sub.isinter==TRUE, binterrogate := sub.bin]
    ##dt.yeet[!is.na(binterrogate)]

    dt.sub = dt.sub[binterrogate!=0]


    ##colnames(dt.sub)[6] = 'pair.value'
    
    dt.sub[, sum.pairwise.contacts := sum(pair.value), by=c('pair.hashes','binterrogate')]

    
    dt.sub[, dist.i := abs(binterrogate - i)]
    dt.sub[, dist.j := abs(binterrogate - j)]


    
    bins.gr$binid = 1:length(bins.gr)
    bins.dt = gr2dt(bins.gr)
    setkey(bins.dt, 'binid')

    dt.sub$chr.i = bins.dt[dt.sub$i]$seqnames
    dt.sub$chr.j = bins.dt[dt.sub$j]$seqnames
    dt.sub$chr.binterrogate = bins.dt[dt.sub$binterrogate]$seqnames

## ###interchromosomal distance

    dt.sub[chr.i != chr.binterrogate, dist.i := 2000]
    dt.sub[chr.j != chr.binterrogate, dist.j := 2000]


    dt.sub[, chr.i := NULL]
    dt.sub[, chr.j := NULL]
    dt.sub[, chr.V1 := NULL]


    
    dt.sub[, diff := j-i]
    dt.sub = dt.sub[diff > 1]

    
    dt.sub[, pair.value := pair.value + 1] ###for log covariate purposes

    print('calculating close and far pairwise contact values')
    
    dt.sub[binterrogate < i & j == sub.bin, value.b := pair.value] ##CLOSER BIN
    dt.sub[binterrogate < i & i == sub.bin, value.a := pair.value] ##FARTHER BIN

    dt.sub[binterrogate > j & i == V1, value.a := pair.value]
    dt.sub[binterrogate > j & j == V1, value.b := pair.value]

##dt.sub[binterrogate < i & j == sub.bin, value..b := value * abs(j - binterrogate)] ##CLOSER BIN
##dt.sub[binterrogate < i & i == sub.bin, value..a := value * abs(i - binterrogate)] ##FARTHER BIN

##dt.sub[binterrogate > j & i == V1, value.ratio.a := value * abs(i - binterrogate)] ##FARTHER BIN
##dt.sub[binterrogate > j & j == V1, value.ratio.b := value * abs(j - binterrogate)] ##CLOSER BIN

    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) < (j - binterrogate)) & i==V1, value.a := pair.value]
    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) < (j - binterrogate)) & j==sub.bin, value.b := pair.value]


    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) >= (j - binterrogate)) & j==sub.bin, value.a := pair.value]
    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) >= (j - binterrogate)) & i==V1, value.b := pair.value]





    

    dt.sub.unique = unique(dt.sub, by=c('pair.hashes','binterrogate'))
    ##dt.sub.unique[, value.ratio.a := NULL]
    ##dt.sub.unique[, value.ratio.b := NULL]





    ##a.ratio = dt.sub[!is.na(value.ratio.a), c('pair.hashes','value.ratio.a','binterrogate')]
    ##b.ratio = dt.sub[!is.na(value.ratio.b), c('pair.hashes','value.ratio.b','binterrogate')]

    ##dt.sub.unique = merge(dt.sub.unique, a.ratio, by=c('pair.hashes','binterrogate'))
    ##dt.sub.unique = merge(dt.sub.unique, b.ratio, by=c('pair.hashes','binterrogate'))




    dt.sub.unique[, value.a := NULL]
    dt.sub.unique[, value.b := NULL]

    a.value = dt.sub[!is.na(value.a), c('pair.hashes','value.a','binterrogate')]
    b.value = dt.sub[!is.na(value.b), c('pair.hashes','value.b','binterrogate')]

    dt.sub.unique = merge(dt.sub.unique, a.value, by=c('pair.hashes','binterrogate'))
    dt.sub.unique = merge(dt.sub.unique, b.value, by=c('pair.hashes','binterrogate'))

    dt.sub.unique[, dist.a := min(dist.i, dist.j), by=c('pair.hashes','binterrogate')]
    dt.sub.unique[, dist.b := max(dist.i, dist.j), by=c('pair.hashes','binterrogate')]

    ##dt.sub.agg = dt.sub[, .(value.ratio.b = value.ratio[1], value.ratio.b = value.ratio[1]), by=c('pair.hashes','binterrogate')]
    ##dt.sub.agg[, value.ratio.a := unlist(V1)[[1]]]
    ##dt.sub.agg[, value.ratio.b := unlist(V1)[2]]

    ##dt.sub.unique[, value.ratio.a := NULL]
    ##dt.sub.unique[, value.ratio.b := NULL]
    
    dt.small = dt.sub.unique[, c('pair.hashes','num.concats','pair.value','binterrogate','sum.pairwise.contacts','dist.a','dist.b','value.a','value.b')]
    dt.small = unique(dt.small, by=c('pair.hashes','binterrogate'))


    
    ##dt.small[dist.a != 1 & dist.b != 1] ###remove bins directly adjacent to 

    dt.small[, value.a.ratio := value.a / dist.a]
    dt.small[, value.b.ratio := value.b / dist.b]
    dt.small = dt.small[dist.a > 1 & dist.b > 1]


    ##saveRDS(dt.small, 'dt_small_25k.rds')
    ##saveRDS(simplicial.complex, 'weighted_simplicial_complexes/simplicial_complex_chr8_25k.rds')
    ##output = list(dt.small, combn.city, pairwise)
    ##return(list(dt.small, combn.city, pairwise))
    ##saveRDS(dt.small, 'dist_decay_windows_training_sampled_faraway.rds')
    ##saveRDS(dt.small, 'dist_decay_windows_training_test_chr8.rds')
    return(dt.small)
}


train_model = function(dt.small) {
    ###making a CHANGE: we are only going to train the model on non-zero values!! this makes the p-value way more useful
    dt.small.model = dt.small[num.concats>0, c('pair.hashes','num.concats','binterrogate','value.a.ratio','value.b.ratio')]
    
    covariates = c('value.a.ratio','value.b.ratio')
    fmstring = paste('num.concats ~', paste(paste0('log(', covariates, ')', collapse = ' + ')))
  ## fmstring = paste0(fmstring, " + ", "offset(log(width))")
    fm = formula(fmstring)

    ##need to think about how to train split 
    model = glm.nb(formula = fm, data = dt.small.model[, control = glm.control(maxit = 500)])
    return(model)
}


score_distance_decay = function(dt.small.model, model){
    dt.small.model$num.concats.pred = (predict(model, type = "response", newdata = dt.small.model))
    pval = dt.small.model[, pnbinom(num.concats -1, mu = num.concats.pred, size = model$theta, lower.tail = F)]
    pval.right = dt.small.model[, pnbinom(num.concats, mu = num.concats.pred, size = model$theta, lower.tail = F)]
    pval.right = ifelse(is.na(pval.right), 1, pval.right)
    pval = ifelse(is.na(pval), 1, pval)
    dt.small.model[, enrichment := num.concats / num.concats.pred]
    dt.small.model$pval = runif(nrow(dt.small.model), min = pval.right, max = pval)
    return(dt.small.model)
}


score_distance_decay_poisson = function(dt.small.model, model){
    dt.small.model$num.concats.pred = (predict(model, type = "response", newdata = dt.small.model))
    pval = ppois(dt.small.model$num.concats -1, lambda = dt.small.model$num.concats.pred, lower.tail = F)
    pval.right = ppois(dt.small.model$num.concats, lambda = dt.small.model$num.concats.pred, lower.tail = F)
    pval.right = ifelse(is.na(pval.right), 1, pval.right)
    pval = ifelse(is.na(pval), 1, pval)
    dt.small.model$pval = runif(nrow(dt.small.model), min = pval.right, max = pval)
    
    dt.small.model[, enrichment := num.concats / num.concats.pred]
    return(dt.small.model)
}

    
make_barplot = function(synergy_results, filename='plot.png', fdr.thresh=0.1) {
    ##synergy_results = synergy_results[cardinality > 3]
    barplot_dat = synergy_results[, total_annot := .N, by=c('annotation')]


    synergy_results[, sig := NULL]
    synergy_results$sig = FALSE
    synergy_results[fdr<fdr.thresh, sig := TRUE]

    barplot_dat = synergy_results[, total := sum(sig), by=c('annotation')]
    barplot_dat[, proportion := total / total_annot]

    ##barplot_dat = barplot_dat[!is.na(sig)]
    ##barplot_dat = unique(barplot_dat, by='proportion')[,c('annotation','sig','proportion', 'total','total_annot')]
    ##barplot_dat = barplot_dat[barplot_dat$sig == 'yes']


    barplot_dat = unique(barplot_dat, by=c('proportion','annotation'))

    binconfs = binconf(barplot_dat$total, barplot_dat$total_annot)
    colnames(binconfs) = c('proportion', 'Lower', 'Upper')
    binconfs = unique(binconfs, by='proportion')

    barplot_dat = unique(barplot_dat, by=c('proportion','annotation'))


    barplot_dat = merge(barplot_dat, binconfs, by='proportion')



    ##browser()

    annots = barplot_dat$annotation %>% unique
    barplot_dat$annotation = factor(barplot_dat$annotation,
                                    levels = annots)

    barplot_dat$annotation = as.character(barplot_dat$annotation)
                                        # Most basic error bar

    ##browser()
    ##setorder(barplot_dat, -annotation)
    ##filename = 'barplot_comp_largerthan3.png'
    barplot_dat = unique(barplot_dat, by=c('proportion','annotation'))

    if(length(unique(barplot_dat$annotation)) > 3) {
        colors = c('blue','orange','red','yellow')
    } else {
        colors = c('blue','orange','red')
    }

    
    max.y = max(max(barplot_dat$Upper), 1)
    ppng(plot(ggplot(barplot_dat[,c('annotation','proportion','Lower','Upper')]) +
              geom_bar(aes(x=annotation, y=proportion), stat="identity", fill=colors, alpha=0.7) +
              geom_errorbar( aes(x=annotation, ymin=Lower, ymax=Upper), width=0.4) + ylim(0,max.y) +
              theme(text = element_text(size=30))), filename)
}




prepare_distance_decay_inputs = function(concatemers, bins, prebinned=FALSE) {

    if(prebinned==FALSE) {
        binned.concats = bin_concatemers(concatemers, bins = bins, max.slice = 1e6, mc.cores=5, verbose=TRUE)
    } else {
        binned.concats = concatemers
    }
    
    concats.gr = unique(binned.concats, by=c('cidi','binid')) %>% dt2gr
    
    ##dt.concats = unique(dt.concats, by=c('cidi','binid'))
    dt.concats = unique(binned.concats[, c('cidi','binid')], by=c('cidi','binid'))

    simplicial.complex = make_weighted_simplicial_complex_simpler(dt.concats)

    iterative.registry = simplicial.complex[[1]]
    simplex.order = 3
    iterative.registry.unlist.filterable = iterative.registry[N==simplex.order, .(unlist(V1), num.concats, N, sidi), by='hashes']

    contact_matrix_unique = cocount(concats.gr, bins = bins, by = 'cid')

    pairwise = contact_matrix_unique$dat

    pairwise[, dist := j - i]
    pairwise$agg = do.call(Map, c(f = c, pairwise[, c('i','j')]))
    pairwise$pair.hashes = map(pairwise$agg, unname) %>% paste0
    colnames(pairwise)[colnames(pairwise) == 'value'] = 'pair.value'
    
    return(list(iterative.registry, iterative.registry.unlist.filterable, contact_matrix_unique, pairwise, binned.concats))
}
####HERE WE HAVE THE INPUTS TO AN ANNOTATE DISTANCE DECAY FUNCTION



annotate_and_score_distance_decay = function(iterative.registry, iterative.registry.unlist.filterable, contact_matrix_unique, pairwise, concatemers, bins, model=NULL, pair.thresh=50) {

    if(is.null(model)){
        dist.decay.train = annotate_distance_decay(concatemers,
                                                   bins,
                                                   min.value = pair.thresh,
                                                   window.size=50,
                                                   simplicial.complex=iterative.registry,
                                                   iterative.registry.unlist.filterable = iterative.registry.unlist.filterable,
                                                   pairwise=pairwise,
                                                   return.training=TRUE)
        print('train done')
    }
    
    dist.decay.test = annotate_distance_decay(concatemers,
                                              bins,
                                              min.value = pair.thresh,
                                              window.size=50,                                                                                                        simplicial.complex=iterative.registry,
                                              iterative.registry.unlist.filterable = iterative.registry.unlist.filterable,
                                              pairwise = pairwise,
                                              return.training=FALSE)

    print('test done')
    
###comm detect



    ####total contacts within window
    if(is.null(model)) {
        total.concats = dist.decay.train[(dist.a < 50 | dist.b < 50), .(total.concats = sum(num.concats)), by=c('pair.hashes')]

        dist.decay.train = merge(dist.decay.train, total.concats, by='pair.hashes')
        dist.decay.test = merge(dist.decay.test, total.concats, by='pair.hashes')
        ##total.concats = dist.decay.test[(dist.a < 50 | dist.b < 50), .(total.concats = sum(num.concats)), by=c('pair.hashes')]
        
###GLM town

        print('training model')
        covariates = c('value.a.ratio','value.b.ratio')
        fmstring = paste('num.concats ~', paste(paste0('log(', covariates, ')', collapse = ' + ')))
        ##fmstring = paste0(fmstring, " + ", "offset(log(total.concats))") this sometimes does 
        fm = formula(fmstring)


        ##browser()
        num.to.sample = 500000
        if(num.to.sample > dim(dist.decay.train)[[1]]) {
            num.to.sample = dim(dist.decay.train)[[1]]
        }
            
        model = glm.nb(formula = fm, data=dist.decay.train[sample(.N, num.to.sample)], control=glm.control(maxit=500))
    }
    ##dist.decay.train = score_distance_decay(dist.decay.train, model)
    ##browser()
    dist.decay.test = score_distance_decay(dist.decay.test, model)
    dist.decay.test[, relative.risk := log2(abs(num.concats - num.concats.pred))]

    print('test scored')
    all.pairs.tested = big.combn(iterative.registry.unlist.filterable, contact_matrix_unique, min.value=pair.thresh)

    ##debug(create_bin_pair_network)


    dist.decay.test.nozero = merge(dist.decay.test, all.pairs.tested[, c('pair.hashes','binterrogate','hashes')], by=c('pair.hashes','binterrogate'))

    ##dist.decay.test.nozero[, relative.risk := log2(abs(num.concats - num.concats.pred))]

    ##undebug(create_bin_pair_network)
    print('creating bin pair network')
    bin_pair_network = create_bin_pair_network(all.pairs.tested, dist.decay.test.nozero, pairwise)

    
    return(list(bin_pair_network, dist.decay.test, model, all.pairs.tested))
}






derive_binsets_from_network = function(G.kant, pairwise, binned.concats, bins, rr.thresh = 0, dist.decay.test=NULL, all.pairs.tested=NULL, num.members=20, pairwise.trimmed=pairwise, expansion.cutoff = 0.5, pval.thresh = 0.05, fdr.thresh=0.1) {
##num.members=2
    ##G.kant = bin.pair.network
    
    pairwise.trimmed$agg = do.call(Map, c(f = c, pairwise.trimmed[, c('i','j')]))
    G.kant.10 = subgraph.edges(G.kant, E(G.kant)[E(G.kant)$weight > 0])
    G.kant.undir = as.undirected(G.kant.10)
    cl.l = cluster_leiden(G.kant.undir, objective_function='modularity')
    cl = cl.l$membership
    dt.membership = data.table(cluster = cl)
    dt.membership$pair = V(G.kant.10)$name
    ##dt.membership$classification = V(G.kant.10)$classification
    dt.membership[, num.memb := .N, by='cluster']

    ##dt.membership = merge(dt.membership, dt.G, by.x='pair', by.y='pair.hashes')





    ## asdf = pbmclapply(rr.threshes, mc.cores = 5, function(thresh) {
    ##     G.sub = subgraph.edges(G.kant.rr, E(G.kant.rr)[E(G.kant.rr)$weight > thresh])
    ##     G.kant.undir = as.undirected(G.sub)
    ##     cl.l = cluster_leiden(G.kant.undir, objective_function='modularity')
    ##     cl = cl.l$membership
    ##     dt.membership = data.table(cluster = cl)
    ##     dt.membership$pair = V(G.sub)$name
    ##     dt.membership$classification = V(G.sub)$classification
    ##     dt.membership[, num.memb := .N, by='cluster']
    ##     dt.membership$thresh = thresh
    ##     dt.membership =  merge(pairwise[, c('pair.hashes','agg')], dt.membership, by.y='pair', by.x='pair.hashes')
    ##     dt.membership = dt.membership[, .(individual.bin = unlist(agg), classification, thresh, cluster), by='pair.hashes']
    ##     pagerank = page_rank(G.sub)
    ##     dt.page = data.table(pair = V(G.sub)$name)
    ##     dt.page$pagerank = pagerank[[1]]
    ##     dt.membership = merge(dt.page, dt.membership, by.x='pair', by.y='pair.hashes', all.y=TRUE)
    ##     dt.membership[is.na(pagerank), pagerank := 0]
    ##     return(dt.membership)
    ## })

    unique.clusters = dt.membership[num.memb > num.members]$cluster %>% unique

    ## binsets = pbmclapply(unique.clusters, mc.cores = 5, function(cluster.id) {
    ##     G.sub = subgraph(G.kant.10, dt.membership[cluster==cluster.id]$pair)
    ##     sub.eig = eigen_centrality(G.sub, directed=TRUE)##, cutoff=-1)##, options)
    ##     dt.subgraph = data.table(pairs = V(G.sub)$name)
    ##     dt.subgraph$eigen = sub.eig[[1]]
    ##     setorder(dt.subgraph, -eigen)
    ##     dt.eig.sub = merge(dt.subgraph, pairwise[, c('pair.hashes','agg')], by.x='pairs', by.y='pair.hashes')
    ##                                     #dt.eig.sub = dt.eig.sub
    ##     dt.eig.sub = dt.eig.sub[, .(unlist(agg), eigen), by='pairs']
    ##     dt.eig.sub = dt.eig.sub[, .(.N, eigen), by='V1']
    ##     unique(dt.eig.sub, by='V1')
    ##     dt.final = dt.eig.sub[, .(sum.eig = sum(eigen), count=N), by='V1'] %>% unique(by='V1')
    ##     ##dt.final[, .(sum.eig / count)]
    ##     setorder(dt.final, -sum.eig)
    ##     mean.eig = dt.final$sum.eig %>% mean
    ##     binset = dt.final[sum.eig > mean.eig]
    ##     binset$cluster = cluster.id
    ##     return(binset)
    ## })


    dt.eig.sub = pbmclapply(unique.clusters, mc.cores = 5, function(cluster.id) {
        G.sub = subgraph(G.kant.10, dt.membership[cluster==cluster.id]$pair)
        sub.eig = eigen_centrality(G.sub, directed=TRUE)##, cutoff=-1)##, options)
        dt.subgraph = data.table(pairs = V(G.sub)$name)
        dt.subgraph$eigen = sub.eig[[1]]
        if(class(pairwise.trimmed$pair.hashes)=='integer'){
            dt.subgraph$pairs = strtoi(dt.subgraph$pairs)
        }
        setorder(dt.subgraph, -eigen)
        dt.eig.sub = merge(dt.subgraph, pairwise.trimmed[, c('pair.hashes','agg')], by.x='pairs', by.y='pair.hashes')
        dt.eig.sub$cluster = cluster.id
        return(dt.eig.sub)
    })

    
    if(length(unique.clusters) == 1) {
        dt.eig.sub.agg = rbindlist(dt.eig.sub[[1]])
        dt.eig.sub.agg$eigen = 1
    } else {
        dt.eig.sub.agg = rbindlist(dt.eig.sub)
    }
        
    dt.eig.sub.agg[, max.eig := max(eigen), by='cluster']
    seeds = dt.eig.sub.agg[eigen == max.eig]

    colnames(seeds)[1] = 'pair.hashes'
    seeds = unique(seeds[max.eig>0], by='cluster')

    print(seeds)
    ##merge(seeds, dist.decay.test, by='pair.hashes')
##OA
    ##colnames(all.pairwise)[4] = 'pair.hashes'
    ##dist.decay.test = dist.decay[pval<0.05]
    binsets.gr = expand_seeds(seeds, dist.decay.test, pairwise, all.pairs.tested=NULL, bins, fdr.thresh=fdr.thresh, pval.thresh=pval.thresh, expansion.cutoff=expansion.cutoff)
    binsets.dt = gr2dt(binsets.gr)
    binsets.dt$chid = binsets.dt$bid

    
    ## bins = dist.decay.test[pair.hashes=='c(2547, 2552)']
    ## setorder(asdf, -relative.risk)
    
    
    
    ## binsets.dt = rbindlist(binsets)
    ## colnames(binsets.dt)[1] = 'binid'
    ## colnames(binsets.dt)[4] = 'bid'
    ## binsets.gr = bins[binsets.dt$binid]
    ## ##binsets.gr$bid = binsets.dt$bid
    ## ##binsets.gr$binid = binsets.dt$binid

###For filtering highly contacted bins
##This is quite important    
    binned.concats[, bincount := .N, by='binid']
    thresh = quantile(unique(binned.concats, by='binid')$bincount, 0.9999) ###Remove very highly mapped bins
    standard.dev = sd(unique(binned.concats, by='binid')$bincount)
    binned.concats = binned.concats[bincount < (thresh + standard.dev)]
    
    
    binned.concats[, bid := NULL]
    binned.concats[, chid := NULL]
    concat.merge = merge(binned.concats, binsets.dt[, c('binid','bid')], by='binid', allow.cartesian=TRUE)
    concat.count = concat.merge[, .N, by=c('cid','bid')]

##    binsets.dt
    concat.contacts = concat.count[N>2]
##    tpe1.dt = gr2dt(tpe1_chr8)
    ##concat.contacts
    
    concat.contacts$cid = concat.contacts$cid %>% as.integer
##    tpe1.chr8.binned$cid = tpe1.chr8.binned$cid %>% as.integer
    binned.concats$cid = binned.concats$cid %>% as.integer
    binned.concats = unique(binned.concats, by=c('cid','binid')) ##we are so BACK

    ###browser()
    contacted.concatemers = merge(binned.concats, concat.contacts, by='cid', all.y=TRUE, allow.cartesian=TRUE)
    contacted.concatemers = dt2gr(contacted.concatemers)
    contacted.concatemers$chid = contacted.concatemers$bid

    
    new.binsets = rebin_community(contacted.concatemers, unique(concat.contacts$bid), resolution=10000, rebin_thresh=0.85)

    
    
    new.binsets.dt = gr2dt(new.binsets)
    new.binsets.dt[, numbins := .N, by='chid']
    new.binsets.dt = new.binsets.dt[width<200000 & numbins<500] ##I will hardcode this as limit: something has gone horribly wrong if you reach this

    new.binsets = dt2gr(new.binsets.dt[numbins > 2])
    chrom = Chromunity(binsets=dt2gr(new.binsets.dt), concatemers=contacted.concatemers, meta=data.table())
    return(chrom)

}


create_bin_pair_network = function(all.pairs.tested, dt.small.model.sub, pairwise){
    
    subsets.fillin = merge(all.pairs.tested[, c('hashes','pair.hashes')], dt.small.model.sub[pval < 0.05], by=c('pair.hashes','hashes'), all.x=TRUE)
    ##subsets.fillin[!is.na(enrichment)]
    ##subsets.fillin[!is.na(enrichment)]
    subsets.fillin[, pair.sig := is.na(enrichment)]
    subsets.fillin[, num.sig.pairs := sum(!pair.sig), by='hashes']
    subsets.fillin = subsets.fillin[num.sig.pairs > 0]

    ##subsets.fillin[, relative.risk := log2(abs(num.concats - num.concats.pred))]
    dt.net = merge(subsets.fillin[!is.na(enrichment), c('pair.hashes','hashes','enrichment','num.concats','relative.risk')], subsets.fillin[is.na(enrichment), c('pair.hashes','hashes')], by='hashes')
    dt.net.2 = merge(subsets.fillin[!is.na(enrichment), c('pair.hashes','hashes','enrichment','num.concats','relative.risk')], subsets.fillin[!is.na(enrichment), c('pair.hashes','hashes')], by='hashes', allow.cartesian=TRUE)[pair.hashes.x != pair.hashes.y]
    dt.net.agg = rbind(dt.net, dt.net.2)


    ##pairwise = contact_matrix$dat
    ##pairwise[, dist := j - i]
    ##pairs.to.test = pairwise[dist > 1 & value > min.value]

    ##pairwise$agg = do.call(Map, c(f = c, pairs.to.test[, c('i','j')]))
    ##pairs.to.test$pair.hashes = map(pairs.to.test$agg, unname) %>% paste0



    colnames(pairwise)[3] = 'pair.value'
    dt.net.agg = merge(pairwise[, c('pair.hashes','pair.value')], dt.net.agg, by.x='pair.hashes',by.y='pair.hashes.x')
    colnames(dt.net.agg)[1:2] = c('pair.hashes.x','pair.value.x')
    dt.net.agg = merge(pairwise[, c('pair.hashes','pair.value')], dt.net.agg, by.x='pair.hashes',by.y='pair.hashes.y')
    colnames(dt.net.agg)[1:2] = c('pair.hashes.y','pair.value.y')
    ##browser()
    G.kant.thresh = graph_from_edgelist(as.matrix(dt.net.agg[,c('pair.hashes.x','pair.hashes.y')]), directed=TRUE)
    E(G.kant.thresh)$weight = dt.net.agg$relative.risk


    ##dt.sub = pairwise[, c('pair.hashes','i','j')]
    ##nodelist = V(G.kant.thresh)$name %>% as.data.table
    ##colnames(nodelist) = 'pair.hashes'
    ##dt.sub.pair = unique(dt.sub[, c('pair.hashes','i','j')], by='pair.hashes')



    ##nodelist = merge(nodelist, dt.sub.pair, by='pair.hashes')
    ##nodelist = merge_nodelist(nodelist, prom.bins.dt.binids, varname='prom')
    ##nodelist = merge_nodelist(nodelist, enh.bins.dt.binids, varname='enh')
    ##nodelist = merge_nodelist(nodelist, sup.bins.dt.binids, varname='sup')



    ##nodelist[, classification := pair_classification(i.sup, j.sup, i.prom, j.prom, i.enh, j.enh), by='pair.hashes']

    ##dt.net.agg = merge(nodelist[, c('pair.hashes','classification')], dt.net.agg, by.x='pair.hashes', by.y='pair.hashes.y')
    ##colnames(dt.net.agg)[1] = 'pair.hashes.y'
    ##G.kant = set_vertex_attr(G.kant.thresh, 'classification', index=nodelist$pair.hashes, nodelist$classification)    
    return(list(G.kant.thresh, dt.net.agg))
}


dist_decay_binsets = function(concatemers, bins, model=NULL, rr.thresh=2.5) {
    output = prepare_distance_decay_inputs(concatemers, bins)
    
    iterative.registry = output[[1]]
    iterative.registry.unlist.filterable = output[[2]]
    contact_matrix_unique = output[[3]]
    pairwise = output[[4]]
    binned.concats = output[[5]]


    annotate.output = annotate_and_score_distance_decay(iterative.registry, iterative.registry.unlist.filterable, contact_matrix_unique, pairwise, concatemers, bins)
    bin.pair.network = annotate.output[[1]]
    dist.decay.test = annotate.output[[2]]
    model = annotate.output[[3]]
    all.pairs.tested = annotate.output[[4]]

    G.kant = bin.pair.network[[1]]

    print('creating binsets from network')
    chrom = derive_binsets_from_network(G.kant, pairwise, binned.concats, bins, rr.thresh=rr.thresh, dist.decay.test=dist.decay.test, all.pairs.tested=all.pairs.tested)

    ##browser()
    return(list(chrom$binsets, bin.pair.network, dist.decay.test, model, chrom))
}

big.combn = function(iterative.registry.unlist.filterable, contact_matrix, min.value=50) {

    pairwise = contact_matrix$dat
    pairwise[, dist := j - i]
    pairs.to.test = pairwise[dist > 1 & value > min.value]

    pairs.to.test$agg = do.call(Map, c(f = c, pairs.to.test[, c('i','j')]))
    pairs.to.test$pair.hashes = map(pairs.to.test$agg, unname) %>% paste0
    
    iterative.registry.unlist.filterable$label = rep(c('i','j','k'), dim(iterative.registry.unlist.filterable)[[1]]/3)    
    ##might genuinely be some of the most elegant code I've written lmao
    casting = data.table::dcast(iterative.registry.unlist.filterable, hashes ~ label, value.var='V1')
    subset.1 = merge(casting, pairs.to.test, by=c('i','j'))[, c('hashes','pair.hashes','k')]
    subset.2 = merge(casting, pairs.to.test, by.x=c('i','k'), by.y=c('i','j'))[, c('hashes','pair.hashes','j')]
    subset.3 = merge(casting, pairs.to.test, by.x=c('j','k'), by.y=c('i','j'))[, c('hashes','pair.hashes','i')]               ###okay this is kind of elegant
    colnames(subset.1)[3] = 'binterrogate'
    colnames(subset.2)[3] = 'binterrogate'
    colnames(subset.3)[3] = 'binterrogate'
    combn.city = rbind(subset.1, subset.2, subset.3)
    return(combn.city)
}








######PIPELINE






generate_synergies = function(concatemers.input, master.directory.input, params, picture.folder='analysis', all.chr=NULL, mc.cores=10, mode='sliding', rr.thresh=3) {
    if (is.null(all.chr)) {
        all.chr = c(paste0("chr", c(as.character(1:22), "X")))
    }
    for(chr in all.chr){
         print(master.directory.input)
         chrom = generate_binsets(concatemers.input %Q% (seqnames == chr), cell.line.dir = master.directory.input, params=params, chromosome.input = chr, mc.cores=mc.cores, mode=mode, rr.thresh = rr.thresh)

         ##for reproducibility save leave out concatemer set
         saveRDS(chrom[[3]], paste0(chrom[[2]], '/leave_out_concatemers.rds'))
         syn = call_synergies(chrom[[1]], leave_out_concatemers=chrom[[3]], chromosome.input=chr, cell.line.dir = master.directory.input, condition.dir=chrom[[2]], params=params, picture.folder = paste0(picture.folder,chr,'/'), mode=mode)
     }
     return(list(chrom, syn))    
}





generate_binsets = function(concatemers, bins=NULL, condition.dir=NULL, params = NULL, cell.line.dir=NULL, mode='sliding', chromosome.input='chr8', mc.cores=10, rr.thresh=3) {
    ##cell.line.dir = paste0(getwd(), cell.line.dir)

    if(is.null(params)){
        params = data.table(k.knn=25, k.min=5, resolution=1e4, seed=198)
    }
    ##browser()
    if(is.null(condition.dir)){
        condition.dir = paste0(cell.line.dir, mode, '_', chromosome.input, '_knn', params$k.knn, '_kmin', params$k.min, '_resolution', params$resolution)
    }

    dir.create(condition.dir)
    saveRDS(params, file.path(condition.dir, 'params.rds'))
    unique_read_idx = concatemers$read_idx %>% unique

    ## if(!file.exists(file.path(cell.line.dir,'/training_concatemers.rds'))){
    ##     ## training
    ##     training_idx = sample(unique_read_idx, length(unique_read_idx)/2)
    ##     training_concatemers = concatemers %Q% (read_idx %in% training_idx)

    ##     ## Set diff training and rest
    ##     leave_out_idx = sample(setdiff(unique_read_idx, training_idx), length(unique_read_idx)/2)
    ##     leave_out_concatemers = concatemers %Q% (read_idx %in% leave_out_idx)

        
    ##     saveRDS(leave_out_concatemers, file.path(cell.line.dir, 'leave_out_concatemers.rds'))  #####
    ##     saveRDS(training_concatemers, file.path(cell.line.dir, 'training_concatemers.rds'))  #####
    ## } else {
    ##     training_concatemers = readRDS(file.path(cell.line.dir, 'training_concatemers.rds'))
    ##     leave_out_concatemers = readRDS(file.path(cell.line.dir, 'leave_out_concatemers.rds'))
    ## }


    ##training_concatemers = concatemers %Q% (read_idx %in% training_idx)

        
    
## 10M for testing    
###I am going to subsample myself this is annoying addy sorry

    
    if(mode=='sliding'){
        ###insert loop
        ##browser()
        chrom = sliding_window_chromunity(concatemers = concatemers,
                                                   take_sub_sample = TRUE,
                                                   chr  = chromosome.input,
                                                   subsample.frac = 0.5,
                                                   k.knn = params$k.knn, k.min = params$k.min,
                                                   resolution = params$resolution,
                                                   window.size = 2e6, mc.cores = 5, min.support = 3, seed=params$seed)

        
        this.chr=chromosome.input
        this.chrom = gr2dt(chrom$concatemers)

        this.sub.parq = concatemers %Q% (seqnames %in% c(this.chr))
        this_gr_testing = dt2gr(gr2dt(this.sub.parq)[!read_idx %in% unique(this.chrom[support > 1]$read_idx)])
##leave_out_concatemers = (colo.concats %Q% (seqnames=='chr10')) %Q% !(cid %in% chrom.25res$concatemers$cid)


    } else if (mode=='re') {
        set.seed(params$seed)
        training_idx = sample(unique_read_idx, length(unique_read_idx)/2)
        training_concatemers = concatemers %Q% (read_idx %in% training_idx)

    ##     ## Set diff training and rest
        leave_out_idx = sample(setdiff(unique_read_idx, training_idx), length(unique_read_idx)/2)
        this_gr_testing = concatemers %Q% (read_idx %in% leave_out_idx)

        chrom = re_chromunity(training_concatemers, k.knn = params$k.knn, k.min = params$k.min, windows = bins, resolution = NULL, piecewise = FALSE, shave = FALSE, bthresh = 5, cthresh = 3)


    } else if (mode=='dist_decay') {
        set.seed(params$seed)
        training_idx = sample(unique_read_idx, length(unique_read_idx)/2)
        training_concatemers = concatemers %Q% (read_idx %in% training_idx)

    ##     ## Set diff training and rest
        leave_out_idx = sample(setdiff(unique_read_idx, training_idx), length(unique_read_idx)/2)
        this_gr_testing = concatemers %Q% (read_idx %in% leave_out_idx)
        bins = gr.tile(hg_seqlengths(genome="BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), params$resolution) %Q% (seqnames==chromosome.input)
        dist_decay_gr = dist_decay_binsets(training_concatemers, bins)
        chrom = dist_decay_gr[[5]]
    }
    
    saveRDS(chrom, file.path(condition.dir, '/chromunity_results.rds'))
    return(list(chrom, condition.dir, this_gr_testing))
}
##oh i messed that up slightly 





call_synergies = function(res, leave_out_concatemers, cell.line.dir, condition.dir, chromosome.input, picture.folder, params, mc.cores=10, mode=NULL) {
    
    ##leave_out_concatemers = readRDS(paste0(cell.line.dir, 'leave_out_concatemers.rds'))
    output.folder = paste0(condition.dir, '/synergy_outputs/')
    dir.create(output.folder)

    synergies = evaluate_synergy(res, leave_out_concatemers, chromosome=chromosome.input, folder=output.folder, resolution=1e4, mc.cores=mc.cores)
    ##auto generate these

    
    make_barplot(synergies, paste0(picture.folder, paste0('synergy_', chromosome.input, '.png')))

    synergies[, log.est := log2(estimate)]
    volc = .make.volcano(synergies, yl=c(0, 10))
    ppdf(plot(volc$chromunity), paste0(picture.folder, 'chromunity.pdf'))
    ppdf(plot(volc$random), paste0(picture.folder, 'random.pdf'))
    ppdf(plot(volc$shuffled), paste0(picture.folder, 'shuffled.pdf'))
    return(synergies)
}



analyze_concats = function(new.concats.dt) {
     ##new.concats.dt = gr2dt(new.concats)
     contact.order.dist = data.table()
     new.concats.dt[, bincount := .N, by='read_idx']
     unique.concats = unique(new.concats.dt, by='read_idx')

    ##unique.concats = unique.rerun
    
    #number of chromosomes a concatemer hits    
     new.concats.dt[, num.chr := length(unique(seqnames)), by='read_idx']

    ##number of times a concatemer hits the given chromosome
     new.concats.dt[, chr.count := .N, by=c('read_idx','seqnames')]

    ##unique.chr.concats
     unique.chr.concats = unique(new.concats.dt, by=c('read_idx','seqnames'))
     unique.chr.concats[, intra.chr.contacts := choose(chr.count, 2)]

     ##unique.concats[, contact.count := choose(bincount, 2), by='read_name']
     unique.chr.concats[, contacts.count := choose(bincount, 2), by='read_idx']
     unique.concats = unique(unique.chr.concats, by=c('read_idx'))

     intra.vs.inter = unique.chr.concats[, .(intra.contacts = sum(intra.chr.contacts), contacts.count), by='read_idx']
     intra.vs.inter = unique(intra.vs.inter, by='read_idx')

     total.inter.contacts = sum(intra.vs.inter$contacts.count) - sum(intra.vs.inter$intra.contacts)
     intra.vs.inter = unique(intra.vs.inter, by='read_idx')

     percent.cis = sum(intra.vs.inter$intra.contacts) / sum(intra.vs.inter$contacts.count)
     contact.order.dist$percent.cis = percent.cis

     ##contact count
     contact.count = sum(unique.concats$contacts.count)
     contact.order.dist$contact.count = contact.count

     ##contacts per gb unique.c
     ##contacts.per.gb = contact.count / ((total.gb) / (10^9))
     ##pairs.new.row$contacts.per.gb = contacts.per.gb

     ##contacts.per.gb = (pairs.table[1]$contact.count) / (pairs.table[1]$total_bases / (10^9))
     ##contact.order.dist$contacts.per.gb = contacts.per.gb
     ##pairs.table[1]$contacts.per.gb = contacts.per.gb

     ###contact order distribution
     
     contact.order.dist$lessone = sum(unique.concats$bincount == 1)

     total.reads = dim(unique.concats)[[1]]

     contact.order.dist$multiway = sum(unique.concats$bincount > 1)
     contact.order.dist$two = sum(unique.concats$bincount == 2) / total.reads
     contact.order.dist$three = sum(unique.concats$bincount == 3)  / total.reads
     contact.order.dist$four = sum(unique.concats$bincount == 4) / total.reads
     contact.order.dist$fivetosix = sum((unique.concats$bincount == 5) | (unique.concats$bincount == 6)) / total.reads
     contact.order.dist$seventoeleven = sum((unique.concats$bincount >= 7) & (unique.concats$bincount <= 11)) / total.reads
     contact.order.dist$twelvetotwentyone = sum((unique.concats$bincount >= 12) & (unique.concats$bincount <= 21)) / total.reads
     contact.order.dist$twentytwotofortynine = sum((unique.concats$bincount >= 22) & (unique.concats$bincount <= 49)) / total.reads
     contact.order.dist$greaterthan50 = sum(unique.concats$bincount >= 50) / total.reads
     return(contact.order.dist)
}




expand_seeds = function(seeds, dist.decay.test, pairwise, all.pairs.tested, bins, expansion.cutoff = 0.5, pval.thresh = 0.05, fdr.thresh=0.1) {
    seeds.archive = seeds %>% copy

    seeds = unique(seeds, by=c('cluster'))
    seeds.loop = seeds[, c('pair.hashes','cluster')]
    
    
    memb.dt = seeds[, .(unlist(agg)), by='cluster']
    memb.dt$memb = TRUE
   
    i=1
    bins=gr.stripstrand(bins)
    if(class(dist.decay.test$pair.hashes)=='integer'){
        seeds.loop$pair.hashes = strtoi(seeds.loop$pair.hashes)
        all.pairs.tested$pair.hashes = strtoi(all.pairs.tested$pair.hashes)
    }

    pairwise.2 = pairwise %>% copy
    pairwise.2$i = pairwise$j
    pairwise.2$j = pairwise$i
    pairwise.sym = rbind(pairwise, pairwise.2)
    ##colnames(pairwise.sym)[colnames(pairwise.sym)=='id'] = 'pair.hashes'
    
    for(i in 1:50){
        if(class(dist.decay.test$pair.hashes)=='integer'){
            seeds.loop$pair.hashes = strtoi(seeds.loop$pair.hashes)
        }

        seed.decays = merge(dist.decay.test, seeds.loop, by.x='pair.hashes', by.y='pair.hashes', allow.cartesian=TRUE)
        ##seed.decays = dist.decay.test[pair.hashes %in% seeds.loop$pair.hashes]

        if(i==1) {
            seed.decays = seed.decays[fdr<fdr.thresh]
        }

        ##seed.decays = merge(seed.decays, seeds.loop[, c('pair.hashes','cluster')], by='pair.hashes', allow.cartesian=TRUE)
        seed.decays[pval < pval.thresh, sum.relative.risk := sum(relative.risk), by=c('binterrogate','cluster')]
        seed.decays[, num.concats := sum(num.concats), by=c('binterrogate','cluster')]
        seed.decays[, num.concats.pred := sum(num.concats.pred), by=c('binterrogate','cluster')]
        
        ##seed.decays = merge(seed.decays, all.pairs.tested, by=c('pair.hashes','binterrogate'))

        ####if we want to be verbose with it

        ##to.drop = merge(memb.dt, seed.decays, by.x=c('cluster','V1'), by.y=c('cluster','binterrogate'), all.y=TRUE, allow.cartesian=TRUE)
        
###also: DROP ANY BIN THAT IS DIRECTLY ADJACENT TO ANOTHER MEMBER BIN
        adj.bins = memb.dt[, .(adj.bins = list((V1-1):(V1+1))), by=c('cluster','V1')]
        adj.bins[, V1 := NULL]
        adj.bins = adj.bins[, .(unlist(adj.bins)), by='cluster']
        adj.bins = unique(adj.bins, by=c('cluster','V1'))
        adj.bins$memb = TRUE
        
        to.drop = merge(adj.bins, seed.decays, by.x=c('cluster','V1'), by.y=c('cluster','binterrogate'), all.y=TRUE, allow.cartesian=TRUE)
        seed.decays = to.drop[is.na(memb)][, memb := NULL]
        seed.decays[, binterrogate := V1]
        
        
        seed.decays[!is.na(sum.relative.risk), max.rr := max(sum.relative.risk), by='cluster']
        
        bin.to.add = seed.decays[sum.relative.risk == max.rr, bin.to.add := V1][!is.na(bin.to.add), c('cluster','bin.to.add')] %>% unique(by='cluster')
        
        colnames(bin.to.add)[2] = 'bin.discovered'
        pairs.discovered = merge(seed.decays, bin.to.add, by='cluster')
        
        pairs.discovered = pairs.discovered[V1 == bin.discovered]

        memb.dt[, numbins := .N, by='cluster']
        ##pairs.discovered[, numbins := .N, by='cluster']
        unique.cluster.numbins = memb.dt %>% unique(by='cluster')
        pairs.discovered = merge(unique.cluster.numbins[, c('cluster','numbins')], pairs.discovered, by='cluster')
        
        check = pairs.discovered[V1 == bin.discovered]

        print('check')
        print(check[cluster==4282])

        prop.check = check[, (choose(numbins, 2) - sum(is.na(bin.to.add))) / choose(numbins, 2), by='cluster']


        cluster.labeled = merge(seeds.loop, dist.decay.test, by='pair.hashes', allow.cartesian=TRUE)
        all.cluster.check = merge(cluster.labeled, bin.to.add, by.x=c('cluster','binterrogate'), by.y=c('cluster','bin.discovered'))

        paircount = seeds.loop[, .(all.count = .N), by='cluster']
        all.cluster.check = merge(all.cluster.check, paircount, by='cluster', all.x=TRUE)

        ##all.cluster.check[pval < pval.thresh & relative.risk>1, enrich.count := .N, by='cluster'] ##I don't think this will be toooooo consequential...
        all.cluster.check[pval < pval.thresh, enrich.count := .N, by='cluster']
        all.cluster.check[is.na(enrich.count), enrich.count := 0]
        all.cluster.check = all.cluster.check[, .(prop = enrich.count / all.count), by='cluster']
        all.cluster.check[, prop := max(prop), by='cluster']
        all.cluster.check[, kill := NA]
        all.cluster.check[prop < expansion.cutoff, kill := TRUE]

        ##prop.check[V1 < 0.5, kill := TRUE]
        print(prop.check)
        
        prop.check = all.cluster.check
        print(prop.check)
        
        pairs.discovered = pairs.discovered[V1 == bin.discovered]


        add.bins.1 = merge(pairs.discovered[, c('i','j','binterrogate','cluster')], pairwise.sym, by.x=c('i','binterrogate'), by.y=c('i','j'), all.x=TRUE)
        add.bins.2 = merge(pairs.discovered[, c('i','j','binterrogate','cluster')], pairwise.sym, by.x=c('j','binterrogate'), by.y=c('i','j'), all.x=TRUE)
        add.bins = rbind(add.bins.1, add.bins.2)


        seeds.loop.chunk = add.bins[, c('pair.hashes','cluster')] %>% unique(by=c('pair.hashes','cluster'))
        print(seeds.loop.chunk)
        
        seeds.loop = rbind(seeds.loop, seeds.loop.chunk) %>% unique(by=c('pair.hashes','cluster'))

        prop.check = unique(prop.check, by='cluster')
        seeds.loop = merge(seeds.loop, prop.check[, c('cluster','kill')], by='cluster')
        seeds.loop = seeds.loop[is.na(kill)] ##stop if you have less than 50% cheunks

        if(dim(seeds.loop)[[1]] == 0){ ##stopping condition
            break
        }
        
        seeds.loop[, kill := NULL]
        
        memb.dt.chunk = merge(seeds.loop, pairwise.sym[, c('pair.hashes','i')], by='pair.hashes')
        ##memb.dt.chunk = memb.dt.chunk[, .(unlist(agg)), by='cluster']
        memb.dt.chunk = unique(memb.dt.chunk, by=c('i','cluster'))
        colnames(memb.dt.chunk)[colnames(memb.dt.chunk) == 'i'] = 'V1'
        memb.dt.chunk$memb = TRUE

        memb.dt.chunk[, numbins.chunk := .N, by='cluster']


        ###kill any bin where the discovered bin is not added to memb dt
        kill.list = merge(memb.dt, unique(memb.dt.chunk[, c('cluster','numbins.chunk')], by='cluster'), by='cluster')
        kill.list.cluster = kill.list[numbins.chunk == numbins]$cluster %>% unique
        seeds.loop = seeds.loop[!(cluster %in% kill.list.cluster)]
        
        memb.dt[, numbins := NULL]
        memb.dt.chunk[, numbins.chunk := NULL]
        memb.dt.chunk[, pair.hashes := NULL]
        memb.dt = rbind(memb.dt, memb.dt.chunk)
        memb.dt = unique(memb.dt, by=c('cluster','V1'))
        print(memb.dt[cluster==4282])
        print(memb.dt[cluster==30])
    }

    binsets.gr = bins[memb.dt$V1]
    binsets.gr$bid = memb.dt$cluster

    print(binsets.gr)
    print(memb.dt)
    binsets.gr$binid = memb.dt$V1
    return(binsets.gr)
   ##merge(dist.decay.test, asdf[, c('pair.hashes','cluster')], by='pair.hashes')
}    



symmetry_check = function(scored.annotate, synergy_results_scored, binsets) {
    synergy_results_scored[, total_annot := .N, by=c('annotation')]
    enriched.3way = scored.annotate[cardinality==3 & pval<0.05 & count>0]
    all.3way = scored.annotate[cardinality==3]

    every.bid = data.table(bid = scored.annotate$bid %>% unique)

    max.dists = merge(every.bid, enriched.3way[, c('bid','max.dist')], by='bid')
    min.dists = merge(every.bid, enriched.3way[, c('bid','min.dist')], by='bid')


    max.dists.all = merge(every.bid, all.3way[, c('bid','max.dist')], by='bid')
    min.dists.all = merge(every.bid, all.3way[, c('bid','min.dist')], by='bid')

    colnames(every.bid)[1] = 'bid.test'
    synergy_results_scored = synergy_results_scored[!is.na(fdr)]
    every.bid = merge(every.bid, synergy_results_scored[, c('bid','fdr')], by.x='bid.test', by.y='bid', all.x=TRUE)
    every.bid = every.bid[fdr<0.1]

    binsets$bid = factor(binsets$bid)
    every.bid = merge(every.bid, binsets[, c('bid','overall.cardinality')], by.x='bid.test', by.y='bid', allow.cartesian=TRUE)
    ##browser()
    every.bid = every.bid[overall.cardinality > 4 & fdr<0.1] %>% unique(by='bid.test')

    
    count.enriched = enriched.3way[, .(enrich.count = .N), by='bid']
    count.all = all.3way[, .(all.count = .N), by='bid']
    coverage = merge(count.enriched, count.all, by='bid')
    coverage[, coverage := enrich.count / all.count]



    ###we cannot
    
    max.dist.comp = every.bid[, muffle(ks.test(max.dists[bid==bid.test]$max.dist, max.dists.all[bid==bid.test]$max.dist, alternative='less'), exact=FALSE), by='bid.test']
    min.dist.comp = every.bid[, muffle(ks.test(min.dists[bid==bid.test]$min.dist, min.dists.all[bid==bid.test]$min.dist), exact=FALSE), by='bid.test']

    symmetry_syn = merge(synergy_results_scored[annotation=='chromunity'], min.dist.comp[, c('bid.test','p.value')], by.x='bid', by.y='bid.test', all.x=TRUE)
    symmetry_syn = merge(symmetry_syn, coverage[, c('bid','coverage', 'enrich.count')], by='bid', all.x=TRUE)

    symmetry_syn[is.na(enrich.count), enrich.count := 0]

    ##browser()
    colnames(symmetry_syn)[colnames(symmetry_syn) == 'p.value'] = 'symmetry.p.value'

    symmetry_syn[is.na(symmetry.p.value), symmetry.p.value := 1]
    asdf = every.bid[, mean(max.dists[bid==bid.test]$max.dist) > mean(max.dists.all[bid==bid.test]$max.dist), by='bid.test']

    ##symmetry_syn = merge(symmetry_syn, unique(binsets[, c('bid','cardinality')], by='bid'), by='bid')

    binsets[, num.chr := length(unique(seqnames)), by='bid']
    symmetry_syn = merge(symmetry_syn, synergy_results_scored[annotation=='shuffled', c('bid','fdr')], by='bid')
    symmetry_syn = merge(symmetry_syn, unique(binsets[, c('bid','num.chr','cardinality')], by='bid'), by='bid')
##    symmetry_syn_chunk1$chunk.number = chunk.number
    
    symmetry_syn[cardinality > 3, is.dispersed := (enrich.count > 1)]
    symmetry_syn[is.na(is.dispersed), is.dispersed := TRUE]


    ##browser()
    symmetry_syn[, strict.synergy := FALSE]
    symmetry_syn[fdr.x<0.1 & symmetry.p.value > 0.1  & coverage > 0.1 & is.dispersed==TRUE & fdr.y > 0.1, strict.synergy := TRUE] ###the only goood good
    ##symmetry_syn[fdr.x<0.1 & symmetry.p.value > 0.1, strict.synergy := TRUE] ###the only goood good

    symmetry_syn[strict.synergy==FALSE, fdr.x := 1]
    sym = symmetry_syn
    sym$annotation='strict.chromunity'


    colnames(sym)[8] = 'fdr'
    synergy_results_scored[, total_annot := .N, by=c('annotation')]
    symmetry_all = rbind(synergy_results_scored, sym[, c('bid','method','p','estimate','ci.lower','ci.upper','effect','fdr','annotation', 'total_annot')])

    
    return(symmetry_all)

}

   
plot_concatemer_histogram = function(concats.gr, filename='plot.pdf') {
    concatemers = gr2dt(concats.gr)
    concatemers[, bincount := .N, by='read_idx']
    unique.concatemers = unique(concatemers, by='read_idx')
    ggsave(paste0('~/public_html/', filename), gghistogram(unique.concatemers[bincount<=10], x='bincount', fill='lightgrey', bins=10) +  theme(axis.title = element_text(hjust = 0.5, size = 20)))
}



## evaluate_synergy_stinky = function(res, leave_out_concatemers, chid.to.test, chromosome = NULL, filter_binsets = TRUE, folder = NULL, rebin_thresh = 0.85, mc.cores = 20, numchunks = mc.cores*200 + 1, dont_rebin=FALSE, sub.binset.order=5, resolution=1e4) {


##     if (!dir.exists(folder)) {
##         stop("output folder does not exist")
##     }
##     if(is.null(chromosome)){
##         chromosome = c(paste0("chr", c(as.character(1:22), "X")))
##     }
##     this.bad = load_bad_regions(chromosome)
    
##     tiles = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), resolution)
##     gc5b = readRDS("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/gc.38.rds")
##     ## create a list of covariates
##     cov_list = list(gc5b)
##     ## Specify what kind of covariate it is. Score will be aggregated over bins while the number of intervals will be calculated otherwise.
##     names(cov_list) <- c("score:gc.cov")
##     ## Make the covariate object
##     gc_cov = covariate(name = c("gc"), type = c("numeric"), field = c("score"), data = cov_list)
##     gc.cov = gc_cov

##     ##handle all that rebinning shit UPSTREAM please thank you 
    
##     ## this.dat.chrom = data.table()
##     ## if(dont_rebin) this.chrom = gr2dt(concatemers)
##     ## if (filter_binsets) {
##     ##     this.chrom = this.chrom[support > summary(unique(this.chrom[, .(support, chid)])$support)[3]]####filtering out all bins below median support
##     ## }
##     ## this.chrom.w = unique(chid.to.test)
##     ## if(dont_rebin) {
##     ##     print('yay')
##     ##     this.chrom.dt = this.chrom
##     ## } else {
##     ##     this.chrom.dt = rebin_community(res$concatemers, this.chrom.w, resolution=resolution)
##     ## }

##     ####

##     ##browser()
##     this.chrom.dt = gr2dt(res$binsets)
    
##     this.chrom.dt[, cardinality := .N, by = chid]
##     this.chrom.dt = na.omit(this.chrom.dt)
##     ##
##     this.all.dat = copy(this.chrom.dt)
##     this.all.dat = this.all.dat[cardinality > 2]
##     this.all.dat[, bid := chid]
##     this.all.dat = this.all.dat[seqnames %in% chromosome]
##     this.all.dat[, overall.cardinality := cardinality, by = bid]
##     this.chrom.card = unique(this.all.dat[, .(overall.cardinality, bid)])
##     ######HELLL NAH
##     this.all.dat = this.all.dat[cardinality < 100]
##     ####
##     this.sub.parq = leave_out_concatemers
##     this.sub.parq$cid = this.sub.parq$read_idx
##     this.all.dat[, binid := .I]
##     ##
##     ###filtering out the bad regions
##     ###this makes sense why it wouldn't be here for RE chromunity cause you're only looking at annotated regions anyway
##     this.all.dat$bid = this.all.dat$chid
##     this.bad.chrom = unique(as.character((dt2gr(this.all.dat) %&% (this.bad))$bid))                                                 
##     this.all.dat = this.all.dat[!bid %in% this.bad.chrom]  
##     #browser()
##     #debug(annotate)
##     chrom.annot.output = annotate(binsets = dt2gr(this.all.dat[, bid := chid]),
##                               k = 5,
##                               concatemers = this.sub.parq,
##                               covariates = gc.cov, resolution = resolution,
##                               mc.cores = mc.cores, numchunks = numchunks)

##     this.chrom.dat = chrom.annot.output[[1]]
##     this.chrom.dat.2 = this.chrom.dat
##     this.chrom.dat = merge(this.chrom.dat, this.chrom.card[, bid := as.factor(bid)], by = "bid")
##     this.chrom.dat[, annotation := "chromunity"]
##     set.seed(198)
##     message("generating background binsets for the model")
##     ##
##     #browser()
##     if(length(chromosome) > 1){
##         back.dt = re_background(binsets = dt2gr(this.all.dat), n = 1000, resolution = resolution)#, mc.cores=mc.cores)
##     } else {
##         back.dt = sliding_window_background(chromosome = chromosome, binsets = dt2gr(this.all.dat), n = 1000, resolution = resolution)#, mc.cores=mc.cores)
##     }
    
##     upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T)
##     setkeyv(back.dt, c("seqnames", "start"))
##     back.dt[, V1 := NULL]
##     back.dt = na.omit(back.dt)
##     back.dt = back.dt[!bid %in% back.dt[width < (resolution-1)]$bid]
##     back.dt = gr2dt(gr.reduce(dt2gr(back.dt), by = "bid"))
##     back.dt$bid <- as.factor(back.dt$bid)
##     back.dt = merge(back.dt, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
##     back.dt = back.dt[end < V2][start < V2]
##     back.dt[, overall.cardinality := .N, by = bid]
##     back.dt = back.dt[overall.cardinality > 1]
##     ##
##     this.card = unique(back.dt[, .(bid, overall.cardinality)])
##     back_gr = dt2gr(back.dt)
##     message("extracting background binsets distances")
##     back.train.output = annotate(binsets = dt2gr(back.dt),
##                                    interchromosomal.table = NULL, #all.hr.dt.mean,
##                                    gg = NULL,
##                                    concatemers = this.sub.parq, k = 5,
##                                    covariates = gc.cov, resolution = resolution,
##                                    mc.cores = mc.cores, numchunks = numchunks)
##     this.back.train.dat = back.train.output[[1]]
##     this.back.train.dat = merge(this.back.train.dat, this.card[, bid := as.factor(bid)], by = "bid")
##     this.back.train.dat[, annotation := "random"]
##     this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][width <= resolution]$bid)]
##     this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][min.dist < resolution]$bid)]
##     back.dt[, binid := .I]
##     this.bad.train = unique(as.character((back_gr %&% (this.bad))$bid))
##     this.back.train.dat = this.back.train.dat[!bid %in% this.bad.train]
##     this.back.train.dat = this.back.train.dat[!bid %in% this.back.train.dat[, .(sum(count)), by = bid][V1 == 0]$bid]
## ####    
##     back.model = fit(na.omit(this.back.train.dat)[sum.counts > 0][, setdiff(names(this.back.train.dat), c('overall.cardinality', 'chr', 'annotation')), with = F])
##     this.chrom.dat = sscore(this.chrom.dat, model = back.model)
##     message("generating random binsets for testing")
##     n.chrom = length(unique(this.chrom.dat$bid))
##     this.all.dat = this.all.dat[seqnames %in% chromosome]

##     if(length(chromosome) > 1){
##         back.test = re_background(binsets = dt2gr(this.all.dat),
##                                           n = n.chrom,
##                                           resolution = resolution)
##      } else {
##         back.test = sliding_window_background(binsets = dt2gr(this.all.dat), chromosome = chromosome,
##                                             n = n.chrom,
##                                             resolution = resolution)
##     }                      

##     back.test[, V1 := NULL]
##     back.test = na.omit(back.test)
##     setkeyv(back.test, c("seqnames", "start"))
##     back.test[, start := ifelse(start < 0, 0, start)]
##     back.test = back.test[!bid %in% back.test[width < (resolution-1)]$bid]
##     back.test = gr2dt(gr.reduce(dt2gr(back.test), by = "bid"))
##     back.test = merge(back.test, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
##     back.test = back.test[end < V2][start < V2]
##     back.test[, overall.cardinality := .N, by = bid]
##     back.test = back.test[overall.cardinality > 1]
##     back_test = dt2gr(back.test)
## ##
##     back.test$bid <- as.factor(back.test$bid)
##     this.back.card = unique(back.test[, .(bid, overall.cardinality)])
##     bid.back = unique(back.test$bid)
##     message("extracting random binsets distances") 
##     back.test.output = annotate(binsets = dt2gr(back.test), k=sub.binset.order,
##                                   concatemers = this.sub.parq,
##                                   covariates = gc.cov, resolution = resolution,
##                                   mc.cores = mc.cores, numchunks=numchunks)
##     this.back.test.dat = back.test.output[[1]]
##     this.back.test.dat = merge(this.back.test.dat, this.back.card, by = "bid")
##     this.back.test.dat[, annotation := "random"]
##     this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][width <= resolution]$bid)]
## ###
##     this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][min.dist < resolution]$bid)]
## ####
## ####
##     this.bad.dat = unique(as.character((back_test %&% (this.bad))$bid))
##     this.back.test.dat = this.back.test.dat[!bid %in% this.bad.dat]
##     back_test = gr2dt(back_test)[bid %in% unique(this.back.test.dat$bid)]
##     #this.back.test.dat = sscore(this.back.test.dat, model = back.model)
## ####
##     back.test[, binid := .I]
## #####
##     set.seed(178)
##     all.bid = unique(this.all.dat$bid)
## ###
##     ##browser()
##     message("starting random walks")
##     message("extracting shuffled binsets distances")

##     ##browser()
##     ##nr = 1
    
##     ##resolution=1e4

##     ##all.bid
##     ##
##     browser()
##     this.all.dat.shuff5 = pbmclapply(1:length(all.bid), function(nr){
##          this.clust = dt2gr(this.all.dat[bid == all.bid[nr]])
##          this.chrs = .chr2str(as.character(unique(seqnames(this.clust))))
##          this.clust.wind = gr.reduce(this.clust+2e6)
##          upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T)
##          this.clust.wind = (gr2dt(this.clust.wind)[, start := ifelse(start < 0, 1, start)])
##          this.clust.wind = merge(this.clust.wind, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
##          this.clust.wind[, end := ifelse(end > V2, V2, end)]
##          this.clust.wind = dt2gr(this.clust.wind)

##          ##this.clust.wind
##          this.sub.win = gr2dt(this.sub.parq %&% this.clust.wind)
##          ##this.sub.win[, new.count := .N, by = read_idx]

##          ###throw this in there
         

         
##          this.tiles.orig = gr.tile(this.clust.wind, resolution)
         
## ##         this.tiles.orig$binid = 1:length(this.tiles.orig)
## ##         this.sub.win = bin_concatemers(dt2gr(this.sub.win), this.tiles.orig)
## ##         this.sub.win = unique(this.sub.win, by=c('read_idx','binid'))
##          this.sub.win[, new.count := .N, by = read_idx]
##          ##
##          card = unique(this.sub.win[, .(read_idx, new.count)])
##          this.steps = sum(card$new.count)
         
##          if (nrow(this.sub.win) > 0){
##              this.tgm = cocount(dt2gr(this.sub.win), bins = (this.tiles.orig), by = 'read_idx', full = T)
##              ##debug(shuffle_concatemers)
##              ##shuff.concats = shuffle_concatemers(this.sub.win, this.tgm)
##              ##shuff.concats$cid = shuff.concats$read_idx

##              ##ppdf(plot(c(shuff.contacts$gtrack(name='shuff', clim=c(0,100)), this.tgm$gtrack(name='normal',clim=c(0,100))), gr.reduce(this.tiles.orig) + 3e5))
##              ##ppdf(plot(this.tgm$gtrack(clim=c(0,100)), gr.reduce(this.tiles.orig) + 3e5))
##              ##shuff.contacts = cocount(dt2gr(shuff.concats), bins = this.tiles.orig, by='read_idx', full=T)
##              A = this.tgm$mat %>% as.matrix
##              rownames(A) <- NULL
##              colnames(A) <- NULL
##              A[cbind(1:nrow(A), 1:nrow(A))] = 0
##              A = A + t(A)
##              An = A 
##              An = round(1+10*An/min(An[An>0]))
##              edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
##              ## ##
##              G = graph.edgelist(edges[, cbind(row, col)])
##              RW = random_walk(G, start = 1, steps = sum(card$new.count)) %>% as.numeric
##              ## ##
##              rm(G)
##              rm(edges)
##              gc()
##              out = this.tgm$gr[RW]%>% gr2dt()
##              out$read_idx = card[, rep(read_idx, new.count)] 
##              out[, bid := all.bid[nr]]
##              out[, cid := read_idx]

##              ##debug(annotate)
##              sh.output = annotate(binsets = this.clust, verbose = F,
##                                           k = 5,
##                                           concatemers = dt2gr(out),
##                                           covariates = gc.cov, resolution = resolution, 
##                                           mc.cores = 5, numchunks = numchunks)
##              ## shuff.concats$cid = shuff.concats$read_idx
##              ## sh.output.2 = annotate(binsets = this.clust, verbose = F,
##              ##                              k = 5,
##              ##                              concatemers = dt2gr(shuff.concats),
##              ##                              covariates = gc.cov, resolution = resolution, 
##              ##                              mc.cores = 10, numchunks = numchunks)
             
##              this.chrom.sh.dat = sh.output[[1]]
##          } else {
##              this.chrom.sh.dat = data.table(NA)
##          }
##          return(this.chrom.sh.dat)
##     }, mc.cores = 5, mc.preschedule = T)
##     #this.shuff.chrom = tryCatch(rbindlist(the.shuff.list, fill = T), error = function(e) NULL)


##     ##browser()

##     this.all.dat.shuff5.ne = this.all.dat.shuff5[sapply(this.all.dat.shuff5, function(x) !inherits(x, "try-error"))]
##     this.all.sh = rbindlist(this.all.dat.shuff5.ne, fill =T)
## ###
##     this.all.sh = sscore(this.all.sh, model = back.model)
##     sh.chrom = tryCatch(synergy(binsets = dt2gr(this.all.dat), annotated.binsets = this.all.sh, model = back.model), error = function(e) NULL)
##     sh.chrom$fdr = signif(p.adjust(sh.chrom$p, "BH"), 2)
##     ##
## #####
##     this.back.test.dat = sscore(this.back.test.dat, model = back.model)
##     theta = back.model$model$theta
##     s.chrom = synergy(binsets = dt2gr(this.all.dat), #theta = back.model$model$theta,
##                       annotated.binsets = this.chrom.dat, model = back.model)
##     s.chrom$fdr = signif(p.adjust(s.chrom$p, "BH"), 2)
##     b.chrom = synergy(binsets = dt2gr(back.test), #theta = back.model$model$theta,
##                       annotated.binsets = na.omit(this.back.test.dat), model = back.model)
##     b.chrom$fdr = signif(p.adjust(b.chrom$p, "BH"), 2)
##     ##
##     synergy.inter.EP = rbind(s.chrom[, annotation := "chromunity"],
##                              b.chrom[, annotation := "random"],

                            
##                              ##make_barplot(synergy.inter.EP)
##     if (!is.null(folder)) {
##         saveRDS(this.chrom.dat, paste0(folder,'chrom_annotate.rds'))
##         saveRDS(this.back.test.dat, paste0(folder,'back_annotate.rds'))
##         saveRDS(this.all.sh, paste0(folder,'shuffled_annotate.rds'))
##         saveRDS(this.all.dat, paste0(folder,'binsets.rds'))
##         saveRDS(back.model, paste0(folder,'back_model.rds'))
##         saveRDS(synergy.inter.EP, paste0(folder,'synergy_results.rds'))
##     }
##     return(synergy.inter.EP)
## }}




plot_model_output = function(dt.small.model, pair.hash, bins, plotname='plot.pdf', max.y=NULL, cluster=NULL, view_range=NULL, fdr.thresh=0.25, pairwise=NULL) {
    grange.out = grange_model_prediction(dt.small.model, pair.hash, bins, pairwise, cluster=cluster, fdr.thresh=fdr.thresh)

    model.gr = grange.out[[1]]
    ##view_range = grange.out[[2]]
   
    viewpoint = grange.out[[3]]

    if(is.null(view_range)) {
        ##view_range = viewpoint[[1]] + 1e6
        view_range = (model.gr %Q% (fdr<fdr.thresh) + 1e6) %>% gr.reduce
        dt.range = rbind(gr2dt(view_range), gr2dt(viewpoint[[1]]), gr2dt(viewpoint[[2]]), fill=TRUE)
        view_range = dt2gr(dt.range)
    }
    if(is.null(max.y)){
        max.y = model.gr$num.concats %>% max
    }

    model.gr$log.pval = -log10(model.gr$pval)

    if(is.null(cluster)){
        ppdf(plot(c(gTrack(model.gr, y.field='relative.risk', bars=T, name='log2(O/E)', y0=0), gTrack(model.gr, y.field='log.pval', bars=T, name='-log10(pval)'), gTrack(model.gr, y.field='num.concats.pred', bars=T, name='prediction', y1=max.y), gTrack(model.gr, y.field='num.concats', bars=T, name='actual', y1=max.y), gTrack(viewpoint, name='viewpoint'), gTrack(features$h27ac, y1 = 100, y.field = "V9", bars = T, height = 2, name='H3K27ac peaks')), view_range+1e6), plotname)
    } else {
        ppdf(plot(c(gTrack(model.gr, y.field='relative.risk', bars=T, name='log2(O/E)', y0=0), gTrack(viewpoint, name='viewpoint')), view_range+5e6), plotname)
    }
    return(list(view_range, max.y))
}


grange_model_prediction = function(dt.small.model, pair.hash, bins, pairwise=NULL, cluster=NULL, fdr.thresh=0.25){

    if(!is.null(cluster)){
        dt.small.model$pair.hashes.archive = dt.small.model$pair.hashes ##goofy
        dt.small.model$pair.hashes = dt.small.model$cluster
        dt.small.model$relative.risk = dt.small.model$sum.relative.risk
        dt.small.model = unique(dt.small.model, by=c('cluster','binterrogate'))
        pair.hash = cluster
    }
    viewpoint.bins = c(dt.small.model[pair.hashes==pair.hash]$i, dt.small.model[pair.hashes==pair.hash]$j) %>% unique
    dt.small.model = dt.small.model[!(binterrogate %in% viewpoint.bins)]
    model.gr = bins[dt.small.model[pair.hashes==pair.hash]$binterrogate]
    model.gr$num.concats = dt.small.model[pair.hashes==pair.hash]$num.concats
    model.gr$num.concats.pred = dt.small.model[pair.hashes==pair.hash]$num.concats.pred
    model.gr$sum.pairwise.contacts = dt.small.model[pair.hashes==pair.hash]$sum.pairwise.contacts
    model.gr$enrichment = dt.small.model[pair.hashes==pair.hash]$enrichment
    model.gr$obs.over.expected = log2(model.gr$num.concats/model.gr$num.concats.pred)
    model.gr$fdr = dt.small.model[pair.hashes==pair.hash]$fdr
    
    dt.small.model.copy = dt.small.model[pair.hashes == pair.hash]# & pval > 0.05, relative.risk := 0]
    dt.small.model.copy[fdr > fdr.thresh, relative.risk := 0]

    ##version with pval thresh
    ##dt.small.model.copy[pval > 0.05, relative.risk := 0]
    
    model.gr$relative.risk = dt.small.model.copy$relative.risk
    ##model.gr$num.concats = dt.small.model[pair.hashes==pair.hash]$num.concats
    model.gr$pval = dt.small.model[pair.hashes==pair.hash]$pval
    
    
    if(!is.null(cluster)){
        all.bins = merge(pairwise[, c('pair.hashes','agg')], dt.small.model[pair.hashes==pair.hash, c('pair.hashes.archive')], by.x='pair.hashes', by.y='pair.hashes.archive')
        all.bins.unlist = all.bins[, .(bins = unlist(agg))]
        all.bins = all.bins.unlist$bins %>% unique
        viewpoint = split(bins[all.bins], bins[all.bins]$binid) %>% unname
    } else {
        bin1 = bins[unique(dt.small.model[pair.hashes==pair.hash]$i)]
        bin2 = bins[unique(dt.small.model[pair.hashes==pair.hash]$j)]
        viewpoint = GRangesList(bin1, bin2)
    }
    view_range = (bins[dt.small.model[pair.hashes==pair.hash]$binterrogate] + 1e6)  %>% gr.reduce
    return(list(model.gr, view_range, viewpoint))
}


aggregate_synergy_results = function(dir, toplevel=TRUE, strict.check=FALSE) {

    dirs = list.dirs(dir, recursive=FALSE)
    print(dirs)
#####aggregate results function here basically

    ##toplevel=TRUE
    agg.synergy = pbmclapply(dirs, mc.cores = 2, mc.preschedule=FALSE, function(dir) {        
        ##if(dir=="GM12878_dist_decay_all_chr/dist_decay_chr21_knn25_kmin5_resolution50000"){ hardcoded trash
        ##    return(NULL)
        ##}
        if(toplevel==TRUE){
            synergy.chunk = readRDS(paste0(dir, '/synergy_outputs/synergy_results.rds'))
            binsets = readRDS(paste0(dir, '/synergy_outputs/binsets.rds'))
        } else {
            synergy.chunk = readRDS(paste0(dir, '/synergy_results.rds'))
            binsets = readRDS(paste0(dir, '/binsets.rds'))
        }
        if(strict.check==TRUE & toplevel==TRUE) {
            annotate = readRDS(paste0(dir, '/synergy_outputs/chrom_annotate.rds'))
            synergy.chunk = symmetry_check(annotate, synergy.chunk, binsets)
        } else if (strict.check == TRUE & toplevel==FALSE) {
            annotate = readRDS(paste0(dir, '/chrom_annotate.rds'))
            synergy.chunk = symmetry_check(annotate, synergy.chunk, binsets)
        }
        return(synergy.chunk)
    })
    agg.synergy = rbindlist(agg.synergy)
    agg.synergy[, .N, by='annotation']


    agg.synsets = pbmclapply(dirs, mc.cores = 2, function(dir) {
        if(toplevel==TRUE){
            synergy.chunk = readRDS(paste0(dir, '/synergy_outputs/synergy_results.rds'))
            binsets = readRDS(paste0(dir, '/synergy_outputs/binsets.rds'))
        } else {
            synergy.chunk = readRDS(paste0(dir, '/synergy_results.rds'))
            binsets = readRDS(paste0(dir, '/binsets.rds'))
        }
        ##synergy.chunk = readRDS(paste0(dir, '/synergy_results.rds'))
        ##binsets = readRDS(paste0(dir, '/binsets.rds'))
        binsets$bid = factor(binsets$bid)
        synsets = merge(binsets, synergy.chunk[annotation=='chromunity', c('bid','fdr')], by='bid')
        return(synsets)
    })


    
    agg.synsets = rbindlist(agg.synsets)
    synergies = dt2gr(agg.synsets[fdr<.1])
    return(list(agg.synergy, synergies))
}





evaluate_synergy_experimental = function(res, leave_out_concatemers, chid.to.test, chromosome = NULL, filter_binsets = TRUE, folder = NULL, rebin_thresh = 0.85, mc.cores = 20, numchunks = mc.cores*200 + 1, dont_rebin=FALSE, sub.binset.order=3, resolution=1e4, remove.bad=TRUE) {

    if (!dir.exists(folder)) {
        stop("output folder does not exist")
    }
    if(is.null(chromosome)){
        chromosome = c(paste0("chr", c(as.character(1:22), "X", "Y"))) ###we should add the Y CHROMOSOME!!!
        all.bad = pbmclapply(chromosome, mc.preschedule=FALSE, function(chr) {
            bad.gr = muffle(load_bad_regions(chr))
            bad.dt = gr2dt(bad.gr)
            return(bad.dt)
        })
        this.bad = rbindlist(all.bad) %>% dt2gr
    } else {
        this.bad = load_bad_regions(chromosome)
    }
    tiles = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), resolution)
    gc5b = readRDS("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/gc.38.rds")
    ## create a list of covariates
    cov_list = list(gc5b)
    ## Specify what kind of covariate it is. Score will be aggregated over bins while the number of intervals will be calculated otherwise.
    names(cov_list) <- c("score:gc.cov")
    ## Make the covariate object
    gc_cov = covariate(name = c("gc"), type = c("numeric"), field = c("score"), data = cov_list)
    gc.cov = gc_cov

    if(!(class(res) == 'data.table'))
        this.chrom.dt = gr2dt(res$binsets)
    else
        this.chrom.dt = res
    
    this.chrom.dt[, cardinality := .N, by = chid]
    this.chrom.dt = na.omit(this.chrom.dt)
    ##
    this.all.dat = copy(this.chrom.dt)
    this.all.dat = this.all.dat[cardinality > 2]
    this.all.dat[, bid := chid]
    this.all.dat = this.all.dat[seqnames %in% chromosome]
    this.all.dat[, overall.cardinality := cardinality, by = bid]
    this.chrom.card = unique(this.all.dat[, .(overall.cardinality, bid)])
    ######HELLL NAH
    this.all.dat = this.all.dat[cardinality < 100]
    ####
    this.sub.parq = leave_out_concatemers
    this.sub.parq$cid = this.sub.parq$read_idx
    this.all.dat[, binid := .I]
    ##
    ###filtering out the bad regions
    ###this makes sense why it wouldn't be here for RE chromunity cause you're only looking at annotated regions anyway

###remove only bins which intersect bad region, not all

    this.all.dat$bid = this.all.dat$chid
    this.bad.chrom = unique(as.character((dt2gr(this.all.dat) %&% (this.bad))$bid))                                                 
    print(this.bad.chrom)
    if (remove.bad)
        this.all.dat = this.all.dat[!bid %in% this.bad.chrom]  
                                        #browser()
    #debug(annotate)
    chrom.annot.output = annotate(binsets = dt2gr(this.all.dat[, bid := chid]),
                              k = sub.binset.order,
                              concatemers = this.sub.parq,
                              covariates = gc.cov, resolution = resolution,
                              mc.cores = mc.cores, numchunks = numchunks)

    this.chrom.dat = chrom.annot.output[[1]]
    this.chrom.dat.2 = this.chrom.dat
    this.chrom.dat = merge(this.chrom.dat, this.chrom.card[, bid := as.factor(bid)], by = "bid")
    this.chrom.dat[, annotation := "chromunity"]
    set.seed(198)
    message("generating background binsets for the model")
    ##
                                        #browser()

    n.chrom = length(unique(this.chrom.dat$bid))
    if(length(chromosome) > 1){
        back.dt = help_bins_2(binsets = dt2gr(this.all.dat), n = (n.chrom*3), resolution = resolution)#, mc.cores=mc.cores)
    } else {

        back.dt = sliding_window_background(chromosome = chromosome, binsets = dt2gr(this.all.dat), n = (n.chrom*3), resolution = resolution)#, mc.cores=mc.cores)
    }
    
    upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T)
    setkeyv(back.dt, c("seqnames", "start"))
    back.dt[, V1 := NULL]
    back.dt = na.omit(back.dt)
    back.dt = back.dt[!bid %in% back.dt[width < (resolution-1)]$bid]
    back.dt = gr2dt(gr.reduce(dt2gr(back.dt), by = "bid"))
    back.dt$bid <- as.factor(back.dt$bid)
    back.dt = merge(back.dt, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
    back.dt = back.dt[end < V2][start < V2]
    back.dt[, overall.cardinality := .N, by = bid]
    back.dt = back.dt[overall.cardinality > 1]
    ##
    this.card = unique(back.dt[, .(bid, overall.cardinality)])
    back_gr = dt2gr(back.dt)
    message("extracting background binsets distances")



    back.train.output = annotate(binsets = dt2gr(back.dt),
                                   interchromosomal.table = NULL, #all.hr.dt.mean,
                                   gg = NULL,
                                   concatemers = this.sub.parq, k = sub.binset.order,
                                   covariates = gc.cov, resolution = resolution,
                                   mc.cores = mc.cores, numchunks = numchunks)
    this.back.train.dat = back.train.output[[1]]
    this.back.train.dat = merge(this.back.train.dat, this.card[, bid := as.factor(bid)], by = "bid")
    this.back.train.dat[, annotation := "random"]
    this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][width <= resolution]$bid)]
    this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][min.dist < resolution]$bid)]
    back.dt[, binid := .I]
    this.bad.train = unique(as.character((back_gr %&% (this.bad))$bid))
    this.back.train.dat = this.back.train.dat[!bid %in% this.bad.train]
    this.back.train.dat = this.back.train.dat[!bid %in% this.back.train.dat[, .(sum(count)), by = bid][V1 == 0]$bid]
####    

    ##browser()
    
    back.model = fit(na.omit(this.back.train.dat)[sum.counts > 0][, setdiff(names(this.back.train.dat), c('overall.cardinality', 'chr', 'annotation')), with = F])


    
    this.chrom.dat = sscore(this.chrom.dat, model = back.model)
    message("generating random binsets for testing")
    n.chrom = length(unique(this.chrom.dat$bid))
    this.all.dat = this.all.dat[seqnames %in% chromosome]

    if(length(chromosome) > 1){
        back.test = help_bins_2(binsets = dt2gr(this.all.dat),
                                          n = n.chrom,
                                          resolution = resolution)
     } else {
        back.test = sliding_window_background(binsets = dt2gr(this.all.dat), chromosome = chromosome,
                                            n = n.chrom,
                                            resolution = resolution)
    }                      

    back.test[, V1 := NULL]
    back.test = na.omit(back.test)
    setkeyv(back.test, c("seqnames", "start"))
    back.test[, start := ifelse(start < 0, 0, start)]
    back.test = back.test[!bid %in% back.test[width < (resolution-1)]$bid]
    back.test = gr2dt(gr.reduce(dt2gr(back.test), by = "bid"))
    back.test = merge(back.test, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
    back.test = back.test[end < V2][start < V2]
    back.test[, overall.cardinality := .N, by = bid]
    back.test = back.test[overall.cardinality > 1]
    back_test = dt2gr(back.test)
###
    back.test$bid <- as.factor(back.test$bid)
    this.back.card = unique(back.test[, .(bid, overall.cardinality)])
    bid.back = unique(back.test$bid)
    message("extracting random binsets distances") 
    back.test.output = annotate(binsets = dt2gr(back.test), k=sub.binset.order,
                                  concatemers = this.sub.parq,
                                  covariates = gc.cov, resolution = resolution,
                                  mc.cores = mc.cores, numchunks=numchunks)
    this.back.test.dat = back.test.output[[1]]
    this.back.test.dat = merge(this.back.test.dat, this.back.card, by = "bid")
    this.back.test.dat[, annotation := "random"]
    this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][width <= resolution]$bid)]
###
    this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][min.dist < resolution]$bid)]
####
####
    this.bad.dat = unique(as.character((back_test %&% (this.bad))$bid))
    this.back.test.dat = this.back.test.dat[!bid %in% this.bad.dat]
    back_test = gr2dt(back_test)[bid %in% unique(this.back.test.dat$bid)]
    #this.back.test.dat = sscore(this.back.test.dat, model = back.model)
####
    back.test[, binid := .I]
#####
    set.seed(178)
    all.bid = unique(this.all.dat$bid)
###
    ##browser()
    message("starting random walks")
    message("extracting shuffled binsets distances")

    ##browser()
    ##nr = 1
    
    ##resolution=1e4

    ##all.bid
    ##
    ##browser()
    s.chrom = synergy(binsets = dt2gr(this.all.dat), #theta = back.model$model$theta,
                      annotated.binsets = this.chrom.dat, model = back.model)
    s.chrom$fdr = signif(p.adjust(s.chrom$p, "BH"), 2)
    print('the hit rate!!!')
    print(s.chrom[, .N])
    print(s.chrom[fdr<0.1, .N])

    saveRDS(s.chrom, paste0(folder,'synergy_results.rds'))
    saveRDS(this.chrom.dat, paste0(folder,'chrom_annotate.rds'))
    
    this.all.dat.shuff5 = pbmclapply(1:length(all.bid), function(nr){
         this.clust = dt2gr(this.all.dat[bid == all.bid[nr]])
         this.chrs = .chr2str(as.character(unique(seqnames(this.clust))))
         this.clust.wind = gr.reduce(this.clust+2e6)
         upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T)
         this.clust.wind = (gr2dt(this.clust.wind)[, start := ifelse(start < 0, 1, start)])
         this.clust.wind = merge(this.clust.wind, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
         this.clust.wind[, end := ifelse(end > V2, V2, end)]
         this.clust.wind = dt2gr(this.clust.wind)

         ##this.clust.wind
         this.sub.win = gr2dt(this.sub.parq %&% this.clust.wind)
         ##this.sub.win[, new.count := .N, by = read_idx]

         ###throw this in there
         

         
         this.tiles.orig = gr.tile(this.clust.wind, resolution)
         
         this.tiles.orig$binid = 1:length(this.tiles.orig)
         this.sub.win = bin_concatemers(dt2gr(this.sub.win), this.tiles.orig)
         this.sub.win = unique(this.sub.win, by=c('read_idx','binid'))
         this.sub.win[, new.count := .N, by = read_idx]
         ##
         card = unique(this.sub.win[, .(read_idx, new.count)])
         this.steps = sum(card$new.count)

         if (nrow(this.sub.win) > 0){
             this.tgm = cocount(dt2gr(this.sub.win), bins = (this.tiles.orig), by = 'read_idx', full = T)
             ##debug(shuffle_concatemers)
             ##shuff.concats = shuffle_concatemers(this.sub.win, this.tgm)
             ##shuff.concats$cid = shuff.concats$read_idx

             ##ppdf(plot(c(shuff.contacts$gtrack(name='shuff', clim=c(0,100)), this.tgm$gtrack(name='normal',clim=c(0,100))), gr.reduce(this.tiles.orig) + 3e5))
             ##ppdf(plot(this.tgm$gtrack(clim=c(0,100)), gr.reduce(this.tiles.orig) + 3e5))
             ##shuff.contacts = cocount(dt2gr(shuff.concats), bins = this.tiles.orig, by='read_idx', full=T)
             A = this.tgm$mat %>% as.matrix
             rownames(A) <- NULL
             colnames(A) <- NULL
             A[cbind(1:nrow(A), 1:nrow(A))] = 0
             A = A + t(A)
             An = A 
             An = round(1+10*An/min(An[An>0]))
             edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
             ## ##
             G = graph.edgelist(edges[, cbind(row, col)])
             RW = random_walk(G, start = 1, steps = sum(card$new.count)) %>% as.numeric
             ## ##
             rm(G)
             rm(edges)
             gc()
             out = this.tgm$gr[RW]%>% gr2dt()
             out = out[1:sum(card$new.count)] ##weird bug
             out$read_idx = card[, rep(read_idx, new.count)] 
             out[, bid := all.bid[nr]]
             out[, cid := read_idx]

             ##debug(annotate)
             sh.output = annotate(binsets = this.clust, verbose = F,
                                          k = sub.binset.order,
                                          concatemers = dt2gr(out),
                                          covariates = gc.cov, resolution = resolution, 
                                          mc.cores = 1, numchunks = numchunks)
             ## shuff.concats$cid = shuff.concats$read_idx
             ## sh.output.2 = annotate(binsets = this.clust, verbose = F,
             ##                              k = 5,
             ##                              concatemers = dt2gr(shuff.concats),
             ##                              covariates = gc.cov, resolution = resolution, 
             ##                              mc.cores = 10, numchunks = numchunks)
             
             this.chrom.sh.dat = sh.output[[1]]
         } else {
             this.chrom.sh.dat = data.table(NA)
         }
         return(this.chrom.sh.dat)
    }, mc.cores = 5, mc.preschedule = T)
    #this.shuff.chrom = tryCatch(rbindlist(the.shuff.list, fill = T), error = function(e) NULL)


    ##browser()

    this.all.dat.shuff5.ne = this.all.dat.shuff5[sapply(this.all.dat.shuff5, function(x) !inherits(x, "try-error"))]
    this.all.sh = rbindlist(this.all.dat.shuff5.ne, fill =T)
###
    this.all.sh = sscore(this.all.sh, model = back.model)
    sh.chrom = tryCatch(synergy(binsets = dt2gr(this.all.dat), annotated.binsets = this.all.sh, model = back.model), error = function(e) NULL)
    sh.chrom$fdr = signif(p.adjust(sh.chrom$p, "BH"), 2)
    ##
#####
    this.back.test.dat = sscore(this.back.test.dat, model = back.model)
    theta = back.model$model$theta
    s.chrom = synergy(binsets = dt2gr(this.all.dat), #theta = back.model$model$theta,
                      annotated.binsets = this.chrom.dat, model = back.model)
    s.chrom$fdr = signif(p.adjust(s.chrom$p, "BH"), 2)
    b.chrom = synergy(binsets = dt2gr(back.test), #theta = back.model$model$theta,
                      annotated.binsets = na.omit(this.back.test.dat), model = back.model)
    b.chrom$fdr = signif(p.adjust(b.chrom$p, "BH"), 2)
    ##

    synergy.inter.EP = rbind(s.chrom[, annotation := "chromunity"],
                             b.chrom[, annotation := "random"], 
                             sh.chrom[, annotation := "shuffled"], fill = T)
    if (!is.null(folder)) {
        saveRDS(this.chrom.dat, paste0(folder,'chrom_annotate.rds'))
        saveRDS(this.back.test.dat, paste0(folder,'back_annotate.rds'))
        saveRDS(this.all.sh, paste0(folder,'shuffled_annotate.rds'))
        saveRDS(this.all.dat, paste0(folder,'binsets.rds'))
        saveRDS(back.model, paste0(folder,'back_model.rds'))
        saveRDS(synergy.inter.EP, paste0(folder,'synergy_results.rds'))
    }
    return(synergy.inter.EP)
}


#' @name interchr_dist_decay_binsets
#' @description
#'
#' This function performs the bin-set nomination procedure based on modeling the distance decay of higher order interaction frequencies from pairs of bins
#' 
#' 
#' @param concatemers GRanges of monomers with fields seqnames, start, end, and $cid specifying concatemer id, which will be counted across each binset
#' @param resolution integer specifying the bin width to use for the distance decay model
#' @param bins GRanges of bins which specify the 
#' @param interchromosomal.dist numeric scalar of "effective" distance for inter chromosomal bins [1e8]
#' @param training.chr Chromosome to use for training distance decay model
#' @param pair.thresh Pairwise contact value used as threshold for considering pairs in analysis
#' @param numchunks Number of chunks to create in annotating higher order distance decay
#' @param mask.bad.regions Will load human telomeres/centromeric regions and remove them from annotated distance decay
#' 
#' @param verbose logical flag
#' @param mc.cores integer how many cores to parallelize
#' @param threads used to set number of data table threads to use with setDTthreads function, segfaults may occur if >1
#' @author Jameson Orvis
#' @export
#' @return data.table of sub-binsets i.e. k-power set of binsets annotated with $count field representing covariates, ready for fitting, **one row per binset

interchr_dist_decay_binsets = function(concatemers, resolution=50000, bins=NULL, interchromosomal.distance = 1e8, training.chr = 'chr8', pair.thresh=50, numchunks=200, mask.bad.regions = TRUE, fdr.thresh=0.1, pval.thresh=0.05, expansion.cutoff=0.5, num.members=10, folder=NULL, chromosome=NULL, model=NULL, pairwise.trimmed=NULL, num.to.sample=500000) {

    if(is.null(chromosome)){
        chromosome = c(paste0("chr", c(as.character(1:22), "X")))
    }

    if(is.null(bins)){
        bins = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), width=resolution) %Q% (seqnames %in% chromosome)
    }
    bins$binid = 1:length(bins)

    if(is.null(model)){
        model = train_dist_decay_model_nozero(concatemers %Q% (seqnames==training.chr), bins %Q% (seqnames==training.chr), numchunks=20)
    }
    
    contact_matrix_unique = cocount(concatemers, bins = bins, by = 'read_idx')
    all.pairwise = contact_matrix_unique$dat
    all.pairwise$id = 1:dim(all.pairwise)[[1]]
    colnames(all.pairwise)[3] = 'pair.value'
    
    concatemers$cid = concatemers$read_idx
    binned.concats = bin_concatemers(concatemers, bins, max.slice=1e6, mc.cores=5)

    dt.concats = unique(binned.concats[, c('cidi','binid')], by=c('cidi','binid'))
    dt.concats.sort = dt.concats[order(binid, cidi)]
    dt.concats.sort[, count := .N, by='cidi']
    

    ##choose subset of bin-pairs S by thresholding
    colnames(all.pairwise)[4] = 'pair.hashes'

    ##browser()
    ##allow option of passing in pairwise trimmed manually 
    if(is.null(pairwise.trimmed)) {
        pairwise.trimmed = all.pairwise[pair.value >= pair.thresh]
        pairwise.trimmed[, dist := j-i]
        pairwise.trimmed = pairwise.trimmed[dist > 1]
    }
     
    if(mask.bad.regions==TRUE) {
        chromosome = c(paste0("chr", c(as.character(1:22), "X")))
        all.bad = pbmclapply(chromosome, mc.preschedule=FALSE, function(chr) {
            bad.gr = muffle(load_bad_regions(chr))
            bad.dt = gr2dt(bad.gr)
            return(bad.dt)
        })
        this.bad = rbindlist(all.bad) %>% dt2gr
        bad.bins = this.bad %*% bins
        pairwise.trimmed = pairwise.trimmed[!(i %in% bad.bins$binid)]
        pairwise.trimmed = pairwise.trimmed[!(j %in% bad.bins$binid)]
    }


    ###Should choose numchunks automatically here
    
    unique.pairs = pairwise.trimmed$pair.hashes %>% unique
    pair.splitting = data.table(unique.pairs, group=ceiling(runif(length(unique.pairs))*numchunks))
    pairwise.trimmed$group = pair.splitting$group

    pairwise.trimmed$id = pairwise.trimmed$pair.hashes
    all.pairwise$id = all.pairwise$pair.hashes

    pairwise.chunks = split(pairwise.trimmed, by='group')

    ##oh this isn't even bad! memory-wise
    ##browser()
    pairwise.chunk = pairwise.chunks[[1]]

    scored.chunks = pbmclapply(pairwise.chunks, mc.cores = 5, function(pairwise.chunk) {
        annot.chunk = count_3way_contacts(pairwise.chunk, dt.concats.sort, all.pairwise, bins, resolution=resolution, interchromosomal.distance = interchromosomal.distance)
        if(mask.bad.regions == TRUE){
            annot.chunk = annot.chunk[!(binterrogate %in% bad.bins$binid)]
        }
        dt.small.scored = score_distance_decay_poisson(annot.chunk, model)
        return(dt.small.scored)
    })

    ##spike.train.set = readRDS('spiked_train_set.rds')
    ##model = glm.nb(formula = fm, data=spike.train.set[, c('num.concats','value.a.ratio','value.b.ratio')], control=glm.control(maxit=500))
    ##genome.wide.dist.decay.2 = score_distance_decay(genome.wide.dist.decay, model)
    ##genome.wide.dist.decay = genome.wide.dist.decay.2
    ##browser()
    ##fdr.thresh=.25
    
    genome.wide.dist.decay = rbindlist(scored.chunks)

    
    ##quant.value = quantile(genome.wide.dist.decay[num.concats>1]$num.concats, 0.99)
    ##genome.wide.dist.decay[, norm.concats := num.concats / quant.value]
    genome.wide.dist.decay[, relative.risk := log2(num.concats/num.concats.pred)]

    ##Original code (strange)
    ##genome.wide.dist.decay[, relative.risk := log2(abs(num.concats - num.concats.pred))]
    ##genome.wide.dist.decay$sidi = 1:dim(genome.wide.dist.decay)[[1]]
    
    genome.wide.dist.decay$fdr = signif(p.adjust(genome.wide.dist.decay$pval, "BH"), 2)
    trimmed.dist.decay = genome.wide.dist.decay[fdr<fdr.thresh]
    
    bin.pair.network = create_bin_pair_network_efficient(all.pairwise, trimmed.dist.decay, rr.thresh=0)

    ##browser()
    ###rr.thresh!!
    if(!is.null(folder))
        saveRDS(genome.wide.dist.decay, paste0(folder,'dist_decay_archive.rds'))

    chrom = derive_binsets_from_network(bin.pair.network, pairwise=all.pairwise, binned.concats=binned.concats, bins=bins, rr.thresh=0, dist.decay.test=genome.wide.dist.decay[pval < pval.thresh], all.pairs.tested=NULL, pairwise.trimmed=pairwise.trimmed, expansion.cutoff = expansion.cutoff, pval.thresh = pval.thresh, fdr.thresh = fdr.thresh, num.members=num.members)
    
    return(chrom)
}

##subsetted to training chr
train_dist_decay_model_nozero = function(concatemers, bins, pair.thresh=50, numchunks=200, num.to.sample=500000, pairs.to.sample = 10000, mode='poisson'){
    contact_matrix_unique = cocount(concatemers, bins = bins, by = 'read_idx')
    all.pairwise = contact_matrix_unique$dat
    all.pairwise$id = 1:dim(all.pairwise)[[1]]
    colnames(all.pairwise)[3] = 'pair.value'
    
    concatemers$cid = concatemers$read_idx
    binned.concats = bin_concatemers(concatemers, bins, max.slice=1e6, mc.cores=5)
    dt.concats = unique(binned.concats[, c('cidi','binid')], by=c('cidi','binid'))
    dt.concats.sort = dt.concats[order(binid, cidi)]
    dt.concats.sort[, count := .N, by='cidi']
    dt.concats.sort

##pairwise = contact_matrix_unique$dat
    
    colnames(all.pairwise)[4] = 'pair.hashes'
    pairwise.trimmed = all.pairwise[pair.value >= pair.thresh]
    pairwise.trimmed[, dist := j-i]
    pairwise.trimmed = pairwise.trimmed[dist > 1]
    
    
    ##pairwise.trimmed$id = 1:dim(pairwise.trimmed)[[1]]
    unique.pairs = pairwise.trimmed$pair.hashes %>% unique
        
    #subsetting to make this more efficient
    if(length(unique.pairs) > pairs.to.sample){
        unique.pairs = unique.pairs[sample(pairs.to.sample)]
        pairwise.trimmed = pairwise.trimmed[pair.hashes %in% unique.pairs]
    }

    pair.splitting = data.table(unique.pairs, group=ceiling(runif(length(unique.pairs))*numchunks))
    pairwise.trimmed$group = pair.splitting$group

    pairwise.trimmed$id = pairwise.trimmed$pair.hashes
    all.pairwise$id = all.pairwise$pair.hashes

    pairwise.chunks = split(pairwise.trimmed, by='group')

##oh this isn't even bad! memory-wise

    scored.chunks = pbmclapply(pairwise.chunks, mc.cores = 5, function(pairwise.chunk) {
        annot.chunk = count_3way_contacts(pairwise.chunk, dt.concats.sort, all.pairwise, bins)
        return(annot.chunk)
    })
    dist.decay.train = rbindlist(scored.chunks)

    ##browser()
    
    print('training model')
    covariates = c('value.a.ratio','value.b.ratio')
    fmstring = paste('num.concats ~', paste(paste0('log(', covariates, ')', collapse = ' + ')))
    ##fmstring = paste0(fmstring, " + ", "offset(log(total.concats))") ##this sometimes does 

    fm = formula(fmstring)


        ##browser()
     ##num.to.sample = 500000
     ##if(num.to.sample > dim(dist.decay.train)[[1]]) {
     ##    num.to.sample = dim(dist.decay.train)[[1]]
     ##}


     ##total.concats = dist.decay.train[(dist.a < 50 | dist.b < 50), .(total.concats = sum(num.concats)), by=c('pair.hashes')]
     ##dist.decay.train = merge(dist.decay.train, total.concats, by='pair.hashes')

    if(num.to.sample > dim(dist.decay.train)[[1]]/2) {
        train.subset = dist.decay.train[dist.a <= 50 & dist.b <= 50]
    } else { 
        close.subset = dist.decay.train[dist.a <= 50 & dist.b <= 50][sample(.N, num.to.sample/2)]
        far.subset = dist.decay.train[dist.a > 50 & dist.b > 50][sample(.N, num.to.sample/8)]
        train.subset = rbind(close.subset, far.subset)
    }

    ##browser()
    if(mode=='poisson'){
        model = glm(formula = fm, data=train.subset[, c('num.concats','value.a.ratio','value.b.ratio')], control=glm.control(maxit=500), family='poisson')
    } else {
        model = glm.nb(formula = fm, data=train.subset[, c('num.concats','value.a.ratio','value.b.ratio')], control=glm.control(maxit=500))
    }

    return(model)
    
}


train_dist_decay_model = function(iterative.registry, iterative.registry.unlist.filterable, contact_matrix_unique, pairwise, concatemers, bins, pair.thresh=50) {

     dist.decay.train = annotate_distance_decay(concatemers,
                                                   bins,
                                                   min.value = pair.thresh,
                                                   window.size=50,
                                                   simplicial.complex=iterative.registry,
                                                   iterative.registry.unlist.filterable = iterative.registry.unlist.filterable,
                                                   pairwise=pairwise,
                                                   return.training=TRUE)



    ####total contacts within window
    
     total.concats = dist.decay.train[(dist.a < 50 | dist.b < 50), .(total.concats = sum(num.concats)), by=c('pair.hashes')]

     dist.decay.train = merge(dist.decay.train, total.concats, by='pair.hashes')
        ##total.concats = dist.decay.test[(dist.a < 50 | dist.b < 50), .(total.concats = sum(num.concats)), by=c('pair.hashes')]
        
###GLM town

     print('training model')
     covariates = c('value.a.ratio','value.b.ratio')
     fmstring = paste('num.concats ~', paste(paste0('log(', covariates, ')', collapse = ' + ')))
        ##fmstring = paste0(fmstring, " + ", "offset(log(total.concats))") this sometimes does 
     fm = formula(fmstring)


        ##browser()
     num.to.sample = 500000
     if(num.to.sample > dim(dist.decay.train)[[1]]) {
         num.to.sample = dim(dist.decay.train)[[1]]
     }
            
     model = glm.nb(formula = fm, data=dist.decay.train[num.concats>0][sample(.N, num.to.sample)], control=glm.control(maxit=500))
     return(model)
}


count_3way_track = function(pairwise.trimmed, dt.concats.sort, all.pairwise, bins) {
    asdf = merge(pairwise.trimmed[i!=j], dt.concats.sort, by.x='i', by.y='binid', allow.cartesian=TRUE)
    asdf.2 = merge(asdf, dt.concats.sort, by.x=c('j','cidi'), by.y=c('binid','cidi'))
    mergerious = merge(asdf.2[, c('id','cidi','i','j')], dt.concats.sort, by='cidi', allow.cartesian=TRUE)
    mergerious = mergerious[binid != i & binid != j]
    dunka = mergerious[, .(num.concats = .N), by=c('id','binid')]
    return(dunka)
}



count_3way_contacts = function(pairwise.trimmed, dt.concats.sort, all.pairwise, bins, resolution=50000, interchromosomal.distance = 1e8) {
    pairwise.trimmed$agg = do.call(Map, c(f = c, pairwise.trimmed[, c('i','j')]))

    asdf = merge(pairwise.trimmed[i!=j], dt.concats.sort, by.x='i', by.y='binid', allow.cartesian=TRUE)
    asdf.2 = merge(asdf, dt.concats.sort, by.x=c('j','cidi'), by.y=c('binid','cidi'))
    mergerious = merge(asdf.2[, c('id','cidi','i','j','agg')], dt.concats.sort, by='cidi')
    mergerious = mergerious[binid != i & binid != j]
    dunka = mergerious[, .(num.concats = .N, agg), by=c('id','i','j','binid')]
    dt.sub = dunka 

    colnames(dunka)[4] = 'V1'
    dt.sub = dunka[, .(sub.bin = unlist(agg)), by=c('id','i','j','num.concats','V1')]
    dt.sub[sub.bin < V1, c('sub.bin','V1') := .(V1, sub.bin)]

    ##dt.yeet = dt.sub[1:100]
    ##dt.yeet[sub.bin < V1, c('sub.bin','V1') := .(V1, sub.bin)]

    ##colnames(dt.sub)[5] = 'pair.value'


    dt.sub = merge(dt.sub, all.pairwise[, c('i','j','pair.value')], by.x=c('V1','sub.bin'), by.y=c('i','j'), all.x=TRUE)

    dt.sub[is.na(pair.value), pair.value := 0]

    dt.sub = dt.sub[V1 != sub.bin]
    ##dt.sub[, sum.pairwise.contact := sum(value), by=c('pair.hashes','V1')]

    ##dt.sub[10000:10100][, .(sum.pairwise.contact = sum(value)), by=c('pair.hashes','V1')]


    ##dt.sub.backup = dt.sub %>% copy
    ##dt.sub = dt.sub.backup %>% copy

    ##dt.sub = merge(dt.sub, pairwise[, c('i','j','pair.hashes','pair.value')], by='pair.hashes', all.x=TRUE)
    dt.sub$V1.isinter = !(dt.sub[, V1 == i] | dt.sub[, V1 == j])
    dt.sub$sub.isinter = !(dt.sub[, sub.bin == i] | dt.sub[, sub.bin == j])


    dt.sub$binterrogate = 0

    ##dt.sub[1:100][V1.isinter==TRUE, binter
    ##dt.subbington = dt.sub[1:100]

    ##dt.subbington[, binterrogate := NULL]
    ##dt.subbington$binterrogate = 0


    dt.sub[V1.isinter==TRUE, binterrogate := V1]
    dt.sub[sub.isinter==TRUE, binterrogate := sub.bin]
    ##dt.yeet[!is.na(binterrogate)]

    dt.sub = dt.sub[binterrogate!=0]


    ##colnames(dt.sub)[6] = 'pair.value'
    dt.sub$pair.hashes = dt.sub$id
    dt.sub[, sum.pairwise.contacts := sum(pair.value), by=c('pair.hashes','binterrogate')]

    
    dt.sub[, dist.i := abs(binterrogate - i)]
    dt.sub[, dist.j := abs(binterrogate - j)]

    bins.gr = bins
    
    bins.gr$binid = 1:length(bins.gr)
    bins.dt = gr2dt(bins.gr)
    setkey(bins.dt, 'binid')

    dt.sub$chr.i = bins.dt[dt.sub$i]$seqnames
    dt.sub$chr.j = bins.dt[dt.sub$j]$seqnames
    dt.sub$chr.binterrogate = bins.dt[dt.sub$binterrogate]$seqnames

## ###interchromosomal distance

    inter.dist = interchromosomal.distance / resolution
    dt.sub[chr.i != chr.binterrogate, dist.i := inter.dist]
    dt.sub[chr.j != chr.binterrogate, dist.j := inter.dist]


    dt.sub[, chr.i := NULL]
    dt.sub[, chr.j := NULL]
    dt.sub[, chr.V1 := NULL]


    
    dt.sub[, diff := j-i]
    dt.sub = dt.sub[diff > 1]

    
    dt.sub[, pair.value := pair.value + 1] ###for log covariate purposes

    print('calculating close and far pairwise contact values')
    
    dt.sub[binterrogate < i & j == sub.bin, value.b := pair.value] ##CLOSER BIN
    dt.sub[binterrogate < i & i == sub.bin, value.a := pair.value] ##FARTHER BIN

    dt.sub[binterrogate > j & i == V1, value.a := pair.value]
    dt.sub[binterrogate > j & j == V1, value.b := pair.value]

##dt.sub[binterrogate < i & j == sub.bin, value..b := value * abs(j - binterrogate)] ##CLOSER BIN
##dt.sub[binterrogate < i & i == sub.bin, value..a := value * abs(i - binterrogate)] ##FARTHER BIN

##dt.sub[binterrogate > j & i == V1, value.ratio.a := value * abs(i - binterrogate)] ##FARTHER BIN
##dt.sub[binterrogate > j & j == V1, value.ratio.b := value * abs(j - binterrogate)] ##CLOSER BIN

    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) < (j - binterrogate)) & i==V1, value.a := pair.value]
    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) < (j - binterrogate)) & j==sub.bin, value.b := pair.value]


    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) >= (j - binterrogate)) & j==sub.bin, value.a := pair.value]
    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) >= (j - binterrogate)) & i==V1, value.b := pair.value]





    

    dt.sub.unique = unique(dt.sub, by=c('pair.hashes','binterrogate'))
    ##dt.sub.unique[, value.ratio.a := NULL]
    ##dt.sub.unique[, value.ratio.b := NULL]





    ##a.ratio = dt.sub[!is.na(value.ratio.a), c('pair.hashes','value.ratio.a','binterrogate')]
    ##b.ratio = dt.sub[!is.na(value.ratio.b), c('pair.hashes','value.ratio.b','binterrogate')]

    ##dt.sub.unique = merge(dt.sub.unique, a.ratio, by=c('pair.hashes','binterrogate'))
    ##dt.sub.unique = merge(dt.sub.unique, b.ratio, by=c('pair.hashes','binterrogate'))




    dt.sub.unique[, value.a := NULL]
    dt.sub.unique[, value.b := NULL]

    a.value = dt.sub[!is.na(value.a), c('pair.hashes','value.a','binterrogate')] %>% unique(by=c('pair.hashes','binterrogate'))
    b.value = dt.sub[!is.na(value.b), c('pair.hashes','value.b','binterrogate')] %>% unique(by=c('pair.hashes','binterrogate'))

    dt.sub.unique = merge(dt.sub.unique, a.value, by=c('pair.hashes','binterrogate'))
    dt.sub.unique = merge(dt.sub.unique, b.value, by=c('pair.hashes','binterrogate'))

    dt.sub.unique[, dist.a := min(dist.i, dist.j), by=c('pair.hashes','binterrogate')]
    dt.sub.unique[, dist.b := max(dist.i, dist.j), by=c('pair.hashes','binterrogate')]

    ##dt.sub.agg = dt.sub[, .(value.ratio.b = value.ratio[1], value.ratio.b = value.ratio[1]), by=c('pair.hashes','binterrogate')]
    ##dt.sub.agg[, value.ratio.a := unlist(V1)[[1]]]
    ##dt.sub.agg[, value.ratio.b := unlist(V1)[2]]

    ##dt.sub.unique[, value.ratio.a := NULL]
    ##dt.sub.unique[, value.ratio.b := NULL]
    
    dt.small = dt.sub.unique[, c('pair.hashes','num.concats','pair.value','binterrogate','sum.pairwise.contacts','dist.a','dist.b','value.a','value.b')]
    dt.small = unique(dt.small, by=c('pair.hashes','binterrogate'))


    
    ##dt.small[dist.a != 1 & dist.b != 1] ###remove bins directly adjacent to 

    dt.small[, value.a.ratio := value.a / dist.a]
    dt.small[, value.b.ratio := value.b / dist.b]
    dt.small = dt.small[dist.a > 1 & dist.b > 1]

    dt.small = merge(dt.small, pairwise.trimmed[, c('i','j','id')], by.x='pair.hashes', by.y='id')
    ##saveRDS(dt.small, 'dt_small_25k.rds')
    ##saveRDS(simplicial.complex, 'weighted_simplicial_complexes/simplicial_complex_chr8_25k.rds')
    ##output = list(dt.small, combn.city, pairwise)
    ##return(list(dt.small, combn.city, pairwise))
    ##saveRDS(dt.small, 'dist_decay_windows_training_sampled_faraway.rds')
    ##saveRDS(dt.small, 'dist_decay_windows_training_test_chr8.rds')

    return(dt.small)
}






####MADE changes and broke it for some reason
## ###FOR INTERCHR
## ###PARALLELIZABLE
## count_3way_contacts = function(pairwise.trimmed, dt.concats.sort, all.pairwise, bins) {
##     pairwise.trimmed$agg = do.call(Map, c(f = c, pairwise.trimmed[, c('i','j')]))

##     asdf = merge(pairwise.trimmed[i!=j], dt.concats.sort, by.x='i', by.y='binid', allow.cartesian=TRUE)
##     asdf.2 = merge(asdf, dt.concats.sort, by.x=c('j','cidi'), by.y=c('binid','cidi'))
##     mergerious = merge(asdf.2[, c('id','cidi','i','j','agg')], dt.concats.sort, by='cidi', allow.cartesian=TRUE)
##     mergerious = mergerious[binid != i & binid != j]
##     dunka = mergerious[, .(num.concats = .N, agg), by=c('id','i','j','binid')]
##     dt.sub = dunka 
## x
##     colnames(dunka)[4] = 'V1'
##     dunka = unique(dunka, by=c('id','V1'))
##     dt.sub = dunka[, .(sub.bin = unlist(agg)), by=c('id','i','j','num.concats','V1')]
##     dt.sub[sub.bin < V1, c('sub.bin','V1') := .(V1, sub.bin)]

##     ##dt.yeet = dt.sub[1:100]
##     ##dt.yeet[sub.bin < V1, c('sub.bin','V1') := .(V1, sub.bin)]

##     ##colnames(dt.sub)[5] = 'pair.value'


##     dt.sub = merge(dt.sub, all.pairwise[, c('i','j','pair.value')], by.x=c('V1','sub.bin'), by.y=c('i','j'), all.x=TRUE)

##     dt.sub[is.na(pair.value), pair.value := 0]

##     dt.sub = dt.sub[V1 != sub.bin]
##     ##dt.sub[, sum.pairwise.contact := sum(value), by=c('pair.hashes','V1')]

##     ##dt.sub[10000:10100][, .(sum.pairwise.contact = sum(value)), by=c('pair.hashes','V1')]


##     ##dt.sub.backup = dt.sub %>% copy
##     ##dt.sub = dt.sub.backup %>% copy

##     ##dt.sub = merge(dt.sub, pairwise[, c('i','j','pair.hashes','pair.value')], by='pair.hashes', all.x=TRUE)
##     dt.sub$V1.isinter = !(dt.sub[, V1 == i] | dt.sub[, V1 == j])
##     dt.sub$sub.isinter = !(dt.sub[, sub.bin == i] | dt.sub[, sub.bin == j])


##     dt.sub$binterrogate = 0

##     ##dt.sub[1:100][V1.isinter==TRUE, binter
##     ##dt.subbington = dt.sub[1:100]

##     ##dt.subbington[, binterrogate := NULL]
##     ##dt.subbington$binterrogate = 0


##     dt.sub[V1.isinter==TRUE, binterrogate := V1]
##     dt.sub[sub.isinter==TRUE, binterrogate := sub.bin]
##     ##dt.yeet[!is.na(binterrogate)]

##     dt.sub = dt.sub[binterrogate!=0]


##     ##colnames(dt.sub)[6] = 'pair.value'
##     dt.sub$pair.hashes = dt.sub$id
##     dt.sub[, sum.pairwise.contacts := sum(pair.value), by=c('pair.hashes','binterrogate')]

    
##     dt.sub[, dist.i := abs(binterrogate - i)]
##     dt.sub[, dist.j := abs(binterrogate - j)]

##     bins.gr = bins
    
##     bins.gr$binid = 1:length(bins.gr)
##     bins.dt = gr2dt(bins.gr)
##     setkey(bins.dt, 'binid')

##     dt.sub$chr.i = bins.dt[dt.sub$i]$seqnames
##     dt.sub$chr.j = bins.dt[dt.sub$j]$seqnames
##     dt.sub$chr.binterrogate = bins.dt[dt.sub$binterrogate]$seqnames

## ## ###interchromosomal distance

##     dt.sub[chr.i != chr.binterrogate, dist.i := 2000]
##     dt.sub[chr.j != chr.binterrogate, dist.j := 2000]


##     dt.sub[, chr.i := NULL]
##     dt.sub[, chr.j := NULL]
##     dt.sub[, chr.V1 := NULL]


    
##     dt.sub[, diff := j-i]
##     dt.sub = dt.sub[diff > 1]

    
##     dt.sub[, pair.value := pair.value + 1] ###for log covariate purposes

##     print('calculating close and far pairwise contact values')
    
##     dt.sub[binterrogate < i & j == sub.bin, value.b := pair.value] ##CLOSER BIN
##     dt.sub[binterrogate < i & i == sub.bin, value.a := pair.value] ##FARTHER BIN

##     dt.sub[binterrogate > j & i == V1, value.a := pair.value]
##     dt.sub[binterrogate > j & j == V1, value.b := pair.value]

## ##dt.sub[binterrogate < i & j == sub.bin, value..b := value * abs(j - binterrogate)] ##CLOSER BIN
## ##dt.sub[binterrogate < i & i == sub.bin, value..a := value * abs(i - binterrogate)] ##FARTHER BIN

## ##dt.sub[binterrogate > j & i == V1, value.ratio.a := value * abs(i - binterrogate)] ##FARTHER BIN
## ##dt.sub[binterrogate > j & j == V1, value.ratio.b := value * abs(j - binterrogate)] ##CLOSER BIN

##     dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) < (j - binterrogate)) & i==V1, value.a := pair.value]
##     dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) < (j - binterrogate)) & j==sub.bin, value.b := pair.value]


##     dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) >= (j - binterrogate)) & j==sub.bin, value.a := pair.value]
##     dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) >= (j - binterrogate)) & i==V1, value.b := pair.value]





    

##     dt.sub.unique = unique(dt.sub, by=c('pair.hashes','binterrogate'))
##     ##dt.sub.unique[, value.ratio.a := NULL]
##     ##dt.sub.unique[, value.ratio.b := NULL]





##     ##a.ratio = dt.sub[!is.na(value.ratio.a), c('pair.hashes','value.ratio.a','binterrogate')]
##     ##b.ratio = dt.sub[!is.na(value.ratio.b), c('pair.hashes','value.ratio.b','binterrogate')]

##     ##dt.sub.unique = merge(dt.sub.unique, a.ratio, by=c('pair.hashes','binterrogate'))
##     ##dt.sub.unique = merge(dt.sub.unique, b.ratio, by=c('pair.hashes','binterrogate'))




##     dt.sub.unique[, value.a := NULL]
##     dt.sub.unique[, value.b := NULL]

##     a.value = dt.sub[!is.na(value.a), c('pair.hashes','value.a','binterrogate')] %>% unique(by=c('pair.hashes','binterrogate'))
##     b.value = dt.sub[!is.na(value.b), c('pair.hashes','value.b','binterrogate')] %>% unique(by=c('pair.hashes','binterrogate'))

##     dt.sub.unique = merge(dt.sub.unique, a.value, by=c('pair.hashes','binterrogate'))
##     dt.sub.unique = merge(dt.sub.unique, b.value, by=c('pair.hashes','binterrogate'))

##     dt.sub.unique[, dist.a := min(dist.i, dist.j), by=c('pair.hashes','binterrogate')]
##     dt.sub.unique[, dist.b := max(dist.i, dist.j), by=c('pair.hashes','binterrogate')]

##     ##dt.sub.agg = dt.sub[, .(value.ratio.b = value.ratio[1], value.ratio.b = value.ratio[1]), by=c('pair.hashes','binterrogate')]
##     ##dt.sub.agg[, value.ratio.a := unlist(V1)[[1]]]
##     ##dt.sub.agg[, value.ratio.b := unlist(V1)[2]]

##     ##dt.sub.unique[, value.ratio.a := NULL]
##     ##dt.sub.unique[, value.ratio.b := NULL]
    
##     dt.small = dt.sub.unique[, c('pair.hashes','num.concats','pair.value','binterrogate','sum.pairwise.contacts','dist.a','dist.b','value.a','value.b')]
##     dt.small = unique(dt.small, by=c('pair.hashes','binterrogate'))


    
##     ##dt.small[dist.a != 1 & dist.b != 1] ###remove bins directly adjacent to 

##     dt.small[, value.a.ratio := value.a / dist.a]
##     dt.small[, value.b.ratio := value.b / dist.b]
##     dt.small = dt.small[dist.a > 1 & dist.b > 1]

##     dt.small = merge(dt.small, pairwise.trimmed[, c('i','j','id')], by.x='pair.hashes', by.y='id')
##     ##saveRDS(dt.small, 'dt_small_25k.rds')
##     ##saveRDS(simplicial.complex, 'weighted_simplicial_complexes/simplicial_complex_chr8_25k.rds')
##     ##output = list(dt.small, combn.city, pairwise)
##     ##return(list(dt.small, combn.city, pairwise))
##     ##saveRDS(dt.small, 'dist_decay_windows_training_sampled_faraway.rds')
##     ##saveRDS(dt.small, 'dist_decay_windows_training_test_chr8.rds')

##     return(dt.small)
## }


create_bin_pair_network_efficient = function(all.pairwise, genome.wide.dist.decay, rr.thresh=0) {


    all.pairwise.trimmed.2 = all.pairwise %>% copy
    all.pairwise.trimmed.2$i = all.pairwise$j
    all.pairwise.trimmed.2$j = all.pairwise$i
    pairwise.sym = rbind(all.pairwise, all.pairwise.trimmed.2)

    ##all.pairwise.trimmed.2 = pairwise.trimmed %>% copy
    ##all.pairwise.trimmed.2$i = pairwise.trimmed$j
    ##all.pairwise.trimmed.2$j = pairwise.trimmed$i
    ##pairwise.sym = rbind(pairwise.trimmed, all.pairwise.trimmed.2)


    neighbor.1 = merge(genome.wide.dist.decay[relative.risk >= rr.thresh & pval < 0.05, c('i','j','binterrogate','relative.risk','pair.hashes')], pairwise.sym, by.x=c('i','binterrogate'), by.y=c('i','j'))
    neighbor.2 = merge(genome.wide.dist.decay[relative.risk >= rr.thresh & pval < 0.05, c('i','j','binterrogate','relative.risk','pair.hashes')], pairwise.sym, by.x=c('j','binterrogate'), by.y=c('i','j'))




    all.neighbors = rbind(neighbor.1, neighbor.2)


##all.neighbors$sidi = do.call(Map, c(f = c, all.neighbors[, c('i','j','binterrogate')]))
##all.neighbors$sidi = map(all.neighbors$sidi, sort)
##all.neighbors$hashes = all.neighbors$sidi %>% unname %>%  paste0
##all.neighbors$hashes = map(all.neighbors$sidi, unname) %>% paste0
##browser()
    ##all.neighbors = unique(all.neighbors, by=c('i','j','binterrogate'))

    all.neighbors$pair.hashes.x = all.neighbors$pair.hashes.x %>% factor
    all.neighbors$pair.hashes.y = all.neighbors$pair.hashes.y %>% factor
    G.kant.10 = graph_from_edgelist(as.matrix(all.neighbors[,c('pair.hashes.x','pair.hashes.y')]), directed=TRUE)
    E(G.kant.10)$weight = all.neighbors$relative.risk
    return(G.kant.10)
}


run_interchr_dist_decay = function(concatemers, folder, bins=NULL, pair.thresh=50, chromosome=NULL, training.chr='chr8', resolution=50000, rr.thresh=0, fdr.thresh=0.1, pval.thresh=0.05, expansion.cutoff=0.5, num.members=10) {

    if(is.null(chromosome)){
        chromosome = c(paste0("chr", c(as.character(1:22), "X")))
    }

    concatemers = concatemers %Q% (seqnames %in% chromosome)
    
    unique_read_idx = concatemers$read_idx %>% unique
    training.idx = sample(unique_read_idx, length(unique_read_idx)/2)
    training.gr = concatemers %Q% (read_idx %in% training.idx)
    test.gr = concatemers %Q% !(read_idx %in% training.idx)

    chrom = interchr_dist_decay_binsets(training.gr, bins=bins, rr.thresh=rr.thresh, pair.thresh=pair.thresh, folder=folder, training.chr=training.chr, resolution=resolution, fdr.thresh=fdr.thresh, pval.thresh=pval.thresh, expansion.cutoff=expansion.cutoff, num.members=num.members)
    
    dir.create(folder)
    saveRDS(chrom, paste0(folder, '/chromunity_results.rds'))
    saveRDS(training.gr, paste0(folder, '/training_concatemers.rds'))

    chid.to.test = (chrom$chid %>% unique)
    syn.interchr = evaluate_synergy_experimental(chrom, test.gr, chid.to.test, folder=folder, chromosome=chromosome)
    syns = aggregate_synergy_results_singledir(folder, strict.check=TRUE)
    return(syns)
}


aggregate_synergy_results_singledir = function(dir, strict.check=FALSE, fdr.thresh=0.1) {

    synergy.chunk = readRDS(paste0(dir, '/synergy_results.rds'))
    synergy.chunk = synergy.chunk[!is.na(fdr)]
    binsets = readRDS(paste0(dir, '/binsets.rds'))

    annotate = readRDS(paste0(dir, '/chrom_annotate.rds'))
    if(strict.check)
        synergy.chunk = symmetry_check(annotate, synergy.chunk, binsets)

    binsets$bid = factor(binsets$bid)
    synsets = merge(binsets, synergy.chunk[annotation=='chromunity', c('bid','fdr')], by='bid')
    synsets = synsets[fdr<fdr.thresh]
    return(list(synergy.chunk, dt2gr(synsets)))
}



calculate_binned_contact_count = function(concats.gr, bins.gr) {
    binned.concats = bin_concatemers(concats.gr, bins = bins.gr, max.slice = 1e6, mc.cores=5, verbose=TRUxE, hyperedge.thresh = 3)
    dedupe.concats = unique(binned.concats, by=c('read_idx','binid'))
    dedupe.concats[, numbins := .N, by='read_idx']
    dedupe.concats[, num.contacts := choose(numbins, 2)]
    unique.concats = unique(dedupe.concats, by='read_idx')
    binned.contact.count = sum(unique.concats$num.contacts)
    return(binned.contact.count)
}

shuffle_concatemers_spiked = function(concatemers, contact_matrix, spike.set, strength=0.1, noise.parameter = 10000, resolution=10000, trim=TRUE, interval.width=1000) {
    A = contact_matrix$mat %>% as.matrix
    rownames(A) <- NULL
    colnames(A) <- NULL
    A[cbind(1:nrow(A), 1:nrow(A))] = 0
    A = A + t(A)
    An = A 
    ##browser()
    An = round(1+10*An/min(An[An>0]))  ##can you remove this background 1 value? does this change anything. peculiar     An = round(1+10*An/min(An[An>0]))
    edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
     ##
    G = graph.edgelist(edges[, cbind(row, col)])
    ##G.2 = graph.adjacency(An)#, directed=FALSE)
    
##    browser()
    concats.binned = bin_concatemers(concatemers, contact_matrix$gr)
    concats.dt = concats.binned
    concats.dt = unique(concats.dt, by=c('binid','cidi')) ###dedupe!!

    concat.counts = concats.dt[, new.count := .N, by='read_idx']


    card = unique(concat.counts[, .(read_idx, new.count)])

    this.steps = sum(card$new.count)

    row.labels = 1:dim(An)[1]
    start.index = row.labels[rowSums(An)>1]
    
    RW = random_walk(G, start = start.index, steps = sum(card$new.count)) %>% as.numeric
    ##browser()
    ##
    ##rm(G)
    ##rm(edges)
    ##gc()

    ###find set of concatemers which already include one of the spike set 
    ##spike.set = c(7200, 7210, 7220, 7230)

    out = contact_matrix$gr[RW] %>% gr2dt()
    ##browser()
    
    out = out[1:sum(card$new.count)]
    out$read_idx = card[, rep(read_idx, new.count)] 

    out[, cardinality := .N, by='read_idx']
##    out = gr2dt(dt2gr(out) - ()) ###Setting width of each monomer to be arbitrarily 1000
    ##browser()
    
    ##spike.concatemers = spike.set
    ##spikable.concats = out[cardinality > 1 & binid %in% spike.set$binid]$read_idx %>% unique  ##no concatemer which only hit a bin-set
    ##concats.to.spike = sample(spikable.concats, number.to.add)

    ##browser()



    ###Split by whatever random binsets you want to add here

    dt = gr2dt(spike.set %&% contact_matrix$gr)[, count := .N, by='bid']
    spike.set = dt2gr(dt[count > 2])
    
    spike.set$binid = gr.match(spike.set, contact_matrix$gr)
    spike.set$chid = spike.set$chid %>% as.character
    binsets.to.add = split(spike.set, spike.set$chid)
    out[, c('query.id','tile.id') := NULL]

###We are going to choose random position within bin to begin monomer from
###And then create fixed width monomer
    
    out$rowid = 1:dim(out)[[1]]
    out[, new.start.interval := sample(resolution - interval.width, 1), by='rowid']
    out[, start := start + new.start.interval]
    out[, end := start + interval.width]
    

    ##spike.set.gr = binsets.to.add$rand_15

    spikes = pbmclapply(binsets.to.add, mc.cores = 5, function(spike.set.gr) {

        ##number.to.add = choose(spike.set.gr[1]$cardinality, 3) * strength
##        spike.set.gr = binsets.to.add$rand_10

        binset.rand = spike.set.gr %>% gr2dt
        binset.rand$bid=1
        binset.rand$index = 1:dim(binset.rand)[[1]]
        power = binset.rand[, powerset(index, 3, 3), by=bid]

        
        num.concats = merge(binset.rand, out, by='binid')$read_idx %>% unique %>% length
        number.to.add = num.concats * strength
        print(number.to.add)
        
        power[, size := .N, by=c('setid','bid')]
        size3 = unique(power[size==3], by=c('bid','setid'))
        size3.samp = sample(size3$setid, number.to.add, replace=TRUE)
        samp = data.table(size3.samp)
        colnames(samp) = 'setid'
        samp.binned = merge(samp, power, by='setid', allow.cartesian=TRUE)
        spike.concats = binset.rand[samp.binned$item]


        

        ###choose spike location randomly from bin-set given
        spike.concats$rowid = 1:dim(spike.concats)[[1]]
        spike.concats[, new.start.interval := sample(resolution - interval.width, 1), by='rowid']
        spike.concats[, start := start + new.start.interval]
        spike.concats[, end := start + interval.width]

        spike.concats.gr = dt2gr(spike.concats)

        dt.spike = gr2dt(spike.concats.gr)
        dt.spike$spike.id = rep(1:dim(samp)[[1]], each=3)
        dt.spike.2 = merge(out, dt.spike, by='binid', allow.cartesian=TRUE)
        sample.spike = dt.spike.2[, .(read_idx = sample(read_idx, 5)), by='spike.id'] ###almost guarantees we don't double select a concatemer
        sample.spike = sample.spike %>% unique(by='read_idx') %>% unique(by='spike.id')
        sample.spike = merge(sample.spike, dt.spike.2[, c('spike.id','read_idx','binid')], by=c('spike.id','read_idx'))
        sample.spike = sample.spike %>% unique(by='read_idx') %>% unique(by='spike.id')
        

        dt.spike[, binid := NULL]
        merge.spike = merge(dt.spike, sample.spike, by='spike.id', allow.cartesian=TRUE)
        merge.spike$cardinality = 3
        
        dt = merge.spike[, c('seqnames','start','end','strand','width','binid','read_idx','cardinality')]


        jitter = rnorm(dim(dt)[[1]], sd=noise.parameter)
        dt$jitter = jitter
        dt[, start := start + jitter]
        dt[, end := end + jitter]
        dt[, jitter := NULL]
        return(dt)
    })
    spikes = spikes %>% rbindlist

####Remove the monomer which allowed concatemer to become spikable in the first place...

    
    out$rowid = 1:dim(out)[[1]]
    out.to.drop = merge(unique(spikes, by='read_idx'), out, by=c('read_idx','binid'))
    out = out[!(rowid %in% out.to.drop$rowid)]
    
    
    ##browser()
    
    spikes$binid = gr.match(dt2gr(spikes), contact_matrix$gr)
    
    spike.out = rbind(out, spikes, fill=TRUE)

    spike.out$cid = spike.out$read_idx
    spike.out[, chid := NULL]

    if(trim == TRUE){
        spiked.concats = spikes[read_idx %in% spikes$read_idx]
        added.pairwise.contacts = cocount(dt2gr(spiked.concats), bins=contact_matrix$gr, by='read_idx')
        dt = added.pairwise.contacts$dat
        to.remove = dt[(i %in% spikes$binid) & (j %in% spikes$binid)][i != j]

        dt.i = merge(out[, c('read_idx','binid')], to.remove, by.x='binid', by.y='i', allow.cartesian=TRUE)
        dt.j = merge(out[, c('read_idx','binid')], dt.i, by.x = c('binid','read_idx'), by.y=c('j','read_idx'))

        dt.j[, count := 1:.N, by='id']
        cid.bins = dt.j[binid != binid.y]
        cid.bins = cid.bins[count <= value]
        
        cid.bins.2 = cid.bins %>% copy

        cid.bins.2$binid = cid.bins$binid.y
        cid.bins.2$binid.y = cid.bins$binid

        dt.bins = rbind(cid.bins, cid.bins.2)[, c('read_idx','binid')] %>% unique(by=c('read_idx','binid'))

        ##browser()
        out$rowid = 1:dim(out)[[1]]
        remove.out = merge(out, dt.bins, by=c('read_idx','binid'))
        remove.out = unique(remove.out, by='read_idx')
        out.trimmed = out[!(rowid %in% remove.out$rowid)]
        out.trimmed[, rowid := NULL]
        
        spike.out.trimmed = rbind(out.trimmed, spikes, fill=TRUE)
        spike.out = spike.out.trimmed
    }
    spike.out[, end := start + interval.width] 
    return(spike.out)
}
