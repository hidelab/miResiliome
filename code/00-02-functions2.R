tarbase.prepare <- function(TarBaseFile = NULL){
    
    tarbase   <- read.csv(Tar.dir, sep = "\t")
    mir.sets  <- tarbase %>% filter(.,species == "Homo sapiens",
                                    positive_negative == "POSITIVE" ) 
    mir.sets  <- mir.sets %>% dplyr::select(.,mirna,geneId) %>% group_split(.,mirna) 
    mir.list  <- sapply(mir.sets, function(X){dplyr::slice(X,1) %>% 
            dplyr::select(.,mirna) %>% pull})
    mir.sets  <- sapply(mir.sets,function(X){dplyr::select(X,geneId) %>% unique})
    
    names(mir.sets) <- mir.list
    return(mir.sets)
}


tarbase.prepare2 <- function(TarBaseFile = NULL){
    
    tarbase   <- read.csv(Tar.dir, sep = "\t")
    mir.sets  <- tarbase %>% filter(.,species == "Homo sapiens",
                                    positive_negative == "POSITIVE" ) 
    
    mir.sets$geneName  <- gsub("\\(hsa\\)","",mir.sets$geneName)
    mir.sets           <- mir.sets %>% dplyr::select(.,geneName,mirna)  
    colnames(mir.sets) <- c("Target.gene", "Representative.miRNA")
    return(mir.sets)
}




target.scan.prepare <- function(tar.path = NULL){
    target.scan.tab <- list()
    for(i in list.files(tar.path,full.names = T)){
        tab <- read.csv(i,sep = "\t")
        tab <- tab[,c("Target.gene","Representative.miRNA")]
        target.scan.tab <- rbind(target.scan.tab,tab)
    }
    return(target.scan.tab)
}

lin.path.model <- function(Y = NA, X1 = NA , DAT = NA, mir.name = NA){
    
    formz        <- paste0(eval(Y), " ~ " , eval(X1))
    out.lm       <- lm(as.formula(formz)  , data = DAT)
    out.lm       <- tidy(out.lm)
    out.lm$transcript <- mir.name
    return(out.lm)
}


cal.cor.tab   <- function(mir.vec = NA, exp.dat = NA){
   
     mir.mr.cor  <- apply(exp.dat, 1,function(X){   
                                cors <- cor.test(mir.vec,X, method = "spearman")
                                return(c(cors$estimate,cors$p.value))})
                            
    mir.mr.cor        <-  as.data.frame(t(mir.mr.cor))
    mir.mr.cor$p.adj  <-  p.adjust(mir.mr.cor$V2, method = "fdr")
    names(mir.mr.cor) <- c("cor", "pval","p.adj")
    return(mir.mr.cor)
}




proportional.summary <- function(Y = NA, X1 = NA , DAT = NA, mir.name = NA){
    
    formz <- paste0(eval(Y), " ~ " , eval(X1))
    fit   <- vglm(as.formula(formz) ,
                  family=cumulative(parallel=TRUE),
                  data=DAT)
    
    temp <- summaryvglm(fit)
    temp <- temp@coef3
    temp <- as.data.frame(temp)
    temp <- rownames_to_column(temp,var =  "term")
    temp$gene <- mir.name
    temp$R2lat <- R2latvar(fit)
    return(temp)
    
}

proportional.summary2 <- function(Y = NA, X1 = NA , DAT = NA, mir.name = NA){
    
    formz <- paste0(eval(Y), " ~ " , eval(X1))
    fit   <- vglm(as.formula(formz) ,
                  family=cumulative,
                  data=DAT)
    
    temp <- summaryvglm(fit)
    temp <- temp@coef3
    temp <- as.data.frame(temp)
    temp <- rownames_to_column(temp,var =  "term")
    temp$gene <- mir.name
    return(temp)
    
}
clean.write.podds    <- function(a.table, an.address = NA){
    
    
    temp.cog <- a.table[a.table$term =="tested_gene",]
    temp.cog$p.adj <- p.adjust(temp.cog$`Pr(>|z|)`,method = "fdr")
    temp.cog       <- temp.cog[order(temp.cog$`Pr(>|z|)`),]
    temp.cog$qval  <- qvalue::qvalue(temp.cog$`Pr(>|z|)`)$qvalues
    temp.cog$Dir   <- ifelse(temp.cog$Estimate<0, "POS","NEG")
    temp.cog       <- temp.cog %>% mutate_if(., is.numeric,signif,digits = 3)
    temp.cog$pval  <- temp.cog$`Pr(>|z|)`
    
    if(!is.na(an.address))
        write.csv(temp.cog,file = an.address)
    return(temp.cog)
    
}


clean.write.podds2    <- function(a.table, an.address = NA){
    
    
    temp.cog <- a.table
    temp.cog$p.adj <- p.adjust(temp.cog$`Pr(>|z|)`,method = "fdr")
    temp.cog       <- temp.cog[order(temp.cog$`Pr(>|z|)`),]
    temp.cog       <- temp.cog[temp.cog$term =="tested_gene",]
    
    temp.cog$Dir   <- ifelse(temp.cog$Estimate<0, "POS","NEG")
    temp.cog       <- temp.cog %>% mutate_if(., is.numeric,signif,digits = 3)
    if(!is.na(an.address))
        write.csv(temp.cog,file = an.address)
    return(temp.cog)
    
}



consensus_cor <- function(r1, m1, r2, m2, k = 0) {
  z1 <- atanh(r1); z2 <- atanh(r2)
  w1 <- m1 - k - 3; w2 <- m2 - k - 3
  stopifnot(w1 > 0, w2 > 0)
  zhat <- (w1*z1 + w2*z2) / (w1 + w2)
  se   <- 1 / sqrt(w1 + w2)
  Z    <- zhat / se
  p    <- 2*pnorm(-abs(Z))
  ci_z <- zhat + c(-1,1) * qnorm(0.975) * se
  list(
    r_hat = tanh(zhat),
    ci = tanh(ci_z),
    Z = Z, p = p,
    #Cochran's Q statistic for heterogeneity
    Q = w1*(z1-zhat)^2 + w2*(z2-zhat)^2,
    w = c(w1=w1, w2=w2)
  )
}

make_consensus <- function(r1, m1, r2, m2, k = 0) {

    consensus_results <- mapply(
    consensus_cor,
    r1 = r1,
    m1 = m1,
    r2 = r2,
    m2 = m2,
    SIMPLIFY = FALSE
    )
    consensus_results_df <- do.call(rbind, lapply(consensus_results, function(res) {
      data.frame(
        r_hat = res$r_hat,
        ci_lower = res$ci[1],
        ci_upper = res$ci[2],
        Z_consensus = res$Z,
        p_consensus = res$p,
        Q = res$Q,
        stringsAsFactors = FALSE
      )
    }))
    consensus_results_df$p_adj_consensus <- p.adjust(consensus_results_df$p_consensus, method = "BH")
    return(consensus_results_df)

}



diff_cor_test <- function(r1, m1, r2, m2) {
  z1 <- atanh(r1); z2 <- atanh(r2)
  se <- sqrt(1/(m1 - 3) + 1/(m2 - 3))
  Z  <- (z1 - z2)/se
  p  <- 2 * pnorm(-abs(Z))
  c(Z = Z, p = p)
}

make_diff_cor <- function(r1, m1, r2, m2) {
  diff_results <- mapply(
    diff_cor_test,
    r1 = r1,
    m1 = m1,
    r2 = r2,
    m2 = m2,
    SIMPLIFY = TRUE
  )

  diff_results_df <- data.frame(
    Z_diff = diff_results["Z", ],
    deltaR = r1 - r2,
    p_diff = diff_results["p", ],
    P_adj_diff = p.adjust(diff_results["p", ], method = "BH"),
    stringsAsFactors = FALSE
  )
  return(diff_results_df)
}
