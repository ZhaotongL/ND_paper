library(readr)
library(data.table)
library(dplyr)
library(gwasvcf)
library(TwoSampleMR)
library(VariantAnnotation)
source('/home/panwei/lin00374/JMR/helper.R')



data_dir = "/home/panwei/lin00374/cML/graph/vcffiles" 
savedata_dir = "/home/panwei/lin00374/ND/data"

temp_dir <- tempdir()

pre1_list = pre2_list = c('ebi-a-GCST005195','ieu-a-89','ebi-a-GCST002222','4080_irnt')
for(pre1 in pre1_list){
    vcf <- readVcf(sprintf('%s/%s.vcf.gz',data_dir,pre1))
    exp_df = vcf_to_granges(vcf) %>% dplyr::as_tibble() %>% 
        dplyr::select(SNP=ID, N=SS, BETA=ES, SE=SE,A1=ALT, A2=REF, LP,EAF=AF) %>% mutate(P=10^-LP) %>% mutate(id=pre1)
    exp_df = format_data(
                      exp_df,
                      type = "exposure",
                      snp_col = "SNP",
                      beta_col = "BETA",
                      se_col = "SE",
                      eaf_col = "EAF",
                      effect_allele_col = "A1",
                      other_allele_col = "A2",
                      pval_col = "P",
                      samplesize_col = "N",
                      id_col = "id"
                      )

    exp_dat = exp_df %>% filter(pval.exposure<5e-8)
    exp_mr_dat = clump_data(exp_dat)
    
    for(pre2 in pre2_list){
        if(pre1==pre2){next;}
       vcf <- readVcf(sprintf('%s/%s.vcf.gz',data_dir,pre2))
        out_df = vcf_to_granges(vcf) %>% dplyr::as_tibble() %>% 
            dplyr::select(SNP=ID, N=SS, BETA=ES, SE=SE,A1=ALT, A2=REF, LP,EAF=AF,id) %>% mutate(P=10^-LP)
        out_df = format_data(
                      out_df,
                      type = "outcome",
                      snp_col = "SNP",
                      beta_col = "BETA",
                      se_col = "SE",
                      eaf_col = "EAF",
                      effect_allele_col = "A1",
                      other_allele_col = "A2",
                      pval_col = "P",
                      samplesize_col = "N",
                      id_col = "id"
                      )

 
        harmo_dat = harmonise_data(exp_mr_dat,out_df,2)
        ## use steiger filtering for bi-directional MR
        steig_dat = steiger_filtering(harmo_dat)
        harmo_dat = steig_dat %>% filter(mr_keep==TRUE) %>% filter(steiger_dir==TRUE)
        saveRDS(harmo_dat,sprintf('%s/%s~%s.rds',savedata_dir,pre1,pre2))
        }
}

#### 
pre1_list = pre2_list = c('ebi-a-GCST005195','ieu-a-89','ebi-a-GCST002222','4080_irnt')
all_IV = list()
k=1
for(pre1 in pre1_list){
    for(pre2 in pre2_list){
        if(pre1==pre2){next;}
        dat = readRDS(sprintf('%s/%s~%s.rds',savedata_dir,pre1,pre2))
        all_IV[[k]] = dat$SNP
        k = k+1
        }
    }
save(all_IV,file=sprintf('%s/4traits_IV.RData',savedata_dir))

## organize IVs used in all 12 analysis in matrix form
library(dplyr)
all_IV = unique(unlist(all_IV))
all_exposure_dat <- extract_outcome_data(all_IV, pre1_list[-4], proxies = 0)
all_IV_df = data.frame(SNP=all_IV)
all_exposure_dat_tmp = all_exposure_dat %>% dplyr::select(SNP,chr,pos,A1=effect_allele.outcome,A2=other_allele.outcome) %>% unique()
all_IV_df = merge(all_IV_df,all_exposure_dat_tmp,by='SNP',sort=FALSE)
all_IV_df = all_IV_df %>% arrange(as.numeric(chr),as.numeric(pos))
b_mat = se_mat = matrix(NA,nrow=nrow(all_IV_df),ncol=length(pre1_list),dimnames=list(all_IV_df$SNP,pre1_list))
for(pre1 in pre1_list){
    for(pre2 in pre2_list){
        if(pre1==pre2){next;}
        dat = readRDS(sprintf('%s/%s~%s.rds',savedata_dir,pre1,pre2))
        dat.tmp = all_IV_df[match(dat$SNP,all_IV_df$SNP),]
        print(all(dat.tmp$A1==dat$effect_allele.exposure)) ## check allele match - all TRUE!
        b_mat[match(dat$SNP,rownames(b_mat)),pre1] = dat$beta.exposure
        b_mat[match(dat$SNP,rownames(b_mat)),pre2] = dat$beta.outcome
        se_mat[match(dat$SNP,rownames(se_mat)),pre1] = dat$se.exposure
        se_mat[match(dat$SNP,rownames(se_mat)),pre2] = dat$se.outcome
#        k = k+1
        }
}

## prepare LD matrix for IVs in data perturbation
ldBlocks = read.table('/home/panwei/shared/zhaotong_share/ldblocks-eur-hg19.txt',header=TRUE)
all_IV_df$pos = as.numeric(all_IV_df$pos)
all_IV_df$chr = as.numeric(all_IV_df$chr)
matList = list()
matListInd = 1
for(i in 1:nrow(ldBlocks)){
  curChr = substring(ldBlocks[i,1],4)
  curMin = ldBlocks[i,2]
  curMax = ldBlocks[i,3]
  curSnps = all_IV_df[all_IV_df$pos > curMin & all_IV_df$pos <= curMax & all_IV_df$chr == curChr,]
  if(nrow(curSnps) > 0){
    mat_res = list()
    if(nrow(curSnps)==1){
        R_flipped = matrix(1,nrow=1)
        colnames(R_flipped) = rownames(R_flipped) = curSnps$SNP
    }else{
        R = ieugwasr::ld_matrix(variants=curSnps$SNP,bfile=sprintf('/home/panwei/lin00374/1kg.ieugwasr/EUR'),plink_bin='/home/panwei/lin00374/plink')
        R_snp_df = do.call(rbind.data.frame,strsplit(colnames(R),split='_'))
        colnames(R_snp_df) = c('SNP','A1','A2')
        R_mydat_df = merge(R_snp_df,curSnps,by='SNP',sort=F)
        R_toflip_ind = which(with(R_mydat_df,A1.x==A2.y & A2.x==A1.y))
        R_flipped = R
        for(j in R_toflip_ind){
            R_flipped[j,] = -R_flipped[j,]
            R_flipped[,j] = -R_flipped[,j]
        }
        colnames(R_flipped) = rownames(R_flipped) = as.character(R_snp_df$SNP)
    }
    mat_res$snp = curSnps$SNP
    mat_res$R = R_flipped
    matList[[matListInd]] = mat_res
    matListInd = matListInd+1
  }
}
saveRDS(list(b_mat=b_mat,se_mat=se_mat,all_IV=all_IV_df,R_list=matList,P=diag(4),n_vec=c(547261,253288,94595,340159)),sprintf('%s/4traits_data.RData',savedata_dir))