library(data.table)
library(reader)
datadir = '/home/panwei/lin00374/cML/graph/vcffiles'
savedir = '/home/panwei/lin00374/ND/revision1/result'
pre1_list = pre2_list = c('ebi-a-GCST005195','ieu-a-89','ebi-a-GCST002222','ebi-a-GCST90018970')


G = G_pval =  matrix(0,nrow=length(pre1_list),ncol=length(pre1_list),dimnames=list(pre1_list,pre2_list))
pv_rg_mat = NULL
rg_jk_mat = NULL
for(i in 1:(length(pre1_list)-1)){
    for(j in (i+1):length(pre1_list)){
        if(j>length(pre1_list)){next;}
        pre1 = pre1_list[i]
        pre2 = pre1_list[j]
        ldsc_cmd = sprintf('source activate ldsc; /home/panwei/lin00374/ldsc/ldsc.py  --rg %s/%s.vcf.gz.txt.sumstats.gz,%s/%s.vcf.gz.txt.sumstats.gz --ref-ld-chr /home/panwei/lin00374/ldsc/eur_w_ld_chr/  --w-ld-chr /home/panwei/lin00374/ldsc/eur_w_ld_chr/ --print-delete-vals --out %s/%s~%s',datadir,pre1,datadir,pre2,savedir,pre1,pre2)
        if(file.exists(sprintf('%s/%s~%s.log',savedir,pre1,pre2))){
            a = n.readLines(sprintf('%s/%s~%s.log',savedir,pre1,pre2),n=100,header=FALSE)
        }else{
            a = system(ldsc_cmd,intern=T)
        }
        rg = as.numeric(unlist(strsplit(x=a[66],split= '\\s+'))[3])
        rg_pval =  as.numeric(unlist(strsplit(x=a[66],split= '\\s+'))[6])
        h2_i = unlist(strsplit(x=a[60],split= '\\s+'))[7]
        h2_i = read.table(substr(h2_i, 1, nchar(h2_i) - 1))
        h2_j = unlist(strsplit(x=a[61],split= '\\s+'))[7]
        h2_j = read.table(substr(h2_j, 1, nchar(h2_j) - 1))
        gcov = unlist(strsplit(x=a[62],split= '\\s+'))[7]
        gcov = read.table(substr(gcov, 1, nchar(gcov) - 1))
        rg_jk = (gcov/sqrt(h2_i)/sqrt(h2_j))[,1]
        rg_jk_mat = cbind(rg_jk_mat,rg_jk)
        n_block = length(rg_jk)
        pv_rg = n_block * rg - (n_block-1) * rg_jk
        pv_rg_mat = cbind(pv_rg_mat,pv_rg)
        sd_rg = sqrt(var(pv_rg)/n_block)
        G[pre1,pre2] = as.numeric(rg)
        G_pval[pre1,pre2] = as.numeric(rg_pval)
    }
}


cov_rg = cov(pv_rg_mat)/n_block
G[lower.tri(G)] = t(G)[lower.tri(G)]
est_rg = G[lower.tri(G)]
B = 500
G_dir_list = G_obs_list = vector(mode='list', length=B)
set.seed(123)
for(i in 1:B){
    G_obs = matrix(0,nrow=nrow(G),ncol=ncol(G))
    G_obs[lower.tri(G_obs)] = MASS::mvrnorm(n=1, mu = est_rg, Sigma=cov_rg)
    G_obs[upper.tri(G_obs)]  = t(G_obs)[upper.tri(G_obs)]
    G_obs_list[[i]] = G_obs
    G_dir_list[[i]] = G_obs %*% solve(diag(nrow(G_obs))+G_obs);
}
ldsc_dir_graph_mean = apply(simplify2array(G_dir_list), 1:2, mean)
ldsc_dir_graph_sd = apply(simplify2array(G_dir_list), 1:2, sd)
ldsc_dir_graph_pval = pnorm(-abs(ldsc_dir_graph_mean/ldsc_dir_graph_sd))*2
ldsc_obs_graph_mean = apply(simplify2array(G_obs_list), 1:2, mean)
ldsc_obs_graph_sd = apply(simplify2array(G_obs_list), 1:2, sd)
ldsc_obs_graph_pval = pnorm(-abs(ldsc_obs_graph_mean/ldsc_obs_graph_sd))*2
colnames(ldsc_obs_graph_pval) = rownames(ldsc_obs_graph_pval) = pre1_list
colnames(ldsc_obs_graph_mean) = rownames(ldsc_obs_graph_mean) = pre1_list
colnames(ldsc_dir_graph_mean) = rownames(ldsc_dir_graph_mean) = pre1_list
colnames(ldsc_dir_graph_pval) = rownames(ldsc_dir_graph_pval) = pre1_list


saveRDS(list(G_obs = ldsc_obs_graph_mean, pval_obs = ldsc_obs_graph_pval,
             G_dir = ldsc_dir_graph_mean, pval_dir = ldsc_dir_graph_pval
            ),'/home/panwei/lin00374/ND/Result_4traits_ldsc.rds')
