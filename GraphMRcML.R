library(TwoSampleMR)
source('/home/panwei/lin00374/cML/cML_O.R')
data_dir = "/home/panwei/lin00374/ND/data"
setwd(data_dir)
pre1_list = pre2_list = c('ebi-a-GCST005195','ieu-a-89','ebi-a-GCST002222','4080_irnt')

array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

## preprocess data in GraphMRcML format - done! ## 
# dat = readRDS('/home/panwei/lin00374/ND/data/4traits_data.RData')
# n_vec = dat$n_vec
# R_list = dat$R_list
# b_mat = dat$b_mat
# se_mat = dat$se_mat
# rho_mat = dat$P

#   n_trait = length(n_vec)
#   N_combination = n_trait * (n_trait - 1) / 2
#   m_block = length(R_list) ## LD blocks

#   IJ_snp_list = vector("list", N_combination)
#   k = 1
# ### Screening for IVs for each pair of traits ###
#   ## i:trait 1; j:trait 2
#   for(i in 1:(n_trait-1)){
#     for(j in (i+1):n_trait){
#         pre1 = pre1_list[i]
#         pre2 = pre1_list[j]
#       dat = readRDS(sprintf('%s/%s~%s.rds',data_dir,pre1,pre2))
#       IJ_snp_list[[k]]$ind_i_new = which(rownames(b_mat) %in% dat$SNP)
#       dat = readRDS(sprintf('%s/%s~%s.rds',data_dir,pre2,pre1))
#       IJ_snp_list[[k]]$ind_j_new = which(rownames(b_mat) %in% dat$SNP)
#       k = k + 1
#     }
#   }
# ### End of screening 1 ###
# ### Screening for LD of b_mat ###
#   DP_mat_list = vector("list", m_block)
#   for(i in 1:m_block){
#       R = R_list[[i]]$R
#       match_order = match(rownames(b_mat),rownames(R))
#       match_order = match_order[!is.na(match_order)]
#       R = R[match_order,match_order,drop=FALSE]
#       snp_in_block = R_list[[i]]$snp
#       PXR = kronecker(rho_mat,R)
#       eigen_decomp = eigen(PXR)
#       D = diag(sqrt(zapsmall(eigen_decomp$value,digits=6)))
#       DP_mat_list[[i]]$V = eigen_decomp$vector %*% D
#       DP_mat_list[[i]]$snp = rownames(R)
#   }
# ### End of screening 2 ###

#   out = list()
#   out$IJ_snp_list = IJ_snp_list
#   out$DP_mat_list = DP_mat_list
#   out$b_mat = b_mat
#   out$se_mat = se_mat
#     out$rho_mat = rho_mat
#     out$n_vec = n_vec
# screen_res = out
#     saveRDS(screen_res,'/home/panwei/lin00374/ND/data/4traits_screen.rds')

screen_res = readRDS('/home/panwei/lin00374/ND/data/4traits_screen.rds')
Generate_Perturb <- function(b_mat,se_mat,n_vec,rho_mat,DP_mat_list){
  n_trait = length(n_vec)
  m_used = nrow(b_mat)
  e_mat_dp = matrix(0,ncol=n_trait,nrow=m_used)
  for(i in 1:length(DP_mat_list)){
      V = DP_mat_list[[i]]$V
      X = rnorm(nrow(V))
      e_vec = V %*% X
      e_mat = matrix(e_vec,ncol=n_trait,byrow=FALSE)
      snp_ind = which(is.element(rownames(b_mat),DP_mat_list[[i]]$snp))
      e_mat_dp[snp_ind,] = se_mat[snp_ind,,drop=FALSE] * e_mat
  }
  b_mat_dp = b_mat + e_mat_dp

  return(b_mat_dp)
}

Graph_Estimate <- function(b_mat,se_mat,n_vec,rho_mat,IJ_snp_list,t=0,random_start=10){
  n_trait = length(n_vec)
  obs_graph = matrix(1,nrow=n_trait,ncol=n_trait)
  obs_graph_pval = obs_graph_se = matrix(0,nrow=n_trait,ncol=n_trait)
  k = 1
  for(i in 1:(n_trait-1)){
    for(j in (i+1):n_trait){
      ind_i_new = IJ_snp_list[[k]]$ind_i_new
      ind_j_new = IJ_snp_list[[k]]$ind_j_new
      rho_ij = rho_mat[i,j]
      ItoJ_cML_O_res = mr_cML_O(b_exp=b_mat[ind_i_new,i],
                                 b_out=b_mat[ind_i_new,j],
                                 se_exp=se_mat[ind_i_new,i],
                                 se_out=se_mat[ind_i_new,j],
                                 n = min(n_vec[i],n_vec[j]),
                                 rho = rho_ij,
                                 random_start = random_start,t=t)
      JtoI_cML_O_res = mr_cML_O(b_exp=b_mat[ind_j_new,j],
                                     b_out=b_mat[ind_j_new,i],
                                     se_exp=se_mat[ind_j_new,j],
                                     se_out=se_mat[ind_j_new,i],
                                     n = min(n_vec[i],n_vec[j]),
                                     rho = rho_ij,
                                     random_start = random_start,t=t)

     obs_graph[i,j] = ItoJ_cML_O_res$BIC_theta
     obs_graph[j,i] = JtoI_cML_O_res$BIC_theta
     obs_graph_se[i,j] = ItoJ_cML_O_res$BIC_se
     obs_graph_se[j,i] = JtoI_cML_O_res$BIC_se
     obs_graph_pval[i,j] = ItoJ_cML_O_res$BIC_p
     obs_graph_pval[j,i] = JtoI_cML_O_res$BIC_p
      k = k + 1

    }
  }
  out = list()
  out$obs_graph = obs_graph
  out$obs_graph_se = obs_graph_se
  out$obs_graph_pval = obs_graph_pval

  return(out)

}

    set.seed(array_id)

## parallel 100 jobs, each job run 5 data perturbations
  out_list = pbmcapply::pbmclapply(1:5,function(i){    b_mat_dp = Generate_Perturb(b_mat=screen_res$b_mat,
                                se_mat=screen_res$se_mat,
                                n_vec=screen_res$n_vec,rho_mat=screen_res$rho_mat,DP_mat_list=screen_res$DP_mat_list);
                                Graph_Estimate(b_mat=b_mat_dp,
                                     se_mat=screen_res$se_mat,
                                     n_vec=screen_res$n_vec,rho_mat=screen_res$rho_mat,
                                     IJ_snp_list=screen_res$IJ_snp_list,
                                     random_start=0)},ignore.interactive=TRUE)
  obs_graph_list = lapply(out_list,function(x){x$obs_graph})
  obs_graph_se_list = lapply(out_list,function(x){x$obs_graph_se})
  obs_graph_pval_list = lapply(out_list,function(x){x$obs_graph_pval})

  out = list()
  out$obs_graph_list = obs_graph_list
  out$obs_graph_se_list = obs_graph_se_list
  out$obs_graph_pval_list = obs_graph_pval_list
  out$trait_vec = pre1_list

save(out,file=paste0('/home/panwei/lin00374/ND/result/dp_list_seed',array_id,'.RData'))
