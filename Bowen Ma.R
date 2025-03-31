#Startup
rm(list=ls()); cat("\f")


#Nr of MC replications
R <-1000
#Values for number of Groups 
range_g<-seq(5, 30, 5)
#Empty list for storing results
result<-list()
#Values of population parameters
beta_0<-0; beta_1<-1
#Number of bootstrap replications
B<-399
#Number of observations within a group
n_g<-60 ## change the ng to 15 and 60
alpha<-0.05
critic<-qnorm(1-0.05/2)

index<-0
for (G in range_g){
  index<-index+1
  #store results for specific G (decide what you want to store and select sufficient columns)
  result_1sim <- matrix(NA, nrow = R, ncol = 13) 
  #total sample size
  n<-n_g*G
  
  set.seed(718020247)
  for (i in 1:R){
    
    #DGP as in Cameron et al. 2007 
    #Store everything in n_g x G matrices
    z_g   <- matrix(rep(rnorm(G, 0, 1), each=n_g), nrow=n_g, ncol = G); 
    eps_g <- matrix(rep(rnorm(G, 0, 1), each=n_g), nrow=n_g, ncol = G);
    z<-matrix(rnorm(n_g*G,0,1), nrow= n_g, ncol = G); eps<- matrix(rnorm(n_g*G,0,1),nrow=n_g,ncol=G);
    Xm<-z_g+z; EPS<-eps_g+eps; Y = beta_0+beta_1*Xm+EPS
    #Some vector notation
    x<-as.vector(Xm); y<-as.vector(Y); 
    
    #Ordinary Least Squares
    X <- cbind(1,x); XX1<-solve(t(X)%*%X); beta_hat <- XX1%*%(t(X)%*%y)
    #only store estimate slope
    result_1sim[i,1]<-beta_hat[2]
    
    e_ols<-y - beta_hat[1] - beta_hat[2]*x; ses_iid<-sqrt(unname(diag((sum(e_ols^2)/(n-2))*XX1 ) ) ); 
    #only store se estimated slope
    result_1sim[i,2]<-ses_iid[2]
    result_1sim[i,3]<-(beta_hat[2]-beta_1)/ses_iid[2]
    result_1sim[i,4]<-abs((beta_hat[2]-beta_1)/ses_iid[2])>critic
    
    #Generate matrix S in Roodman et al. 2018
    S <- matrix( rep(as.vector(diag(G)),each=n_g), nrow = G, ncol = n, byrow = T)
    # u_hat:*S and sqrt(om_hat)*X
    tmp<-t(S)*e_ols; sq_om_hatX <- t(tmp)%*%X; 
    
    ses_crv<- sqrt((G/(G-1)) * unname(diag(XX1%*%(t(sq_om_hatX)%*%sq_om_hatX)%*%XX1 ) ) )
    
    result_1sim[i,5]<-ses_crv[2];
    
    t_crv<- (beta_hat[2]-beta_1)/ses_crv[2];
    result_1sim[i,6]<-t_crv;
    result_1sim[i,7]<-abs(t_crv)>critic
    ## 
    
    ########vectors for storing t-statistics for bootstrap######
    t_p_boot <- rep(NA, B)
    t_w_boot <- rep(NA, B)
    
    #Restricted residuals needed for Wild bootstrap
    b0_r <- mean(y-1*x)
    yhat_r<-b0_r+beta_1*x
    residual_r<-y-yhat_r
    
    #loop bootstrap for both pairs and wild
    for (b in 1:B) {
      #for pairs
      draw <- sample(1:G, G, replace = T) 
      X_pairs<-Xm[,draw];Y_pairs<-Y[,draw]
      x_pairs<-cbind(1,as.vector(X_pairs));y_pairs<-as.vector(Y_pairs)
      ## regression
      XX_pairs<-solve(t(x_pairs)%*%x_pairs); beta_hat_pairs <- XX_pairs%*%(t(x_pairs)%*%y_pairs)
      residuals_pairs<-as.vector(y_pairs-x_pairs%*%beta_hat_pairs)
      
      # u_hat:*S and sqrt(om_hat)*X
      tmp<-t(S)*residuals_pairs; sq_om_hatX <- t(tmp)%*%x_pairs; 
      
      ses_crv<- sqrt((G/(G-1)) * unname(diag(XX_pairs%*%(t(sq_om_hatX)%*%sq_om_hatX)%*%XX_pairs ) ) )
      
      t_crv_p<- (beta_hat_pairs[2]-beta_hat[2])/ses_crv[2];
      t_p_boot[b] <-t_crv_p 
      

      #for wild
      wild<-sample(c(-1,1), G, replace = T)
      wild_b<-rep(wild,each=n_g)
      residual_wild<-wild_b*residual_r
      
      y_wild<-yhat_r+residual_wild
      beta_hat_wild <- XX1%*%(t(X)%*%y_wild)
      residuals_wild2<-as.vector(y_wild-X%*%beta_hat_wild)
      
      tmp<-t(S)*residuals_wild2; sq_om_hatX <- t(tmp)%*%X; 
      
      ses_crv<- sqrt((G/(G-1)) * unname(diag(XX1%*%(t(sq_om_hatX)%*%sq_om_hatX)%*%XX1 ) ) )
      
      t_crv_w<- (beta_hat_wild[2]-beta_1)/ses_crv[2];
      t_w_boot[b] <-t_crv_w 
      
    }
    #sort t-statistics to determine alternative critical values
    #it is probably a good idea to store lower and upper critical value for each MC sim (not symmetric)
    t_p_boot<-sort(t_p_boot)
    t_w_boot<-sort(t_w_boot)
    
    ## pairs
    t_p_boot_low<-t_p_boot[(alpha/2)*(B+1)]
    t_p_boot_high<-t_p_boot[(1-alpha/2)*(B+1)]
    result_1sim[i,8]<-t_p_boot_low
    result_1sim[i,9]<-t_p_boot_high
    result_1sim[i,10]<-t_crv>t_p_boot_high|t_crv<t_p_boot_low
    ## Wild
    t_w_boot_low<-t_w_boot[(alpha/2)*(B+1)]
    t_w_boot_high<-t_w_boot[(1-alpha/2)*(B+1)]
    result_1sim[i,11]<-t_w_boot_low
    result_1sim[i,12]<-t_w_boot_high
    result_1sim[i,13]<-t_crv>t_w_boot_high|t_crv<t_w_boot_low
    }
    result[[index]]<-result_1sim  
}
## result store matrix result[[]]
means_list <- result
for (i in 1:6) {
  # Extract the data frame
  df <- result[[i]]
  
  # Calculate means for columns 4, 7, 10, 13
  means <- colMeans(as.matrix(df[, c(4, 7, 10, 13)]))
  
  # Store the means in the list
  means_list[[i]] <- means
}

# Print the means for each data frame
for (i in 1:6) {
  cat("Means for result[[", i, "]]:\n", sep = "")
  print(means_list[[i]])
  cat("\n")
}

