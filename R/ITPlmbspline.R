ITPlmbspline <-
function(formula,order=2,nknots=dim(model.response(model.frame(formula)))[2],B=10000,method='residuals'){
  fisher_cf_L <- function(L){ #fisher on rows of the matrix L
    return(-2*rowSums(log(L)))
  }
  fisher_cf <- function(lambda){ #fisher on vector lambda
    return(-2*sum(log(lambda)))
  }
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
    corrected.pval <- numeric(p)
    for(var in 1:p){
      pval_var <- matrice_pval_2_2x[p,var]
      inizio <- var
      fine <- var #inizio fisso, fine aumenta salendo nelle righe
      for(riga in (p-1):1){
        fine <- fine + 1
        pval_cono <- matrice_pval_2_2x[riga,inizio:fine]
        pval_var <- max(pval_var,pval_cono)
      }
      corrected.pval[var] <- pval_var
    }
    corrected.pval <- corrected.pval[p:1]
    return(corrected.pval)
  }
  stat_lm_glob <- function(regr){
    result <- summary(regr)$f[1]
    return(result)
  }
  stat_lm_part <- function(regr){ # risultato: vettore con stat test t per ogni termine
    result <- abs(summary(regr)$coefficients[,3])
    return(result)
  }
  extract.residuals = function(regr){
    return(regr$residuals)
  }
  extract.fitted = function(regr){
    return(regr$fitted)
  }
  
  variables = all.vars(formula)
  y.name = variables[1]
  #data.all = model.frame(formula)
  cl <- match.call()
  design.matrix = model.matrix(formula)
  mf = model.frame(formula)
  data = model.response(mf)
  
  nvar = dim(design.matrix)[2] - 1
  var.names = colnames(design.matrix)
  n <- dim(data)[1]
  J <- dim(data)[2]
  
  
  print('First step: basis expansion')
  #splines coefficients:
  eval <- data
  bspl.basis <- create.bspline.basis(c(1,J),norder=order,breaks=seq(1,J,length.out=nknots))
  ascissa <- seq(1,J,1)
  
  data.fd <- Data2fd(t(data),ascissa,bspl.basis)
  coeff <- t(data.fd$coef)
  p <- dim(coeff)[2]
  
  #functional data
  npt <- 1000
  ascissa.2 <- seq(1,J,length.out=npt)
  bspl.eval.smooth <- eval.basis(ascissa.2,bspl.basis)
  data.eval <- t(bspl.eval.smooth %*% t(coeff))
  
  print('Second step: joint univariate tests')
  #univariate permutations
  #attach(mf)
  formula.const <- deparse(formula[[3]],width.cutoff = 500L) #estraggo dalla formula la parte dopo ~. se è lunga più di 500 caratteri non va!!!!
  coeffnames <- paste('coeff[,',as.character(1:p),']',sep='')
  formula.coeff <- paste(coeffnames,'~',formula.const) #ha dentro tutte le formule da usare
  formula.coeff <- sapply(formula.coeff,as.formula)
  #print(formula.coeff[[1]])
  
  formula.temp = coeff ~ design.matrix
  mf.temp = cbind(model.frame(formula.temp)[-((p+1):(p+nvar+1))],as.data.frame(design.matrix[,-1]))
  
  regr0 = lapply(formula.coeff,lm,data=mf.temp)
  
  #test statistics:
  T0_part = simplify2array(lapply(regr0,stat_lm_part))
  
  if(nvar >0){
    T0_glob <- as.numeric(lapply(regr0,stat_lm_glob))
  }else{
    method = 'responses' # se ho solo intercetta permutare i residui o le y è la stessa roba!
    T0_glob = numeric(p)
    T0_part = t(as.matrix(T0_part))
  }
  
  
  
  #calculate residuals
  if(method=='residuals'){
    #calcolo i residui. 
    #un vettore di residui per ogni variabile (1:p) 
    #e per ogni test parziale da fare (nvar+1) 
    #metto tutto in un array di dim (nvar+1,n,p)
    design.matrix.names2 = design.matrix
    var.names2 = var.names
    if(length(grep('factor',formula.const))>0){
      index.factor = grep('factor',var.names)
      replace.names = paste('group',(1:length(index.factor)),sep='')
      var.names2[index.factor] = replace.names
      colnames(design.matrix.names2) = var.names2
    }
    
    residui = array(dim=c(nvar+1,n,p))
    fitted_part = array(dim=c(nvar+1,n,p)) # per ogni variabile che devo testare ho i coeff del modello ridotto (senza quella variabile)
    formula.coeff_part = vector('list',nvar+1)
    regr0_part = vector('list',nvar+1)
    #coeff.perm_part = array(dim=c(nvar+1,n,p))
    for(ii in 2:(nvar+1)){ #il primo è l'intercetta la tratto a parte
      var.ii = var.names2[ii]
      variables.reduced = var.names2[-c(1,which(var.names2==var.ii))] #nomi di tutte le variabili tranne intercetta e quella da testare
      if(nvar>1){
        formula.temp = paste(variables.reduced,collapse=' + ')
      }else{
        formula.temp = '1' #se ho una sola variabile il modello giusto ha solo l'intercetta
      }
      
      formula.temp2 = coeff ~ design.matrix.names2
      mf.temp2 = cbind(model.frame(formula.temp2)[-((p+1):(p+nvar+1))],as.data.frame(design.matrix.names2[,-1]))
      
      
      formula.coeff.temp <- paste(coeffnames,'~',formula.temp) #ha dentro tutte le formule da usare
      formula.coeff_part[[ii]] <- sapply(formula.coeff.temp,as.formula)
      regr0_part[[ii]] = lapply(formula.coeff_part[[ii]],lm,data=mf.temp2)
      
      residui[ii,,] = simplify2array(lapply(regr0_part[[ii]],extract.residuals))
      fitted_part[ii,,] = simplify2array(lapply(regr0_part[[ii]],extract.fitted))
      
      #coeff.perm_part[ii,,] = fitted_part[ii,,] + residui[ii,,]
      
    }
    ii = 1 # intercetta
    formula.temp = paste(formula.const,' -1', sep='')
    formula.coeff.temp <- paste(coeffnames,'~',formula.temp)
    formula.coeff_part[[ii]] <- sapply(formula.coeff.temp,as.formula)
    regr0_part[[ii]] = lapply(formula.coeff_part[[ii]],lm,data=mf.temp)
    
    residui[ii,,] = simplify2array(lapply(regr0_part[[ii]],extract.residuals))
    fitted_part[ii,,] = simplify2array(lapply(regr0_part[[ii]],extract.fitted))
    
  }
  
  
  T_glob <- matrix(ncol=p,nrow=B)
  T_part = array(dim=c(B,nvar+1,p))
  #creo formula con coeff permutati
  coeffpermnames <- paste('coeff_perm[,',as.character(1:p),']',sep='')
  formula.coeff_perm <- paste(coeffpermnames,'~',formula.const) #ha dentro tutte le formule da usare
  formula.coeff_perm <- sapply(formula.coeff_perm,as.formula)
  
  for (perm in 1:B){
    # the F test is the same for both methods
    if(nvar >0){
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni,]
    }else{ # sto facendo un test sull'intercetta e basta. cambio i segni
      signs <- rbinom(n,1,0.5)*2 - 1
      coeff_perm <- coeff*signs
    }
    
    formula.temp.perm = coeff_perm ~ design.matrix
    mf.temp.perm = cbind(model.frame(formula.temp.perm)[-((p+1):(p+nvar+1))],as.data.frame(design.matrix[,-1]))
    regr_perm = lapply(formula.coeff_perm,lm,data=mf.temp.perm)
    if(nvar > 0)
      T_glob[perm,] <- as.numeric(lapply(regr_perm,stat_lm_glob))
    
    # partial tests: differ depending on the method
    if(method=='responses'){
      T_part[perm,,] = simplify2array(lapply(regr_perm,stat_lm_part))
    }else if(method=='residuals'){
      residui_perm = residui[,permutazioni,]
      regr_perm_part = vector('list',nvar+1)
      for(ii in 1:(nvar+1)){ #il primo è l'intercetta la tratto a parte
        coeff_perm = fitted_part[ii,,] + residui_perm[ii,,]  #sono le nuove y
        formula.temp.perm = coeff_perm ~ design.matrix
        mf.temp.perm = cbind(model.frame(formula.temp.perm)[-((p+1):(p+nvar+1))],as.data.frame(design.matrix[,-1]))
        regr_perm_part[[ii]] = lapply(formula.coeff_perm,lm,data=mf.temp.perm)
        T_part[perm,ii,] = simplify2array(lapply(regr_perm_part[[ii]] ,stat_lm_part))[ii,]
      }
    }
  }
  
  pval_glob <- numeric(p)
  pval_part = matrix(nrow=nvar+1,ncol=p)
  for(i in 1:p){
    pval_glob[i] <- sum(T_glob[,i]>=T0_glob[i])/B
    pval_part[,i] = colSums(T_part[,,i]>=matrix(T0_part[,i],nrow=B,ncol=nvar+1,byrow=TRUE))/B
  }
  
  #combination
  print('Third step: interval-wise combination and correction')
  q <- numeric(B)
  L_glob <- matrix(nrow=B,ncol=p)
  for(j in 1:p){
    ordine <- sort.int(T_glob[,j],index.return=T)$ix
    q[ordine] <- (B:1)/(B)
    L_glob[,j] <- q
  }
  
  L_part <- array(dim=c(B,nvar+1,p))
  for(j in 1:p){
    for(i in 1:(nvar+1)){
      ordine <- sort.int(T_part[,i,j],index.return=T)$ix
      q[ordine] <- (B:1)/(B)
      L_part[,i,j] <- q
    }
  }
  
  #asymmetric combination matrix:
  matrice_pval_asymm_glob <- matrix(nrow=p,ncol=p)
  matrice_pval_asymm_glob[p,] <- pval_glob[1:p]
  pval_2x_glob <- c(pval_glob,pval_glob)
  L_2x_glob <- cbind(L_glob,L_glob)
  
  matrice_pval_asymm_part <- array(dim=c(nvar+1,p,p))
  pval_2x_part <- cbind(pval_part,pval_part)
  L_2x_part = array(dim=c(B,nvar+1,p*2))
  for(ii in 1:(nvar+1)){
    matrice_pval_asymm_part[ii,p,] <- pval_part[ii,1:p]
    L_2x_part[,ii,] <- cbind(L_part[,ii,],L_part[,ii,])
  }
  
  for(i in (p-1):1){
    for(j in 1:p){
      inf <- j
      sup <- (p-i)+j
      T0_temp <- fisher_cf(pval_2x_glob[inf:sup])
      T_temp <- fisher_cf_L(L_2x_glob[,inf:sup])
      pval_temp <- sum(T_temp>=T0_temp)/B
      matrice_pval_asymm_glob[i,j] <- pval_temp
      for(ii in 1:(nvar + 1)){
        T0_temp <- fisher_cf(pval_2x_part[ii,inf:sup])
        T_temp <- fisher_cf_L(L_2x_part[,ii,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm_part[ii,i,j] <- pval_temp
      }
      
    }
    print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
  }
  
  #symmetric combination matrix
  matrice_pval_symm_glob <- matrix(nrow=p,ncol=4*p)
  matrice_pval_symm_part <- array(dim=c(nvar+1,p,4*p))
  
  for(i in 0:(p-1)){
    for(j in 1:(2*p)){
      matrice_pval_symm_glob[p-i,j+i+p] <- matrice_pval_asymm_glob[p-i,(j+1)%/%2]
      if(j+i>2*p-i){
        matrice_pval_symm_glob[p-i,j+i-p] <- matrice_pval_asymm_glob[p-i,(j+1)%/%2]
      }
      for(ii in 1:(nvar+1)){
        matrice_pval_symm_part[ii,p-i,j+i+p] <- matrice_pval_asymm_part[ii,p-i,(j+1)%/%2]
        if(j+i>2*p-i){
          matrice_pval_symm_part[ii,p-i,j+i-p] <- matrice_pval_asymm_part[ii,p-i,(j+1)%/%2]
        }
      }
    }
  }
  
  corrected.pval_glob <- pval.correct(matrice_pval_asymm_glob)
  corrected.pval_part = matrix(nrow=nvar+1,ncol=p)
  for(ii in 1:(nvar+1)){
    corrected.pval_part[ii,] = pval.correct(matrice_pval_asymm_part[ii,,])
  }
  
  coeff.regr = matrix(nrow=nvar+1,ncol=p)
  for(i in 1:(p)){
    coeff.regr[,i] = regr0[[i]]$coeff
  }
  coeff.t <- t(bspl.eval.smooth %*% t(coeff.regr))
  
  fitted.regr = simplify2array(lapply(regr0,extract.fitted))
  fitted.t <- t(bspl.eval.smooth %*% t(fitted.regr))
  
  rownames(corrected.pval_part) = var.names
  rownames(coeff.t) = var.names
  rownames(coeff.regr) = var.names
  rownames(pval_part) = var.names
  
  residuals.t = data.eval - fitted.t
  ybar.t = colMeans(data.eval)
  R2.t = colSums((fitted.t - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)/colSums((data.eval - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)
  
  print('Interval Testing Procedure completed')
  
  #   printresult = vector('list')
  #   printresult$call = cl
  #   #class(printresult) <- "lm"
  #   printresult$coefficients = matrix(data=apply(corrected.pval_part,1,min),ncol=1)
  #   rownames(printresult$coefficients) = var.names
  #   colnames(printresult$coefficients) = 'Minimum p-value'
  #   printresult$R2 = as.matrix(range(R2.t))
  #   colnames(printresult$R2) = 'Range of functional R-squared'
  #   rownames(printresult$R2) = c('Min R-squared', 'Max R-squared')
  #   printresult$ftest = as.matrix(min(corrected.pval_glob))
  #   colnames(printresult$ftest) = 'Minimum p-value'
  #   
  #   signif = rep('',length(var.names))
  #   signif[which(printresult$coefficients[,1] <0.001)] = '***'
  #   signif[which(printresult$coefficients[,1] <0.01 & printresult$coefficients[,1] >= 0.001)] = '**'
  #   signif[which(printresult$coefficients[,1] <0.05 & printresult$coefficients[,1] >= 0.01)] = '*'
  #   signif[which(printresult$coefficients[,1] <0.1 & printresult$coefficients[,1] >= 0.05)] = '.'
  #   printresult$coefficients[,2] = signif
  #   colnames(printresult$coefficients) = c('Minimum p-value','')
  
  #printresult
  
  
  
  ITPresult <- list(call=cl,design.matrix=design.matrix,basis='B-spline',coeff=coeff,coeff.regr=coeff.regr,
                    pval.F=pval_glob,pval.matrix.F=matrice_pval_asymm_glob,corrected.pval.F=corrected.pval_glob,
                    pval.t=pval_part,pval.matrix.t=matrice_pval_asymm_part,corrected.pval.t=corrected.pval_part,
                    data.eval=data.eval,coeff.regr.eval=coeff.t,fitted.eval=fitted.t,residuals.eval=residuals.t,
                    R2.eval=R2.t,
                    heatmap.matrix.F=matrice_pval_symm_glob,
                    heatmap.matrix.t=matrice_pval_symm_part)
  class(ITPresult) = 'ITPlm'
  return(ITPresult)
}
