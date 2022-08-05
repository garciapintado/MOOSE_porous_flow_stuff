# R script to create a table from several chunks of tabulated fluid properties
# then, extrapolate away of the input range
# the inputs are assummed to homogeneously repeat the same T sequence for every given pressure
# batch run as
# Rscript extrapolate_EOS.R
#
#library(matrixStats)
#library(raster) # just for plotting
rotate90 <- function(A) {
    # rotate matrix 90º counterclockwise
    A <- t(A)
    A[nrow(A):1,]
}
pascal2mH2O <- 1/9806.65 # [mH2O.Pa-1]
c2k <- 273.15
bulk_modulus <- 2e+09    # [Pa]

#fnames <- c('water97_tabulated_asPROST_head_chunk_from1000.csv','water97_tabulated_asPROST_main_chunk.csv')
fnames <- "water_IAPWS95.csv"               # string or vector of strings
fnameo <- "water_IAPWS95_extrap.csv"        

dT  <- 40 # increment for patched extrapolation over temperature
dP  <- 5E06
#Tmin <- 100.15
#Tmax <- 5000+c2k    # to guarantee no temperature exceed this value in my rift2ridge2D domain
                     # this range goes well beyond the expected solution but initial oscillation in the solver may hit
                     # very high values. So, the higher this value, the less likely the simulation should crash because of this
Tmin <- NULL         # do not extrapolate
Tmax <- NULL         # do not extrapolate
Pmin <- -3.0e7
#Pmax <- 18.0e8      # [Pa] ~180 km depth [to cover our standard domain up to 150 km depth
Pmax <- NULL

# extrapolation alternatives
altT <- 0             # [0:zero-derivatives, 1:slope at boundaries, 2:minimum slope at boundaries]
altP <- 3             # [as above, plus 3:bulk modulus pressure-dependence density extrapolation, and 0-slope otherwise]

alt2_fac <- 0.1       # multiplier of minimum slope for alternative 2

# --- process ---
tbl <- read.table(fnames[1], header=TRUE, sep=',')
if (length(fnames) > 1) {
    tbl <- rbind(tbl, read.table(fnames[2], header=TRUE, sep=','))
}

Tval <- unique(tbl$temperature)
Pval <- unique(tbl$pressure)
num_T <- length(Tval)
num_P <- length(Pval)
vnames <- names(tbl)[-c(1:2)]

nTb <- ceiling((min(Tval)-Tmin)/dT)
nTa <- ceiling((Tmax-max(Tval))/dT)         # values to be appended
nPb <- ceiling((min(Pval)-Pmin)/dP)         # values to be prepended towards lower values
nPa <- ceiling((Pmax-max(Pval))/dP)         # values to be appended towards higher values

# extrapolate towards lower and higher [P,T] values [all vectors in ascending order]
if (length(nTb)==0) {
    nTb <- 0
    Tvalb <- NULL
} else {
    Tvalb <- Tval[1] - (nTb:1)*dT             # values to be pre-pended
}
if (length(nTa)==0) {
    nTa <- 0
    Tvala <- NULL
} else {
    Tvala <- Tval[num_T] + (1:nTa)*dT          # values to be appended
}
if (length(nPb)==0) {
    nPb <- 0
    Pvalb <- NULL
} else {
    Pvalb <- Pval[1] - (nPb:1)*dP             # values to be pre-pended
}
if (length(nPa)==0) {
    nPa <- 0
    Pvala <- NULL
} else {
    Pvala <- Pval[num_P] + (1:nPa)*dP         # values to be appended
}

extmat <- matrix(0,nrow=(nPb+num_P+nPa)*(nTb+num_T+nTa), ncol=ncol(tbl))  # output matrix including extrapolated values
colnames(extmat) <- names(tbl)
Tval_all <- c(Tvalb,Tval,Tvala)
Pval_all <- c(Pvalb,Pval,Pvala)
extmat[,1] <- rep(Pval_all, each=nTb+num_T+nTa)
extmat[,2] <- rep(Tval_all, nPb+num_P+nPa)

for (i in 1:length(vnames)) {
    vname <- vnames[i]
    valmat <- matrix(tbl[[vname]],nrow=num_T,ncol=num_P) # [num_T,num_P]
    vrange0 <- range(tbl[[vname]])
    if (vname %in% c("cp","cv")) {
        vrange0[2] <- max(valmat[,num_P])
    }
    # image(valmat, main=vname)
    # dev.new(); plot(raster(rotate90(valmat)), main=vname, xlab = "Temperature [K]", ylab = "Pressure [Pa]")
    
    # extrapolate over temperature
    if (altT == 0) {                             # 0-derivative
        newblkb <- t(valmat[1,] + matrix(0.0, num_P,nTa))                 # [nTa,num_P]
    } else if (altT == 1 || altT == 2) {         # boundary slope || minimum boundary slope
        dVdTb <- (valmat[2,] - valmat[1,]) / diff(Tval[1:2])              # [num_P] boundary slope 
        if (altT == 2) {
            minid <- which.min(abs(dVdTb))
            dVdTb <- alt2_fac * rep(dVdTb[minid],num_P)                   # [num_P]
        }
        newblkb <- t(valmat[1,] + dVdTb %o% (Tvalb - min(Tval)))          # [nTb,num_P]
    } else {
        stop("unknown alternative for T-dependence extrapolation")
    } 

    if (altT == 0) {                     # 0-derivative
        newblka <- t(valmat[num_T,] + matrix(0.0, num_P,nTa))             # [nTa,num_P]
    } else if (altT == 1 || altT == 2) { # boundary slope || minimum boundary slope
        dVdTa <- (valmat[num_T,] - valmat[num_T-1,]) / diff(tail(Tval,2)) # [num_P] boundary slope
        if (altT == 2) {
            minid <- which.min(abs(dVdTa))
            dVdTa <- alt2_fac * rep(dVdTa[minid],num_P)                   # [num_P]
        }
        #dVdTa <- rep(median(dVdTa),num_P)
        newblka <- t(valmat[num_T,] + dVdTa %o% (Tvala - max(Tval)))      # [nTa,num_P]
    } else {
        stop("unknown alternative for T-dependence extrapolation")
    }
    #if (vname %in% c("density","enthalpy","internal_energy","viscosity","k","cp","cv")) {
    #    newblk[newblk < vrange0[1]] <- vrange0[1]
    #    newblk[newblk > vrange0[2]] <- vrange0[2]
    #}
    valmat <- rbind(newblkb, valmat, newblka) # [nTb+num_T+nTa,num_P]
    
    #extrapolate over pressure
    if (altP == 0 || vname != "density") {
        newblkb <- valmat[,1] + matrix(0.0, nTb+num_T+nTa,nPb)            # [nTb+num_T+nTa,nPb]
    } else if (altP == 1 || altP == 2) {
        dVdPb <-  (valmat[,2] - valmat[,1]) / diff(Pval[1:2])             # [nTb+num_T+nTa]
        if (altP == 2) {
            minid <- which.min(abs(dVdPb))
            dVdPb <- alt2_fac * rep(dVdPb[minid],nTb+num_T+nTa)
        }
        newblkb <- valmat[,1] + dVdPb %o% (Pvalb - min(Pval))  
    } else if (altP == 3 && vname == "density") {
        rho_0 <- valmat[,1]                                               # [nTb+num_T+nTa]
        newblkb <- rho_0 %o% exp((Pvalb - min(Pval)) / bulk_modulus)      # [nTb+num_T+nTa,nPb]
    } else {
        stop("unknown alternative for P-dependence extrapolation")
    }

    if (altP == 0 || vname != "density") {
        newblka <- valmat[,num_P] + matrix(0.0, nTb+num_T+nTa,nPa)        # [nTb+num_T+nTa,nPb]
    } else if (altP == 1 || altP == 2) {
        dVdPa <- (valmat[,num_P] - valmat[,num_P-1]) / diff(tail(Pval,2)) # [nTb+num_T+nTa]
        if (altP == 2) {
            minid <- which.min(abs(dVdPa))
            dVdPa <- alt2_fac * rep(dVdPa[minid],nTb+num_T+nTa)
        }
        newblka <- valmat[,num_P] + dVdPa %o% (Pvala - max(Pval))         # [nTb+num_T+nTa,nPa] <- [num_T+nTa] + [num_T+nTa,nPa]
    } else if (altP == 3 && vname == "density") {
        rho_0 <- valmat[,num_P]                                           # [nTb+num_T+nTa]
        newblka <- rho_0 %o% exp((Pvala - max(Pval)) / bulk_modulus)      # [nTb+num_T+nTa,nPb]
    } else {
        stop("unknown extrapolation for P-dependence extrapolation")
    }
    
    #if (vname %in% c("density","enthalpy","internal_energy","viscosity","k","cp","cv")) {
    #    newblk[newblk < vrange0[1]] <- vrange0[1]
    #    newblk[newblk > vrange0[2]] <- vrange0[2]
    #}
    valmat <- cbind(newblkb,valmat,newblka)              # [nTb+num_T+nTa,nPb+num_P+nPa]
    extmat[,i+2] <- valmat
    #stop('wait, Walter')

    # plot extrapolated tabulated property
    # dev.new(); image(valmat, main=vname)
    if (1 > 2) {
        dev.new(); plot(raster(rotate90(valmat)), main=vname, xlab = "Temperature [ºC]", ylab = "Pressure [mH2O]", xaxt="n", yaxt="n")
        nlab = 6;
        xlab_at = seq(0,1,len=length(Tval_all))
        at_ids = seq(1,length(Tval_all),by=floor(length(Tval_all)/nlab))
        axis(side=1, at=xlab_at[at_ids], labels=round(Tval_all[at_ids]-c2k,2))
        ylab_at = seq(0,1,len=length(Pval_all))
        at_ids= seq(1,length(Pval_all),by=floor(length(Pval_all)/nlab))
        axis(side=2, at=ylab_at[at_ids], labels=round(Pval_all[at_ids]*pascal2mH2O), las=1)
        
        pid <- 90
        mainlab <- paste("P=",Pval[pid],"[Pa] ~",round(Pval[pid]*pascal2mH2O,2),"[mH2O]")
        dev.new(); plot(c(Tval,Tvala),valmat[,pid], type='b', main=mainlab,
                        xlab="temperature [K]", ylab=vname)
        tid <- 15
        mainlab <- paste("T=",Tval[tid],"[K] ==",Tval[tid]-c2k,"[ºC]")
        dev.new(); plot(c(Pvalb,Pval,Pvala),valmat[tid,], type='b', main=mainlab,
                        xlab="Pressure [Pa]", ylab=vname)
    }
}
write.table(extmat, file=fnameo, row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep=",")
