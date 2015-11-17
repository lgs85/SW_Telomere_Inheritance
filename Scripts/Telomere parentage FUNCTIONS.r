#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## DEFINE FUNCTIONS
#################################################################################*






# Model and table functions -----------------------------------------------

#Order interaction terms alphabetically (needed for function below)
orderinteraction <- function(x)
{terms <- unlist(strsplit(x,':'))
newterms <- terms[order(terms)]
paste(newterms,collapse=':')
}



#Randomise a model and get baseline relative importance
baselineimp <- function(mymodel,n.permutations = 1000)
{
  my.exp.terms <- attr(terms(formula(mymodel)),'term.labels')
  my.exp.terms <- unlist(lapply(my.exp.terms,orderinteraction))
  myresponse <- paste(formula(mymodel)[[2]])
  mydata <- mymodel$model         
  output <- mat.or.vec(n.permutations,length(my.exp.terms))
  colnames(output) <- my.exp.terms
  
  for(i in 1:n.permutations)
  {
    tempdata <- mydata
    tempdata[,myresponse] <- sample(tempdata[,myresponse])
    newmodel <- update(mymodel,data = tempdata)
    allmodels <- dredge(newmodel)
    mymodelavg <- model.avg(allmodels,subset = delta <7)
    myimp <- importance(mymodelavg)
    
    for(j in 1:length(myimp))
    {
      currentterm <- names(myimp)[j]
      currentcol <- which(colnames(output) == currentterm)
      output[i,paste(currentterm)] <- myimp[j]
    }
  }
  return(output)
}








#Get terms
getterms <- function(x,cols)
{
  temp <- which(!(is.na(x[cols])))
  return(paste(names(temp),collapse = " + "))
}






#Correct relative importance according to baseline importance
correctedimp <- function(importance,baselineimportance)
{
  bl <- baselineimportance[,colnames(baselineimportance) %in% names(importance)]
  bl <- bl[,order(colnames(bl))]
  imp <- importance[order(names(importance))]
  blmean <- apply(bl,2,mean)
  return(imp-blmean)
}



#Combine summary table from modelavereging with relative importance values
#to produce a publication-ready table
imptable <- function(modelavg,imp)
{  temp <- cbind(summary(modelavg)$coefmat.subset[,1:2],confint(modelavg))
temp <- temp[rownames(temp) != '(Intercept)',]
temp <- cbind(temp,rep(NA,nrow(temp)))
colnames(temp)[5] <- 'imp'

if(length(grep(':',names(imp)))>0)
{
  tempints <- temp[grep(':',rownames(temp)),, drop=FALSE]
  tempmains <- temp[-(grep(':',rownames(temp))),, drop = FALSE]
  impints <- imp[grep(':',names(imp))]
  impmains <- imp[-(grep(':',names(imp)))]
} else
{
  tempmains <- temp
  impmains <- imp
}



for(i in 1:length(impmains))
{
  topaste <- grep(names(impmains)[i],rownames(tempmains))
  tempmains[topaste,'imp'] <- impmains[i]
}


if(length(grep(':',names(imp)))>0)
{
  
  for(i in 1:length(impints))
  {
    impterms <- unlist(strsplit(names(impints)[i],':'))
    avgterms <- strsplit(rownames(tempints),':')
    topaste <- which(count(unlist(lapply(impterms,grep,avgterms)))[,'freq'] ==2)
    tempints[topaste,'imp'] <- impints[i]
  }
}
ifelse(exists('tempints'),
       output <- rbind(tempmains,tempints),
       output <- tempmains)

output <- output[order(output[,'imp'],decreasing = T),]


return(output)
}






# Standard error ----------------------------------------------------------

se <- function(x) sd(x,na.rm=T)/sqrt(length(x))



# Find telomere function --------------------------------------------------

findTL <- function(bird,TLfile = ParTL)
{
  if(is.na(bird))
  {
    return(NA)
  } else
  {
    if(!(bird %in% TLfile$BirdID))
    {
      return(NA)
    } else
    {
      {
        return(TLfile[TLfile$BirdID == bird,'TL'])
      }
    }
  }
}


# Knitr options -----------------------------------------------------------

chunkref <- local({
  function(chunklabel) {
    sprintf('[%s](#%s)', chunklabel, chunklabel )
  }  
})

secref <- local({
  function(seclabel) {
    sprintf('[%s](#%s)', seclabel, seclabel )
  }  
})

pgref <- local({
  function(n)
    sprintf('[Page-%i](#Page-%i)', n, n)
})

sec <- local({
  function(seclabel) {
    sprintf('# <a name="%s"/> %s', seclabel, seclabel )
  }  
})

pgcount <- local({
  pg <- 0
  function(inc=T) {
    if( inc ) { pg <<- pg + 1 }
    return( pg )
  }
})

pganchor <- local({
  function(doLabel=T) {
    if( doLabel) {
      sprintf('\n-----\nPage-%i\n<a name="Page-%i"/>\n', pgcount(inc=F), pgcount() )
    } else {
      sprintf('\n<a name="Page-%i"/>\n', pgcount() )
    }
  }
})
# 
# knit_hooks$set( anchor = function(before, options, envir) {
#   if ( before ) {
#     sprintf('<a name="%s"/>\n', options$label )
#   }
# })
# 
# knit_hooks$set( echo.label = function(before, options, envir) {
#   if ( before ) {
#     sprintf('> %s', options$label )
#   }
# })
# 
# knit_hooks$set( pgbreak = function(before, options, envir) {
#   if ( !before ) {
#     pganchor();
#   }
# })
# 






# Plotting functions ------------------------------------------------------



#Plot theme

theme_lgs <- function(addlegend=FALSE)
{
  theme_classic() +
    theme(axis.text.x = element_text(size = 14), 
          axis.text.y = element_text(size = 14), 
          axis.title.y = element_text(size = 16,vjust=0.8),
          axis.title.x = element_text(size = 16),
          legend.position = if(addlegend == TRUE) c(0.9,0.9) else "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 13))
}




# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#Percentage males in a sample
propmales <- function(x)
{
  (length(x[x=='Males'])/length(x))*100
}



#Forestplot

forestplot <- function(mytable)
{
  
  ggplot(mytable,
         aes(x = myx,
             y=Estimate)) +
    geom_point() + 
    geom_pointrange(aes(ymin=X2.5.,ymax=X97.5.)) +
    ylab('Estimate') + 
    xlab('') + 
    theme_classic() +
    theme(axis.text.x = element_text(angle=90,size = 14,hjust=1,vjust=0.5), 
          axis.text.y = element_text(size = 14), 
          axis.title.y = element_text(size = 16,vjust=0.8),
          axis.title.x = element_text(size = 16)) +
    geom_abline(slope=0,intercept=0,lty=2) +
    annotate('text',x=mytable$myx,y=min(mytable$X2.5)*1.1,label=paste('(',mytable$RI,')',sep=''))
  
}

ggsurv <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                   cens.col = 'red', lty.est = 1, lty.ci = 2,
                   cens.shape = 3, back.white = F, xlab = 'Time',
                   ylab = 'Survival', main = ''){
  
  library(ggplot2)
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = ''){
    
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
    
    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)
    
    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)
    
    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = '') {
    n <- s$strata
    
    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
    
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))
    
    col <- if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl <- if(surv.col[1] != 'gg.def'){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
    
    line <- if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
    
    pl <- pl + line
    
    pl <- if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop('Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata')
      }else if((length(surv.col) > 1 | surv.col == 'gg.def')[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
          geom_step(aes(y = low,lty = group), col = surv.col)}
    } else {pl}
    
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl <- if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                  cens.col, lty.est, lty.ci,
                                  cens.shape, back.white, xlab,
                                  ylab, main)
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}


# In-text functions -------------------------------------------------------

#Word to number
num2word <- function(x,x.units = 'none',capstart = F)
{
  if(x > 9 | x - round(x,0) != 0) 
    {
    output <- paste(x)
  } else
  {
    nums <- c(0:9)
    words <- c('zero','one','two','three','four','five','six','seven','eight','nine')
    if(capstart == T) words <- c('Zero','One','Two','Three','Four','Five','Six','Seven','Eight','Nine')
    output <- words[x == nums]
  }
    if(x.units!='none')
  {
    if(x == 1)
    {
      return(paste(output,x.units))
    } else
    {
      return(paste(output,paste0(x.units,'s')))
    }
  } else
  {
    return(output)
  }
}

#Get statistic from model
getstat <- function(model,variable,stat,standardise = T)
{
  
  
  if(class(model) %in% c('lmerMod','glmerMod'))
  {
    if(standardise == T) model <- standardize(model)
    modelterms <- rownames(summary(model)$coefficients)
    if(paste('z.',variable,sep='') %in% modelterms) variable <- paste('z.',variable,sep='')
    
    if(stat == 'CI')
    {
      est <- summary(model)$coefficients[variable,'Estimate']
      se1 <- summary(model)$coefficients[variable,'Std. Error']
      UCI <- round2(est+(1.96*se1),lessthan = F)
      LCI <- round2(est-(1.96*se1),lessthan = F)
      return(paste(LCI,', ',UCI,sep = ''))
    } else
    {
      if(stat == 'est')
      {
        return(round2(summary(model)$coefficients[variable,'Estimate'], lessthan=F))
      } else stop('you need to add this stat to the function')
    }
  } else 
  {
    if(class(model) == 'lm')
    {
      if(standardise == T) model <- standardize(model)
      modelterms <- rownames(summary(model)$coefficients)
      if(paste('z.',variable,sep='') %in% modelterms) variable <- paste('z.',variable,sep='')
      
      if(stat == 'R2')
      {
        return(round2(summary(model)$r.squared))
      } else
      {
        if(stat == 't')
        {
          return(round2(summary(model)$coef[variable,'t value']))
        } else
        {
          if(stat == 'P')
          {
            return(round2(summary(model)$coef[variable,'Pr(>|t|)']))
          } else
            {
              if(stat == 'est')
              {
                return(round2(summary(model)$coef[variable,'Estimate']))
              } else stop('you need to add this stat to the function')
            }
        }
      }
    } else
      {
        if(class(model) == 'htest')
        {
          if(stat == 'est')
          {
            return(round2(model$estimate,lessthan = F))            
          } else
          {
            if(stat == 'CI')
            {
              return(paste(round2(model$conf.int[1],lessthan = F),round2(model$conf.int[2], lessthan = F),sep = ', '))
            } else stop('you need to add this stat to this type of model')
          }

        } else stop('you need to add this type of model to the function')
      }
  }
}





#Round to 2 or 3 dp and keep trailing zeros
round2 <- function(x,lessthan = T)
{
  if(lessthan == T && x < 0.01)
  {
    return('< 0.01')
  } else
  {
    sprintf(round(x,2), fmt="%.2f")
  }
}

round3 <- function(x,lessthan = T)
{
  if(lessthan == T && x < 0.001)
  {
    return('< 0.001')
  } else
  {
    sprintf(round(x,3), fmt="%.3f")
  }
}



# function for number of observations 
give.n <- function(x)
{
  return(c(y = median(x)*0.96, label = length(x)))
}




# Find parent function ----------------------------------------------------

findpar <- function(x,parfile = pars,parent = 'mother')
{
  bird <- paste(x)
  return(parfile[which(parfile$offspring == bird),parent])
}







